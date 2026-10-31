#include "rep_panel.h"
#include <mkl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
typedef union {size_t bytes;long double align;void*p;} Allocation;
typedef struct {int begin,n,q[3];int*cols[3];double*u[4];} Endpoint;
typedef struct {int block,left,right,sector;double weight;int first_route,routes;} Panel;
typedef struct {int panel,dest,row,col,nr,nc;size_t offset;} Route;
typedef struct {int sector,nr,nc;int*rows,*cols;size_t offset;} Incoming;
struct RpaPlan {
 MPI_Comm comm;int rank,np,nb,pr,pc,row,col,blocks,ne,npanel,nroute,nincoming,width[3],nr[3],nc[3];
 RpaContext context;size_t limit,used,peak,entries,ns,nrecv,max_panel,max_temp;
 Endpoint*endpoint;Panel*panel;Route*route;Incoming*incoming;RpaBlock*block;
 int *rankmap,*sc,*sd,*rc,*rd;double*send,*recv,*temp,*product;
 double proof;uint64_t last_h,last_s;int last_step,s_done;
};
static void msg(char*e,size_t z,const char*s){if(e&&z)snprintf(e,z,"%s",s);}
static int agree(RpaPlan*p,int bad,char*e,size_t z,const char*s){int all;MPI_Allreduce(&bad,&all,1,MPI_INT,MPI_MAX,p->comm);if(all)msg(e,z,bad?s:"another rank rejected representative panel");return all;}
static void *take(RpaPlan*p,size_t n,size_t w){if(w&&n>(SIZE_MAX-sizeof(Allocation))/w)return NULL;size_t z=n*w+sizeof(Allocation);if(p->used>p->limit||z>p->limit-p->used)return NULL;Allocation*a=calloc(1,z);if(!a)return NULL;a->bytes=z;p->used+=z;if(p->used>p->peak)p->peak=p->used;return a+1;}
static void drop(RpaPlan*p,void*x){if(x){Allocation*a=(Allocation*)x-1;p->used-=a->bytes;free(a);}}
static int owner(int i,int b,int p){return(i/b)%p;}
static int local(int i,int b,int p){return((i/b)/p)*b+i%b;}
static int count(int n,int b,int r,int p){int whole=n/b;return(whole/p)*b+(r<whole%p?b:0)+(r==whole%p?n%b:0);}
static int same_context(RpaContext a,RpaContext b){return a.geometry==b.geometry&&a.basis==b.basis&&a.distribution==b.distribution&&a.projector_map==b.projector_map;}
static uint64_t hash_bytes(uint64_t h,const void*data,size_t bytes){const unsigned char*x=data;for(size_t i=0;i<bytes;i++){h^=x[i];h*=UINT64_C(1099511628211);}return h;}
uint64_t Rpa_ValueHash(const double*v,size_t n){const unsigned char*x=(const unsigned char*)v;uint64_t h=UINT64_C(1469598103934665603);for(size_t i=0;i<n*sizeof(double);i++){h^=x[i];h*=UINT64_C(1099511628211);}return h;}
/* Sparse all-column AO proof. The A-character distinction is stronger than
 * the density-kernel certificate and is essential for dropping A crossblocks. */
typedef struct {const RpaBasis*b;int *ao,*atom,*orb,*roff,*seen,*touch,epoch,nt;double*value;const int*perm;const double*rotation;} Proof;
static void add(Proof*a,int row,double v){if(a->seen[row]!=a->epoch){a->seen[row]=a->epoch;a->value[row]=0.;a->touch[a->nt++]=row;}a->value[row]+=v;}
static void action(Proof*a,int column){a->nt=0;if(a->epoch==INT_MAX){memset(a->seen,0,(size_t)a->b->n*4);a->epoch=0;}a->epoch++;for(int64_t z=a->b->cp[column];z<a->b->cp[column+1];z++){int row=a->b->ri[z],at=a->atom[row],orb=a->orb[row],s=a->b->types[at],w=a->b->norb[s],dest=a->perm[at];const double*d=a->rotation+a->roff[s];for(int i=0;i<w;i++)add(a,a->ao[dest]+i,d[i*w+orb]*a->b->cv[z]);}}
static double dot(Proof*a,int column){double v=0;for(int64_t z=a->b->cp[column];z<a->b->cp[column+1];z++)if(a->seen[a->b->ri[z]]==a->epoch)v+=a->b->cv[z]*a->value[a->b->ri[z]];return v;}
static void subtract(Proof*a,int column,double weight){for(int64_t z=a->b->cp[column];z<a->b->cp[column+1];z++)add(a,a->b->ri[z],-weight*a->b->cv[z]);}
static int intertwiner(RpaPlan*p,const RpaBasis*b){
 Proof a={0};a.b=b;int okay=1;double maximum=0,E[6][4],A[6][2];int detsign[6],off[5]={0};for(int k=0;k<4;k++)off[k+1]=off[k]+b->sectors[k];
 a.ao=take(p,(size_t)b->natoms+1,4);a.atom=take(p,b->n,4);a.orb=take(p,b->n,4);a.roff=take(p,b->species,4);a.seen=take(p,b->n,4);a.touch=take(p,b->n,4);a.value=take(p,b->n,8);
 if(!a.ao||!a.atom||!a.orb||!a.roff||!a.seen||!a.touch||!a.value){okay=0;goto done;}
 for(int s=1;s<b->species;s++)a.roff[s]=a.roff[s-1]+b->norb[s-1]*b->norb[s-1];
 for(int at=0;at<b->natoms;at++){int s=b->types[at];if(s<0||s>=b->species){okay=0;break;}int w=b->norb[s];if(w<1||w>128||a.ao[at]>b->n-w){okay=0;break;}a.ao[at+1]=a.ao[at]+w;for(int i=0;i<w;i++){a.atom[a.ao[at]+i]=at;a.orb[a.ao[at]+i]=i;}}
 if(!okay||a.ao[b->natoms]!=b->n){okay=0;goto done;}
 for(int g=0;g<6&&okay;g++){
  a.perm=b->permutations[g];a.rotation=b->rotations[g];if(!a.perm||!a.rotation){okay=0;break;}
  memset(a.seen,0,(size_t)b->n*4);for(int at=0;at<b->natoms;at++){int d=a.perm[at];if(d<0||d>=b->natoms||b->types[d]!=b->types[at]||a.seen[d])okay=0;else a.seen[d]=1;}memset(a.seen,0,(size_t)b->n*4);a.epoch=0;if(!okay)break;
  action(&a,0);E[g][0]=dot(&a,0);E[g][2]=dot(&a,off[1]);action(&a,off[1]);E[g][1]=dot(&a,0);E[g][3]=dot(&a,off[1]);
  for(int i=0;i<2;i++)for(int j=0;j<2;j++){double v=E[g][i]*E[g][j]+E[g][2+i]*E[g][2+j]-(i==j);if(!isfinite(v)||fabs(v)>1e-10)okay=0;}
  double determinant=E[g][0]*E[g][3]-E[g][1]*E[g][2];detsign[g]=determinant>0?1:-1;
  for(int k=2;k<4;k++){action(&a,off[k]);A[g][k-2]=dot(&a,off[k]);if(!isfinite(A[g][k-2])||fabs(fabs(A[g][k-2])-1)>1e-10)okay=0;}
  for(int k=0;k<4&&okay;k++)for(int j=0;j<b->sectors[k]&&okay;j++){
   action(&a,off[k]+j);if(k<2){subtract(&a,j,E[g][k]);subtract(&a,off[1]+j,E[g][2+k]);}else subtract(&a,off[k]+j,A[g][k-2]);
   for(int z=0;z<a.nt;z++){double v=fabs(a.value[a.touch[z]]);if(!isfinite(v)||v>1e-10)okay=0;if(v>maximum)maximum=v;}
  }
 }
 if(okay){
  int identity=0,c3=0,c2=0,plus0=1,plus1=1,minus0=1,minus1=1;
  for(int g=0;g<6;g++){double trace=E[g][0]+E[g][3];if(detsign[g]>0){if(fabs(trace-2)<1e-10)identity++;else if(fabs(trace+1)<1e-10)c3++;else okay=0;}else{if(fabs(trace)<1e-10)c2++;else okay=0;}
   if(fabs(A[g][0]-1)>1e-10)plus0=0;if(fabs(A[g][1]-1)>1e-10)plus1=0;if(fabs(A[g][0]-detsign[g])>1e-10)minus0=0;if(fabs(A[g][1]-detsign[g])>1e-10)minus1=0;
   for(int h=0;h<g;h++){double diff=0;for(int z=0;z<4;z++)diff+=fabs(E[g][z]-E[h][z]);if(diff<1e-10)okay=0;}
  }
  if(identity!=1||c3!=2||c2!=3||!((plus0&&minus1)||(plus1&&minus0)))okay=0;
  for(int g=0;g<6&&okay;g++)for(int h=0;h<6&&okay;h++){double product[4];for(int i=0;i<2;i++)for(int j=0;j<2;j++)product[2*i+j]=E[g][2*i]*E[h][j]+E[g][2*i+1]*E[h][2+j];int found=-1;for(int k=0;k<6;k++){double d=0;for(int z=0;z<4;z++)d+=fabs(product[z]-E[k][z]);if(d<1e-9)found=k;}if(found<0)okay=0;else for(int k=0;k<2;k++)if(fabs(A[g][k]*A[h][k]-A[found][k])>1e-10)okay=0;}
 }
done:p->proof=maximum;drop(p,a.ao);drop(p,a.atom);drop(p,a.orb);drop(p,a.roff);drop(p,a.seen);drop(p,a.touch);drop(p,a.value);return okay;
}
static int find_endpoint(RpaPlan*p,int begin,int rows,const RpaBasis*b,const int*rp,const int*ci,const double*uv,int*mark){
 for(int k=0;k<p->ne;k++)if(p->endpoint[k].begin==begin&&p->endpoint[k].n==rows)return k;
 int k=p->ne++;Endpoint*x=p->endpoint+k;x->begin=begin;x->n=rows;int off[5]={0};for(int s=0;s<4;s++)off[s+1]=off[s]+b->sectors[s];
 for(int s=0;s<3;s++){
  int width=p->width[s];memset(mark,0,(size_t)width*4);
  for(int i=begin;i<begin+rows;i++)for(int z=rp[i];z<rp[i+1];z++){int c=ci[z],sec=c<off[1]?0:c<off[2]?1:c<off[3]?2:3,target=sec<2?0:sec-1;if(target==s&&uv[z]!=0)mark[c-off[sec]]=1;}
  for(int j=0;j<width;j++)x->q[s]+=mark[j];x->cols[s]=take(p,x->q[s],4);if(!x->cols[s])return -1;int at=0;for(int j=0;j<width;j++)if(mark[j])x->cols[s][at++]=j;
 }
 for(int s=0;s<4;s++){int t=s<2?0:s-1;x->u[s]=take(p,(size_t)rows*x->q[t],8);if(!x->u[s])return -1;
  for(int j=0;j<x->q[t];j++)mark[x->cols[t][j]]=j;
  for(int i=0;i<rows;i++)for(int z=rp[begin+i];z<rp[begin+i+1];z++){int c=ci[z];if(c>=off[s]&&c<off[s+1]&&uv[z]!=0)x->u[s][(size_t)i*x->q[t]+mark[c-off[s]]]+=uv[z];}
 }
 return k;
}
RpaPlan *Rpa_Create(MPI_Comm comm,const RpaBasis*b,int blocks,const RpaBlock*block,int nb,int pr,int pc,const int*coords,size_t cap,char*e,size_t ez){
 RpaPlan*p=calloc(1,sizeof(*p));int bad=p==NULL,all;MPI_Allreduce(&bad,&all,1,MPI_INT,MPI_MAX,comm);if(all){free(p);msg(e,ez,"plan allocation");return NULL;}
 p->comm=MPI_COMM_NULL;MPI_Comm_dup(comm,&p->comm);MPI_Comm_rank(comm,&p->rank);MPI_Comm_size(comm,&p->np);p->used=p->peak=sizeof(*p);p->limit=cap;p->nb=nb;p->pr=pr;p->pc=pc;p->blocks=blocks;
 bad=!b||blocks<0||(blocks&&!block)||nb<1||pr<1||pc<1||(int64_t)pr*pc!=p->np||!coords||cap<sizeof(*p);
 if(!bad)bad=b->n<4||b->natoms<1||b->species<1||b->ng!=6||b->sectors[0]<1||b->sectors[0]!=b->sectors[1]||b->sectors[2]<1||b->sectors[3]<1||(int64_t)b->sectors[0]+b->sectors[1]+b->sectors[2]+b->sectors[3]!=b->n||!b->cp||!b->ri||!b->cv||!b->types||!b->norb||!b->permutations||!b->rotations||!b->context.geometry||!b->context.basis||!b->context.distribution||!b->context.projector_map;
 if(agree(p,bad,e,ez,"requires complete four-sector D3 scalar Gamma and nonzero generation context"))goto fail;
 p->context=b->context;p->width[0]=b->sectors[0];p->width[1]=b->sectors[2];p->width[2]=b->sectors[3];p->row=coords[2*p->rank];p->col=coords[2*p->rank+1];
 if(b->cp[0]!=0||b->cp[b->n]<1||b->cp[b->n]>INT_MAX)bad=1;
 for(int j=0;j<b->n&&!bad;j++){if(b->cp[j]<0||b->cp[j]>b->cp[j+1]||b->cp[j+1]>b->cp[b->n]){bad=1;break;}for(int64_t z=b->cp[j];z<b->cp[j+1];z++)if(b->ri[z]<0||b->ri[z]>=b->n||!isfinite(b->cv[z]))bad=1;}
 for(int s=0;s<b->species;s++)if(b->norb[s]<1||b->norb[s]>128)bad=1;
 for(int a=0;a<b->natoms;a++)if(b->types[a]<0||b->types[a]>=b->species)bad=1;
 for(int g=0;g<6;g++)if(!b->permutations[g]||!b->rotations[g])bad=1;
 if(agree(p,bad,e,ez,"invalid CSC basis"))goto fail;
 uint64_t metadata_hash=UINT64_C(1469598103934665603);metadata_hash=hash_bytes(metadata_hash,b->cp,(size_t)(b->n+1)*8);metadata_hash=hash_bytes(metadata_hash,b->ri,(size_t)b->cp[b->n]*4);metadata_hash=hash_bytes(metadata_hash,b->cv,(size_t)b->cp[b->n]*8);metadata_hash=hash_bytes(metadata_hash,b->types,(size_t)b->natoms*4);metadata_hash=hash_bytes(metadata_hash,b->norb,(size_t)b->species*4);metadata_hash=hash_bytes(metadata_hash,coords,(size_t)2*p->np*4);
 size_t rotation_values=0;for(int s=0;s<b->species;s++)rotation_values+=(size_t)b->norb[s]*b->norb[s];
 for(int g=0;g<6;g++){metadata_hash=hash_bytes(metadata_hash,b->permutations[g],(size_t)b->natoms*4);metadata_hash=hash_bytes(metadata_hash,b->rotations[g],rotation_values*8);}
 uint64_t shared[13]={b->context.geometry,b->context.basis,b->context.distribution,b->context.projector_map,(uint64_t)b->n,(uint64_t)b->natoms,(uint64_t)b->species,(uint64_t)nb,(uint64_t)pr,(uint64_t)pc,(uint64_t)b->sectors[0],(uint64_t)b->sectors[2],metadata_hash},low[13],high[13];
 MPI_Allreduce(shared,low,13,MPI_UINT64_T,MPI_MIN,p->comm);MPI_Allreduce(shared,high,13,MPI_UINT64_T,MPI_MAX,p->comm);for(int i=0;i<13;i++)if(low[i]!=high[i])bad=1;
 if(agree(p,bad,e,ez,"basis, actual AO metadata, context or distribution differs across ranks"))goto fail;
 /* Validate U orthogonality using row incidence, not a dense Gram matrix. */
 int *rp=take(p,(size_t)b->n+1,4),*ci=take(p,(size_t)b->cp[b->n],4),*pos=take(p,b->n,4),*mark=take(p,b->n,4);double*uv=take(p,(size_t)b->cp[b->n],8),*gram=take(p,b->n,8);
 bad=!rp||!ci||!pos||!mark||!uv||!gram;if(agree(p,bad,e,ez,"basis scratch exceeds cap")){drop(p,rp);drop(p,ci);drop(p,pos);drop(p,mark);drop(p,uv);drop(p,gram);goto fail;}
 for(int j=0;j<b->n;j++)for(int64_t z=b->cp[j];z<b->cp[j+1];z++)rp[b->ri[z]+1]++;for(int i=0;i<b->n;i++)rp[i+1]+=rp[i];memcpy(pos,rp,(size_t)b->n*4);
 for(int j=0;j<b->n;j++)for(int64_t z=b->cp[j];z<b->cp[j+1];z++){int at=pos[b->ri[z]]++;ci[at]=j;uv[at]=b->cv[z];}
 for(int j=0;j<b->n&&!bad;j++){int touched=0;for(int64_t z=b->cp[j];z<b->cp[j+1];z++)for(int t=rp[b->ri[z]];t<rp[b->ri[z]+1];t++){int k=ci[t];if(!mark[k]){mark[k]=1;pos[touched++]=k;gram[k]=0.;}gram[k]+=b->cv[z]*uv[t];}for(int t=0;t<touched;t++){int k=pos[t];if(fabs(gram[k]-(k==j))>1e-10)bad=1;mark[k]=0;}}
 if(!bad&&!intertwiner(p,b))bad=1;
 if(agree(p,bad,e,ez,"AO intertwiner, irreducible D3 E, distinct A characters or U orthogonality failed")){drop(p,rp);drop(p,ci);drop(p,pos);drop(p,mark);drop(p,uv);drop(p,gram);goto fail;}
 p->endpoint=take(p,(size_t)blocks*2,sizeof(Endpoint));p->panel=take(p,(size_t)blocks*3,sizeof(Panel));p->block=take(p,blocks,sizeof(RpaBlock));p->rankmap=take(p,p->np,4);p->sc=take(p,p->np,4);p->sd=take(p,p->np,4);p->rc=take(p,p->np,4);p->rd=take(p,p->np,4);
 bad=!p->endpoint||!p->panel||!p->block||!p->rankmap||!p->sc||!p->sd||!p->rc||!p->rd;
 if(!bad){for(int r=0;r<p->np;r++)p->rankmap[r]=-1;for(int r=0;r<p->np;r++){int a=coords[2*r],c=coords[2*r+1];if(a<0||a>=pr||c<0||c>=pc||p->rankmap[a*pc+c]!=-1)bad=1;else p->rankmap[a*pc+c]=r;}for(int s=0;s<3;s++){p->nr[s]=count(p->width[s],nb,p->row,pr);p->nc[s]=count(p->width[s],nb,p->col,pc);}}
 for(int i=0;i<blocks&&!bad;i++){
  const RpaBlock*x=block+i;size_t area=(size_t)x->rows*x->cols;if(x->row<0||x->col<0||x->rows<1||x->cols<1||x->rows>128||x->cols>128||x->row>b->n-x->rows||x->col>b->n-x->cols||(x->multiplicity!=1&&x->multiplicity!=2&&x->multiplicity!=3&&x->multiplicity!=6)||x->offset!=p->entries||area>SIZE_MAX-p->entries){bad=1;break;}p->block[i]=*x;p->entries+=area;
  int left=find_endpoint(p,x->row,x->rows,b,rp,ci,uv,mark),right=find_endpoint(p,x->col,x->cols,b,rp,ci,uv,mark);if(left<0||right<0){bad=1;break;}
  for(int s=0;s<3;s++){Endpoint*a=p->endpoint+left,*c=p->endpoint+right;if(!a->q[s]||!c->q[s])continue;if(a->q[s]>256||c->q[s]>256){bad=1;break;}p->panel[p->npanel++]=(Panel){i,left,right,s,x->multiplicity*(s==0?.5:1.),0,0};size_t sz=(size_t)a->q[s]*c->q[s],tmp=(size_t)a->n*c->q[s];if(sz>p->max_panel)p->max_panel=sz;if(tmp>p->max_temp)p->max_temp=tmp;}
 }
 drop(p,rp);drop(p,ci);drop(p,pos);drop(p,mark);drop(p,uv);drop(p,gram);
 if(agree(p,bad,e,ez,"representative block, endpoint panel or graph budget invalid"))goto fail;
 int *rows=take(p,pr,4),*cols=take(p,pc,4);bad=!rows||!cols;
 if(!bad)for(int q=0;q<p->npanel;q++){Panel*x=p->panel+q;Endpoint*a=p->endpoint+x->left,*c=p->endpoint+x->right;memset(rows,0,(size_t)pr*4);memset(cols,0,(size_t)pc*4);for(int i=0;i<a->q[x->sector];i++)rows[owner(a->cols[x->sector][i],nb,pr)]++;for(int j=0;j<c->q[x->sector];j++)cols[owner(c->cols[x->sector][j],nb,pc)]++;for(int i=0;i<pr;i++)for(int j=0;j<pc;j++)if(rows[i]&&cols[j]){if(p->nroute==INT_MAX)bad=1;else p->nroute++;}}
 p->route=take(p,p->nroute,sizeof(Route));int*mc=take(p,p->np,4),*md=take(p,p->np,4),*mrc=take(p,p->np,4),*mrd=take(p,p->np,4);bad|=!p->route||!mc||!md||!mrc||!mrd;
 if(!bad){int at=0;for(int q=0;q<p->npanel;q++){Panel*x=p->panel+q;Endpoint*a=p->endpoint+x->left,*c=p->endpoint+x->right;x->first_route=at;memset(rows,0,(size_t)pr*4);memset(cols,0,(size_t)pc*4);for(int i=0;i<a->q[x->sector];i++)rows[owner(a->cols[x->sector][i],nb,pr)]++;for(int j=0;j<c->q[x->sector];j++)cols[owner(c->cols[x->sector][j],nb,pc)]++;for(int i=0;i<pr;i++)for(int j=0;j<pc;j++)if(rows[i]&&cols[j]){int dest=p->rankmap[i*pc+j],z=rows[i]*cols[j],meta=3+rows[i]+cols[j];if(p->sc[dest]>INT_MAX-z||mc[dest]>INT_MAX-meta){bad=1;break;}p->route[at++]=(Route){q,dest,i,j,rows[i],cols[j],(size_t)p->sc[dest]};p->sc[dest]+=z;mc[dest]+=meta;x->routes++;}}}
 drop(p,rows);drop(p,cols);
 if(agree(p,bad,e,ez,"panel routes exceed budget or MPI count range")){drop(p,mc);drop(p,md);drop(p,mrc);drop(p,mrd);goto fail;}
 MPI_Alltoall(p->sc,1,MPI_INT,p->rc,1,MPI_INT,p->comm);MPI_Alltoall(mc,1,MPI_INT,mrc,1,MPI_INT,p->comm);size_t sm=0,rm=0;
 for(int r=0;r<p->np;r++){if(p->ns>INT_MAX||p->nrecv>INT_MAX||sm>INT_MAX||rm>INT_MAX)bad=1;else{p->sd[r]=(int)p->ns;p->rd[r]=(int)p->nrecv;md[r]=(int)sm;mrd[r]=(int)rm;}p->ns+=(size_t)p->sc[r];p->nrecv+=(size_t)p->rc[r];sm+=(size_t)mc[r];rm+=(size_t)mrc[r];}if(p->ns>INT_MAX||p->nrecv>INT_MAX||sm>INT_MAX||rm>INT_MAX)bad=1;
 int*sent=take(p,sm,4),*received=take(p,rm,4),*cursor=take(p,p->np,4);p->send=take(p,p->ns,8);p->recv=take(p,p->nrecv,8);p->temp=take(p,p->max_temp,8);p->product=take(p,p->max_panel,8);
 bad|=!sent||!received||!cursor||!p->send||!p->recv||!p->temp||!p->product;
 if(agree(p,bad,e,ez,"panel buffers/schedules exceed explicit cap")){drop(p,sent);drop(p,received);drop(p,cursor);drop(p,mc);drop(p,md);drop(p,mrc);drop(p,mrd);goto fail;}
 memcpy(cursor,md,(size_t)p->np*4);
 for(int q=0;q<p->nroute;q++){Route*x=p->route+q;Panel*pan=p->panel+x->panel;Endpoint*a=p->endpoint+pan->left,*c=p->endpoint+pan->right;int s=pan->sector,z=cursor[x->dest];sent[z++]=s;sent[z++]=x->nr;sent[z++]=x->nc;for(int i=0;i<a->q[s];i++)if(owner(a->cols[s][i],nb,pr)==x->row)sent[z++]=a->cols[s][i];for(int j=0;j<c->q[s];j++)if(owner(c->cols[s][j],nb,pc)==x->col)sent[z++]=c->cols[s][j];cursor[x->dest]=z;x->offset+=(size_t)p->sd[x->dest];}
 MPI_Alltoallv(sent,mc,md,MPI_INT,received,mrc,mrd,MPI_INT,p->comm);
 for(int r=0;r<p->np;r++){size_t at=(size_t)mrd[r],end=at+mrc[r];while(at<end){if(end-at<3){bad=1;break;}int s=received[at],a=received[at+1],c=received[at+2];if(s<0||s>=3||a<1||c<1||(size_t)(3+a+c)>end-at){bad=1;break;}p->nincoming++;at+=(size_t)3+a+c;}}
 p->incoming=take(p,p->nincoming,sizeof(Incoming));if(!p->incoming)bad=1;
 if(!bad){int q=0;for(int r=0;r<p->np;r++){size_t at=(size_t)mrd[r],end=at+mrc[r],value=(size_t)p->rd[r];while(at<end){Incoming*x=p->incoming+q++;x->sector=received[at++];x->nr=received[at++];x->nc=received[at++];x->offset=value;value+=(size_t)x->nr*x->nc;x->rows=take(p,x->nr,4);x->cols=take(p,x->nc,4);if(!x->rows||!x->cols){bad=1;break;}for(int i=0;i<x->nr;i++){int g=received[at++];if(g<0||g>=p->width[x->sector]||owner(g,nb,pr)!=p->row)bad=1;x->rows[i]=local(g,nb,pr);}for(int j=0;j<x->nc;j++){int g=received[at++];if(g<0||g>=p->width[x->sector]||owner(g,nb,pc)!=p->col)bad=1;x->cols[j]=local(g,nb,pc);}}if(value!=(size_t)p->rd[r]+p->rc[r])bad=1;}}
 drop(p,sent);drop(p,received);drop(p,cursor);drop(p,mc);drop(p,md);drop(p,mrc);drop(p,mrd);
 if(agree(p,bad,e,ez,"received panel route violates true BLACS owner/budget"))goto fail;
 return p;
fail:Rpa_Destroy(p);return NULL;
}
size_t Rpa_InputSize(const RpaPlan*p){return p?p->entries:0;}
size_t Rpa_LocalSize(const RpaPlan*p,int s){return p&&s>=0&&s<3?(size_t)p->nr[s]*p->nc[s]:0;}
int Rpa_Apply(RpaPlan*p,const RpaStamp*stamp,const double*v,double*const*out,const size_t*capacity,char*e,size_t ez){
 if(!p){msg(e,ez,"null panel plan");return -1;}int bad=!stamp||!out||!capacity||(p->entries&&!v);
 if(!bad){bad=!same_context(p->context,stamp->context)||(stamp->kind!=1&&stamp->kind!=2)||!stamp->revision||stamp->value_hash!=Rpa_ValueHash(v,p->entries);if(stamp->kind==1)bad|=stamp->step<1||stamp->step<=p->last_step||stamp->revision<=p->last_h;else bad|=p->s_done||stamp->step!=0;for(int s=0;s<3;s++)if(capacity[s]<Rpa_LocalSize(p,s)||(Rpa_LocalSize(p,s)&&!out[s]))bad=1;for(size_t k=0;k<p->entries;k++)if(!isfinite(v[k]))bad=1;}
 if(agree(p,bad,e,ez,"stale/nonprojected field generation, changed values, wrong context or output"))return -1;
 uint64_t field[3]={stamp->revision,(uint64_t)stamp->kind,(uint64_t)stamp->step},lo[3],hi[3];MPI_Allreduce(field,lo,3,MPI_UINT64_T,MPI_MIN,p->comm);MPI_Allreduce(field,hi,3,MPI_UINT64_T,MPI_MAX,p->comm);for(int i=0;i<3;i++)if(lo[i]!=hi[i])bad=1;
 if(agree(p,bad,e,ez,"field revision/kind/step differs across ranks"))return -1;
 for(int q=0;q<p->npanel;q++){
  Panel*x=p->panel+q;Endpoint*a=p->endpoint+x->left,*b=p->endpoint+x->right;RpaBlock*block=p->block+x->block;int s=x->sector,ka=a->q[s],kb=b->q[s],first=s==0?0:s+1,last=s==0?1:first;
  for(int partner=first;partner<=last;partner++){
   cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,a->n,kb,b->n,1.,v+block->offset,b->n,b->u[partner],kb,0.,p->temp,kb);
   cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,ka,kb,a->n,x->weight,a->u[partner],ka,p->temp,kb,partner==first?0.:1.,p->product,kb);
  }
  for(int r=x->first_route;r<x->first_route+x->routes;r++){Route*route=p->route+r;size_t at=route->offset;for(int i=0;i<ka;i++)if(owner(a->cols[s][i],p->nb,p->pr)==route->row)for(int j=0;j<kb;j++)if(owner(b->cols[s][j],p->nb,p->pc)==route->col)p->send[at++]=p->product[(size_t)i*kb+j];}
 }
 for(size_t k=0;k<p->ns;k++)if(!isfinite(p->send[k]))bad=1;
 if(agree(p,bad,e,ez,"nonfinite representative panel contraction"))return -1;
 MPI_Alltoallv(p->send,p->sc,p->sd,MPI_DOUBLE,p->recv,p->rc,p->rd,MPI_DOUBLE,p->comm);
 for(int s=0;s<3;s++)if(Rpa_LocalSize(p,s))memset(out[s],0,Rpa_LocalSize(p,s)*8);
 for(int q=0;q<p->nincoming;q++){Incoming*x=p->incoming+q;size_t at=x->offset;for(int i=0;i<x->nr;i++)for(int j=0;j<x->nc;j++)out[x->sector][(size_t)x->cols[j]*p->nr[x->sector]+x->rows[i]]+=p->recv[at++];}
 for(int s=0;s<3;s++)for(size_t k=0;k<Rpa_LocalSize(p,s);k++)if(!isfinite(out[s][k]))bad=1;
 if(agree(p,bad,e,ez,"nonfinite owner panel sum"))return -1;
 if(stamp->kind==1){p->last_step=stamp->step;p->last_h=stamp->revision;}else{p->s_done=1;p->last_s=stamp->revision;}return 0;
}
void Rpa_GetStats(const RpaPlan*p,RpaStats*s){if(!s)return;memset(s,0,sizeof(*s));if(p){s->retained_bytes=p->used;s->peak_bytes=p->peak;s->local_representatives=p->blocks;s->endpoint_count=p->ne;s->panels=p->npanel;s->routes=p->nroute;s->send_values=p->ns;s->recv_values=p->nrecv;s->max_panel_values=p->max_panel;s->intertwiner_residual=p->proof;s->value_collectives=1;}}
void Rpa_Destroy(RpaPlan*p){if(!p)return;for(int k=0;k<p->ne;k++){Endpoint*x=p->endpoint+k;for(int s=0;s<3;s++)drop(p,x->cols[s]);for(int s=0;s<4;s++)drop(p,x->u[s]);}if(p->incoming)for(int k=0;k<p->nincoming;k++){drop(p,p->incoming[k].rows);drop(p,p->incoming[k].cols);}drop(p,p->endpoint);drop(p,p->panel);drop(p,p->route);drop(p,p->incoming);drop(p,p->block);drop(p,p->rankmap);drop(p,p->sc);drop(p,p->sd);drop(p,p->rc);drop(p,p->rd);drop(p,p->send);drop(p,p->recv);drop(p,p->temp);drop(p,p->product);if(p->comm!=MPI_COMM_NULL)MPI_Comm_free(&p->comm);free(p);}
