#ifndef OPENMX_CFS_SELECTED_RETURN_H
#define OPENMX_CFS_SELECTED_RETURN_H
#include <mpi.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* Exact value-only redistribution: Y[permutation(state),AO]=X[AO,state].
 * Column-major block-cyclic X and Y use the same validated physical grid.
 * AO tiles bound communication buffers; no full replicated coefficient matrix.
 * Output rows >=wanted are zero; leading-dimension padding is untouched. */
typedef struct {unsigned long long sent_values,received_values;size_t peak_bytes;int batches;double pack,collective,unpack;} CfsStats;
static int CfsLocalCount(int n,int b,int r,int p){int blocks=n/b,full=blocks/p,rem=blocks%p;return full*b+(r<rem?b:0)+(r==rem?n%b:0);}
static int CfsGlobal(int i,int b,int r,int p){return(int)((int64_t)(i/b)*b*p+(int64_t)r*b+i%b);}
static int CfsLocal(int i,int b,int p){return(i/b/p)*b+i%b;}
static void *CfsAlloc(MPI_Comm comm,size_t n,size_t width){if(width&&n>SIZE_MAX/width){MPI_Abort(comm,96);abort();}void*p=calloc(n?n:1,width);if(!p){MPI_Abort(comm,96);abort();}return p;}
static unsigned long long CfsHash(unsigned long long h,unsigned long long x){for(int j=0;j<8;j++){h^=(x>>(8*j))&255ULL;h*=1099511628211ULL;}return h;}
static int CfsOffsets(int n,const int*counts,int*off){long long sum=0;for(int j=0;j<n;j++){off[j]=(int)sum;sum+=counts[j];if(sum>INT_MAX)return -1;}off[n]=(int)sum;return(int)sum;}
static int CfsTransposeSelected(MPI_Comm comm,int n,int wanted,int nb,int pr,int pc,int row,int col,const int*rankmap,const int*inv,const double*x,int ldx,double*y,int ldy,int tile,CfsStats*stats){
 int p,rank;MPI_Comm_size(comm,&p);MPI_Comm_rank(comm,&rank);if(stats)memset(stats,0,sizeof(*stats));
 int valid=n>0&&wanted>0&&wanted<=n&&nb>0&&pr>0&&pc>0&&(long long)pr*pc==p&&row>=0&&row<pr&&col>=0&&col<pc&&tile>0&&rankmap&&inv&&x&&y&&x!=y;
 int nr=valid?CfsLocalCount(n,nb,row,pr):0,nc=valid?CfsLocalCount(n,nb,col,pc):0;
 if(valid&&(ldx<(nr?nr:1)||ldy<(nr?nr:1)))valid=0;
 int allvalid;MPI_Allreduce(&valid,&allvalid,1,MPI_INT,MPI_MIN,comm);if(!allvalid)return 1;
 int *seen=CfsAlloc(comm,(size_t)(n>p?n:p),sizeof(int));unsigned long long h=1469598103934665603ULL;
 h=CfsHash(h,n);h=CfsHash(h,wanted);h=CfsHash(h,nb);h=CfsHash(h,pr);h=CfsHash(h,pc);h=CfsHash(h,tile);
 for(int i=0;i<n;i++){if(inv[i]<0||inv[i]>=n){valid=0;}else if(seen[inv[i]]++)valid=0;h=CfsHash(h,(unsigned)inv[i]);}
 memset(seen,0,(size_t)p*sizeof(int));for(int i=0;i<p;i++){if(rankmap[i]<0||rankmap[i]>=p)valid=0;else if(seen[rankmap[i]]++)valid=0;h=CfsHash(h,(unsigned)rankmap[i]);}
 if(rankmap[row*pc+col]!=rank)valid=0;
 unsigned long long local_hash[2]={h,~h},common_hash[2];MPI_Allreduce(local_hash,common_hash,2,MPI_UNSIGNED_LONG_LONG,MPI_MIN,comm);
 if(common_hash[0]!=h||common_hash[1]!=~h)valid=0;
 MPI_Allreduce(&valid,&allvalid,1,MPI_INT,MPI_MIN,comm);free(seen);if(!allvalid)return 2;
 int *sc=CfsAlloc(comm,pr,sizeof(int)),*so=CfsAlloc(comm,pr+1,sizeof(int)),*rc=CfsAlloc(comm,pc,sizeof(int)),*ro=CfsAlloc(comm,pc+1,sizeof(int));
 int *sj=CfsAlloc(comm,nc,sizeof(int));int kept_rows=CfsLocalCount(wanted,nb,row,pr);int *rr=CfsAlloc(comm,kept_rows,sizeof(int));
 int *cursor=CfsAlloc(comm,p,sizeof(int));
 for(int j=0;j<nc;j++){int q=CfsGlobal(j,nb,col,pc);if(inv[q]<wanted)sc[(inv[q]/nb)%pr]++;}
 for(int q=0;q<n;q++)if(inv[q]<wanted&&(inv[q]/nb)%pr==row)rc[(q/nb)%pc]++;
 int sk=CfsOffsets(pr,sc,so),rk=CfsOffsets(pc,rc,ro);if(sk<0||rk!=kept_rows){MPI_Abort(comm,96);abort();}
 memcpy(cursor,so,pr*sizeof(int));for(int j=0;j<nc;j++){int q=CfsGlobal(j,nb,col,pc);if(inv[q]<wanted)sj[cursor[(inv[q]/nb)%pr]++]=j;}
 memcpy(cursor,ro,pc*sizeof(int));for(int q=0;q<n;q++)if(inv[q]<wanted&&(inv[q]/nb)%pr==row)rr[cursor[(q/nb)%pc]++]=CfsLocal(inv[q],nb,pr);
 int maxao=tile<n?tile:n;
 int *asc=CfsAlloc(comm,pc,sizeof(int)),*aso=CfsAlloc(comm,pc+1,sizeof(int)),*arc=CfsAlloc(comm,pr,sizeof(int)),*aro=CfsAlloc(comm,pr+1,sizeof(int));
 int *asi=CfsAlloc(comm,maxao,sizeof(int)),*ari=CfsAlloc(comm,maxao,sizeof(int));
 int *sendcounts=CfsAlloc(comm,p,sizeof(int)),*recvcounts=CfsAlloc(comm,p,sizeof(int)),*sd=CfsAlloc(comm,p+1,sizeof(int)),*rd=CfsAlloc(comm,p+1,sizeof(int));
 double *send=NULL,*recv=NULL;size_t sendcap=0,recvcap=0;
 size_t meta=((size_t)4*(pr+pc)+4+nc+kept_rows+2*(size_t)maxao+5*(size_t)p+2)*sizeof(int);
 unsigned long long transmitted=0,received=0;int bad=0;
 for(int j=0;j<nc;j++)for(int i=0;i<nr;i++)y[(size_t)j*ldy+i]=0.;
 for(int begin=0;begin<n;){int end=n-begin<tile?n:begin+tile;double t0=MPI_Wtime();
  memset(asc,0,pc*sizeof(int));memset(arc,0,pr*sizeof(int));
  for(int ao=begin;ao<end;ao++){int r=(ao/nb)%pr,c=(ao/nb)%pc;if(r==row)asc[c]++;if(c==col)arc[r]++;}
  CfsOffsets(pc,asc,aso);CfsOffsets(pr,arc,aro);
  memcpy(cursor,aso,pc*sizeof(int));for(int ao=begin;ao<end;ao++)if((ao/nb)%pr==row)asi[cursor[(ao/nb)%pc]++]=CfsLocal(ao,nb,pr);
  memcpy(cursor,aro,pr*sizeof(int));for(int ao=begin;ao<end;ao++)if((ao/nb)%pc==col)ari[cursor[(ao/nb)%pr]++]=CfsLocal(ao,nb,pc);
  for(int r=0;r<pr;r++)for(int c=0;c<pc;c++){
   int peer=rankmap[r*pc+c];long long a=(long long)sc[r]*asc[c],b=(long long)rc[c]*arc[r];
   if(a>INT_MAX||b>INT_MAX){MPI_Abort(comm,96);abort();}sendcounts[peer]=(int)a;recvcounts[peer]=(int)b;
  }
  int ns=CfsOffsets(p,sendcounts,sd),nt=CfsOffsets(p,recvcounts,rd);if(ns<0||nt<0){MPI_Abort(comm,96);abort();}
  if(!send||ns>sendcap){free(send);send=CfsAlloc(comm,ns,sizeof(double));sendcap=ns;}
  if(!recv||nt>recvcap){free(recv);recv=CfsAlloc(comm,nt,sizeof(double));recvcap=nt;}
  for(int r=0;r<pr;r++)for(int c=0;c<pc;c++){
   int peer=rankmap[r*pc+c],at=sd[peer];
   for(int a=so[r];a<so[r+1];a++){const double*src=x+(size_t)sj[a]*ldx;for(int b=aso[c];b<aso[c+1];b++)send[at++]=src[asi[b]];}
   if(at!=sd[peer]+sendcounts[peer])bad=1;
  }
  double t1=MPI_Wtime();MPI_Alltoallv(send,sendcounts,sd,MPI_DOUBLE,recv,recvcounts,rd,MPI_DOUBLE,comm);double t2=MPI_Wtime();
  for(int r=0;r<pr;r++)for(int c=0;c<pc;c++){
   int peer=rankmap[r*pc+c],at=rd[peer];
   for(int a=ro[c];a<ro[c+1];a++){int destrow=rr[a];for(int b=aro[r];b<aro[r+1];b++){double v=recv[at++];if(!isfinite(v))bad=1;y[(size_t)ari[b]*ldy+destrow]=v;}}
   if(at!=rd[peer]+recvcounts[peer])bad=1;
  }
  if(stats){stats->batches++;stats->pack+=t1-t0;stats->collective+=t2-t1;stats->unpack+=MPI_Wtime()-t2;size_t bytes=meta+(sendcap+recvcap)*sizeof(double);if(bytes>stats->peak_bytes)stats->peak_bytes=bytes;}
  transmitted+=(unsigned)ns;received+=(unsigned)nt;begin=end;
 }
 if(stats){stats->sent_values=transmitted;stats->received_values=received;}
 int anybad;MPI_Allreduce(&bad,&anybad,1,MPI_INT,MPI_MAX,comm);
 free(sc);free(so);free(rc);free(ro);free(sj);free(rr);free(cursor);free(asc);free(aso);free(arc);free(aro);free(asi);free(ari);free(sendcounts);free(recvcounts);free(sd);free(rd);free(send);free(recv);
 return anybad?3:0;
}
#endif
