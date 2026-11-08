/* Independent dense oracle for distributed D3 panel assembly. No production
 * contraction/routing helper is used to construct expected sector entries. */
#include "rep_panel.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define K 8
#define N (6*K)
static double U[N][N];
static int mul(int a,int b){int k=(a%3+(a/3?-1:1)*(b%3)+3)%3;return k+3*((a/3+b/3)%2);}
static double value(int i,int j){return sin(.17*(i+1)*(j+1))+cos(.23*(i+j+2));}
static double oracle(int sector,int i,int j,double scale){
 int off=sector==0?0:sector==1?4*K:5*K;double sum=0.;
 for(int a=0;a<N;a++)for(int b=0;b<N;b++){
  double u=U[a][off+i]*U[b][off+j];
  if(!sector)u=.5*(u+U[a][2*K+i]*U[b][2*K+j]);
  sum+=scale*value(a,b)*u;
 }return sum;
}
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
 const char*budget=argc>1?argv[1]:"512";setenv("OPENMX_CPU_NSOC_PANEL_BUFFER_BYTES",budget,1);
 int pr=argc>2?atoi(argv[2]):1;assert(pr>0&&np%pr==0);int pc=np/pr;
 int *coords=calloc(2*np,sizeof(int));for(int r=0;r<np;r++){coords[2*r]=(np-1-r)/pc;coords[2*r+1]=(np-1-r)%pc;}
 int64_t cp[N+1];int ri[N*N],types[N],norbs[1]={1},perm[6][N];double cv[N*N],rot[6][1];const int*perms[6];const double*rots[6];int nz=0;
 for(int c=0;c<K;c++)for(int g=0;g<6;g++){
  int a=6*c+g;types[a]=0;double t=2*acos(-1.)*(g%3)/3,co=cos(t),si=sin(t),sg=g/3?-1.:1.;
  U[a][2*c]=co/sqrt(3.);U[a][2*c+1]=-si*sg/sqrt(3.);
  U[a][2*K+2*c]=si/sqrt(3.);U[a][2*K+2*c+1]=co*sg/sqrt(3.);
  U[a][4*K+c]=1/sqrt(6.);U[a][5*K+c]=sg/sqrt(6.);
  for(int h=0;h<6;h++)perm[h][a]=6*c+mul(h,g);
 }
 for(int g=0;g<6;g++){rot[g][0]=1.;perms[g]=perm[g];rots[g]=rot[g];}
 for(int j=0;j<N;j++){cp[j]=nz;for(int i=0;i<N;i++)if(fabs(U[i][j])>1e-14){ri[nz]=i;cv[nz++]=U[i][j];}}cp[N]=nz;
 RpaBasis b={0};b.n=N;b.natoms=N;b.species=1;b.ng=6;b.sectors[0]=b.sectors[1]=2*K;b.sectors[2]=b.sectors[3]=K;b.cp=cp;b.ri=ri;b.cv=cv;b.types=types;b.norb=norbs;b.permutations=perms;b.rotations=rots;b.context=(RpaContext){1,2,3,4};
 RpaBlock*blocks=calloc(N*N,sizeof(*blocks));double*v=calloc(N*N,sizeof(double));int nb=0;
 for(int i=0;i<N;i++)for(int j=0;j<N;j++){int owner=np==1?0:1+(i*7+j*3)%(np-1);if(owner==rank){blocks[nb]=(RpaBlock){i,1,j,1,1,(size_t)nb};v[nb++]=value(i,j);}}
 char error[256];int kind;RpaPlan*p=Rpa_CreateChecked(MPI_COMM_WORLD,&b,nb,blocks,2,pr,pc,coords,64*1024*1024,error,sizeof(error),&kind);
 if(!strcmp(budget,"0")||!strcmp(budget,"invalid")){int expected=!strcmp(budget,"0")?RPA_CREATE_RESOURCE:RPA_CREATE_INVALID;int okay=!p&&kind==expected,all;MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);if(!rank)printf("%s collective budget=%s kind=%d\n",all?"PASS":"FAIL",budget,kind);Rpa_Destroy(p);free(blocks);free(v);free(coords);MPI_Finalize();return all?0:1;}
 if(!p){fprintf(stderr,"create rank%d kind%d %s\n",rank,kind,error);MPI_Abort(MPI_COMM_WORLD,2);}
 RpaStats stats;Rpa_GetStats(p,&stats);int failed=0;
 if(strtoull(budget,NULL,10)==512&&stats.value_collectives<2){if(!rank)fprintf(stderr,"FAIL bounded payload requested but only %d batch\n",stats.value_collectives);failed=1;}
 if(8*(stats.send_capacity_values+stats.recv_capacity_values)>strtoull(budget,NULL,10))failed=1;
 double*out[3];size_t cap[3];for(int s=0;s<3;s++){cap[s]=Rpa_LocalSize(p,s);out[s]=calloc(cap[s]+1,8);}
 double maxerr=0.;
 for(int step=0;step<3;step++){
  double scale=1.+step*.125;for(int z=0;z<nb;z++)v[z]=scale*value(blocks[z].row,blocks[z].col);
  RpaStamp stamp={b.context,(uint64_t)(step+1),Rpa_ValueHash(v,nb),step==1?2:1,step==1?0:step+1};
  if(Rpa_Apply(p,&stamp,v,out,cap,error,sizeof(error))){fprintf(stderr,"apply: %s\n",error);MPI_Abort(MPI_COMM_WORLD,3);}
  for(int s=0;s<3;s++){int w=s?K:2*K,lc=0,nr=0;for(int i=0;i<w;i++)if((i/2)%pr==coords[2*rank])nr++;for(int j=0;j<w;j++)if((j/2)%pc==coords[2*rank+1]){int lr=0;for(int i=0;i<w;i++)if((i/2)%pr==coords[2*rank]){double err=fabs(out[s][lc*nr+lr++]-oracle(s,i,j,scale));if(err>maxerr)maxerr=err;}lc++;}}
 }
 double global;MPI_Allreduce(&maxerr,&global,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);if(global>1e-11)failed=1;
 int all;MPI_Allreduce(&failed,&all,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 if(!rank)printf("%s ranks=%d batches=%d dense_oracle_max=%.3e peak_rank0=%zu retained_rank0=%zu\n",all?"FAIL":"PASS",np,stats.value_collectives,global,stats.peak_bytes,stats.retained_bytes);
 Rpa_Destroy(p);for(int s=0;s<3;s++)free(out[s]);free(blocks);free(v);free(coords);MPI_Finalize();return all?1:0;
}
