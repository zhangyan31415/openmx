#include "cpu_vna_shared.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);int n=np,owners[n+1],fan[n+1],norb[n+1],*neighbors[n+1],*cells[n+1],image0[4]={0},*images[1]={image0};
 for(int a=1;a<=n;a++){owners[a]=a-1;fan[a]=1;norb[a]=1;neighbors[a]=calloc(2,4);cells[a]=calloc(2,4);neighbors[a][0]=a;neighbors[a][1]=a%n+1;}
 setenv("OPENMX_CPU_VNA_SHARED","1",1);
 for(int deferred=1;deferred>=0;deferred--){
  CpuVnaShared_Init(MPI_COMM_WORLD,n,owners,fan,neighbors,cells,images,norb,3,1,1);assert(!CpuVnaShared_WorthDeferring(0));assert(CpuVnaShared_WorthDeferring(1024));CpuVnaShared_ConfigureDeferred(deferred);
  if(CpuVnaShared_Deferred()!=deferred){if(!rank)fputs("FAIL deferred VNA policy not activated\n",stderr);CpuVnaShared_Destroy();MPI_Finalize();return 1;}
  int a=rank+1;for(int h=0;h<2;h++)for(int j=0;j<3;j++){CpuVnaShared_Row(0,a,h,0)[j]=(float)(10*a+3*h+j);if(!deferred)for(int k=1;k<4;k++)CpuVnaShared_Row(k,a,h,0)[j]=(float)(100*k+10*a+3*h+j);}
  CpuVnaShared_Publish(0);
  float*****v=calloc(4,sizeof(*v));for(int k=0;k<4;k++){v[k]=calloc(1,sizeof(**v));v[k][0]=calloc(2,sizeof(***v));for(int h=0;h<2;h++)v[k][0][h]=calloc(1,sizeof(****v));}
  if(argc>1&&deferred)CpuVnaShared_BeginWeighted(); /* must abort before derivative completion */
  int target=a%n+1;CpuVnaShared_BindAtom(v,0,target);assert(v[0][0][0][0][0]==10*target);for(int k=1;k<4;k++)assert(deferred?!v[k][0][0][0]:v[k][0][0][0]!=NULL);
  if(deferred){for(int h=0;h<2;h++)for(int j=0;j<3;j++)for(int k=1;k<4;k++)CpuVnaShared_Row(k,a,h,0)[j]=(float)(100*k+10*a+3*h+j);CpuVnaShared_DerivativesReady();}
  CpuVnaShared_BeginWeighted();for(int h=0;h<2;h++)for(int j=0;j<3;j++)for(int k=1;k<4;k++)CpuVnaShared_Row(k,a,h,0)[j]*=2;CpuVnaShared_Publish(1);
  CpuVnaShared_BindAtom(v,0,target);for(int k=0;k<4;k++)for(int h=0;h<2;h++)for(int j=0;j<3;j++)assert(v[k][0][h][0][j]==(k?2*(100*k+10*target+3*h+j):10*target+3*h+j));
  for(int k=0;k<4;k++){for(int h=0;h<2;h++)free(v[k][0][h]);free(v[k][0]);free(v[k]);}free(v);CpuVnaShared_Destroy();
 }
 for(int a=1;a<=n;a++){free(neighbors[a]);free(cells[a]);}if(!rank)puts("PASS deferred VNA raw/value/weighted transitions and legacy reinitialization");MPI_Finalize();return 0;
}
