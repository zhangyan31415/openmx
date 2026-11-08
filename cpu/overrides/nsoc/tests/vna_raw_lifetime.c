#define _DEFAULT_SOURCE
#include "cpu_vna_shared.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#ifdef __linux__
#include <sys/mman.h>
#ifdef CPU_TEST_WRAP_MADVISE
int __wrap_madvise(void *base,size_t bytes,int advice){(void)base;(void)bytes;(void)advice;errno=ENOSYS;return -1;}
#endif
static size_t resident(void *base,size_t bytes){
 size_t page=(size_t)sysconf(_SC_PAGESIZE);uintptr_t lo=((uintptr_t)base+page-1)/page*page,hi=((uintptr_t)base+bytes)/page*page;
 if(hi<=lo)return 0;size_t pages=(hi-lo)/page;unsigned char *vec=malloc(pages);assert(vec);assert(!mincore((void*)lo,hi-lo,vec));size_t n=0;for(size_t i=0;i<pages;i++)if(vec[i]&1)n+=page;free(vec);return n;
}
#else
static size_t resident(void *base,size_t bytes){(void)base;(void)bytes;return 0;}
#endif
int main(int argc,char **argv){
 MPI_Init(&argc,&argv);int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);assert(np>=2);
 int n=np,p=65536,owners[n+1],fan[n+1],norb[n+1],*neighbor[n+1],*cell[n+1],image0[4]={0},*images[1]={image0};
 for(int a=1;a<=n;a++){owners[a]=a-1;fan[a]=n-1;norb[a]=1;neighbor[a]=calloc(n,sizeof(int));cell[a]=calloc(n,sizeof(int));for(int h=0;h<n;h++)neighbor[a][h]=(a+h-1)%n+1;}
 setenv("OPENMX_CPU_VNA_SHARED","1",1);CpuVnaShared_Init(MPI_COMM_WORLD,n,owners,fan,neighbor,cell,images,norb,p,n-1,1);CpuVnaShared_ConfigureDeferred(1);
 int a=rank+1;for(int h=0;h<n;h++)for(int j=0;j<p;j++)CpuVnaShared_Row(0,a,h,0)[j]=(float)(100*a+h+j%13);
 CpuVnaShared_Publish(0);size_t bytes=CpuVnaShared_OwnedBytes(),before=bytes?resident(CpuVnaShared_Row(0,1,0,0),bytes):0;
 size_t removed=CpuVnaShared_DropRaw0();assert(CpuVnaShared_Raw0Dropped());
 size_t after=bytes?resident(CpuVnaShared_Row(0,1,0,0),bytes):0;
 if(bytes&&removed){assert(before>0);assert(after<before/4);printf("PASS shared backing release before=%zu after=%zu advised=%zu\n",before,after,removed);}
#ifdef CPU_TEST_WRAP_MADVISE
 assert(!removed);if(bytes)assert(after>=before);assert(CpuVnaShared_Row(0,a,0,0)[0]==100*a);
#endif
 if(argc>1)CpuVnaShared_DerivativesReady(); /* must reject incomplete raw restoration */
 for(int h=0;h<n;h++)for(int j=0;j<p;j++)for(int k=0;k<4;k++)CpuVnaShared_Row(k,a,h,0)[j]=(float)(1000*k+100*a+h+j%13);
 CpuVnaShared_RestoreRaw0();assert(!CpuVnaShared_Raw0Dropped());CpuVnaShared_DerivativesReady();
 for(int b=1;b<=n;b++)for(int h=0;h<n;h++)for(int j=0;j<p;j++)assert(CpuVnaShared_Row(0,b,h,0)[j]==(float)(100*b+h+j%13));
 CpuVnaShared_BeginWeighted();for(int h=0;h<n;h++)for(int j=0;j<p;j++)for(int k=1;k<4;k++)CpuVnaShared_Row(k,a,h,0)[j]*=2;CpuVnaShared_Publish(1);
 for(int b=1;b<=n;b++)for(int h=0;h<n;h++)for(int j=0;j<p;j++)for(int k=0;k<4;k++)assert(CpuVnaShared_Row(k,b,h,0)[j]==(float)(k?2*(1000*k+100*b+h+j%13):100*b+h+j%13));
 CpuVnaShared_Destroy();for(int b=1;b<=n;b++){free(neighbor[b]);free(cell[b]);}if(!rank)puts("PASS raw0 drop/restore, needed-node values and weighted epoch");MPI_Finalize();return 0;
}
