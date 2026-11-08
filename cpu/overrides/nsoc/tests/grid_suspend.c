#define _GNU_SOURCE
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sched.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sys/mman.h>
#include "openmx_common.h"
#include "cpu_grid_shared.h"

static size_t interior(void *pointer,size_t bytes,uintptr_t *first)
{
  uintptr_t page=(uintptr_t)sysconf(_SC_PAGESIZE),p=(uintptr_t)pointer;
  *first=((p+page-1)/page)*page;
  uintptr_t end=((p+bytes)/page)*page;
  return end>*first?(size_t)(end-*first):0;
}
static size_t resident(void *pointer,size_t bytes)
{
  uintptr_t first;size_t length=interior(pointer,bytes,&first);
  size_t pages=length/(size_t)sysconf(_SC_PAGESIZE),sum=0;
  unsigned char *vec=malloc(pages?pages:1);assert(vec);
  assert(mincore((void*)first,length,vec)==0);
  for(size_t k=0;k<pages;k++)sum+=(vec[k]&1)!=0;
  free(vec);return sum;
}
static int zero_interior(void *pointer,size_t bytes)
{
  uintptr_t first;size_t length=interior(pointer,bytes,&first);
  volatile unsigned char *data=(volatile unsigned char*)first;
  for(size_t k=0;k<length;k++)if(data[k])return 0;
  return 1;
}
static float pattern(int epoch,int atom,int grid,int orbital)
{
  return (float)(epoch*100000+atom*1000+grid%997+orbital);
}
int main(int argc,char **argv)
{
  MPI_Init(&argc,&argv);int rank,np,local_rank;MPI_Comm local;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);
  MPI_Comm_rank(local,&local_rank);assert(argc>1);
  char *list=strdup(argv[1]),*part=strtok(list,",");
  for(int k=0;k<local_rank&&part;k++)part=strtok(NULL,",");assert(part);
  cpu_set_t affinity;CPU_ZERO(&affinity);CPU_SET(atoi(part),&affinity);
  assert(sched_setaffinity(0,sizeof(affinity),&affinity)==0);free(list);
  unsigned cpu,node;assert(syscall(SYS_getcpu,&cpu,&node,NULL)==0);
  unsigned minimum,maximum;MPI_Allreduce(&node,&minimum,1,MPI_UNSIGNED,MPI_MIN,local);
  MPI_Allreduce(&node,&maximum,1,MPI_UNSIGNED,MPI_MAX,local);assert(minimum!=maximum);
  int n=4*np,*owners=calloc((size_t)n+1,sizeof(int));
  int *grids=calloc((size_t)n+1,sizeof(int)),*species=calloc((size_t)n+1,sizeof(int));
  int norbs[2]={2,3};FNAN=calloc((size_t)n+1,sizeof(int));natn=calloc((size_t)n+1,sizeof(int*));
  Type_Orbs_Grid **canonical=calloc((size_t)n+1,sizeof(*canonical));
  Type_Orbs_Grid **read=calloc((size_t)n+1,sizeof(*read));assert(owners&&grids&&species&&FNAN&&natn&&canonical&&read);
  for(int a=1;a<=n;a++){
    owners[a]=(a-1)/4;grids[a]=16384+17*a;species[a]=a%2;FNAN[a]=n-1;
    natn[a]=malloc((size_t)n*sizeof(int));assert(natn[a]);natn[a][0]=a;
    int h=1;for(int b=1;b<=n;b++)if(b!=a)natn[a][h++]=b;
  }
  setenv("OPENMX_CPU_GRID_SHARED","1",1);setenv("OPENMX_CPU_GRID_NUMA","1",1);
  mpi_comm_level1=MPI_COMM_WORLD;
  CpuGridShared_Init(MPI_COMM_WORLD,n,owners,grids,species,norbs);
  int bad=0;size_t observed_before=0,observed_after=0;
  for(int epoch=1;epoch<=3;epoch++){
    if(epoch>1){CpuGridShared_BeginRebuild();if(CpuGridShared_IsSuspended())bad=1;}
    for(int a=1;a<=n;a++)if(owners[a]==rank)for(int g=0;g<grids[a];g++){
      Type_Orbs_Grid *row=CpuGridShared_Row(a,g);if(!row)MPI_Abort(MPI_COMM_WORLD,2);
      for(int o=0;o<norbs[species[a]];o++)row[o]=pattern(epoch,a,g,o);
    }
    if(epoch>1)for(int a=1;a<=n;a++)if(CpuGridShared_ReadRow(a,0))bad=1;
    CpuGridShared_Publish();
    for(int a=1;a<=n;a++){
      Type_Orbs_Grid *owner=CpuGridShared_Row(a,0),*value=CpuGridShared_ReadRow(a,0);
      if(epoch==1){canonical[a]=owner;read[a]=value;}
      if(owner!=canonical[a]||value!=read[a]||!value)bad=1;
      if(value)for(int g=0;g<grids[a];g++)for(int o=0;o<norbs[species[a]];o++)
        if(value[(size_t)g*norbs[species[a]]+o]!=pattern(epoch,a,g,o))bad=1;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if(epoch==3)break;
    for(int a=1;a<=n;a++){
      size_t bytes=(size_t)grids[a]*norbs[species[a]]*sizeof(Type_Orbs_Grid);
      observed_before+=resident(canonical[a],bytes);
      if(read[a]!=canonical[a])observed_before+=resident(read[a],bytes);
    }
    int removed=CpuGridShared_Suspend();
#if defined(CPU_GRID_TEST_REMOVE_UNSUPPORTED)
    if(removed)bad=1;
#else
    if(!removed)bad=1;
#endif
    if(!CpuGridShared_IsSuspended()||CpuGridShared_Suspend()!=removed)bad=1;
    /* Collect residency before reading any zero payload. */
    size_t remaining=0;
    for(int a=1;a<=n;a++){
      if(CpuGridShared_ReadRow(a,0))bad=1;
      size_t bytes=(size_t)grids[a]*norbs[species[a]]*sizeof(Type_Orbs_Grid);
      remaining+=resident(canonical[a],bytes);
      if(read[a]!=canonical[a])remaining+=resident(read[a],bytes);
    }
    observed_after+=remaining;
#if !defined(CPU_GRID_TEST_REMOVE_UNSUPPORTED)
    if(remaining)bad=1;
    for(int a=1;a<=n;a++){
      size_t bytes=(size_t)grids[a]*norbs[species[a]]*sizeof(Type_Orbs_Grid);
      if(!zero_interior(canonical[a],bytes)||!zero_interior(read[a],bytes))bad=1;
    }
#endif
    MPI_Barrier(MPI_COMM_WORLD);
  }
  CpuGridShared_Destroy();if(CpuGridShared_IsSuspended())bad=1;
  for(int a=1;a<=n;a++)free(natn[a]);free(natn);free(FNAN);
  free(owners);free(grids);free(species);free(canonical);free(read);
  int any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  if(!rank)printf("%s grid_suspend ranks=%d NUMA=%u..%u resident_pages_before=%zu after=%zu pointer_stability=1 native_values_after_rebuild=1\n",any?"FAIL":"PASS",np,minimum,maximum,observed_before,observed_after);
  MPI_Comm_free(&local);MPI_Finalize();return any?1:0;
}
