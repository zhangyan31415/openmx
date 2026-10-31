#ifndef OPENMX_CPU_WORKSPACE_BOUNDS_H
#define OPENMX_CPU_WORKSPACE_BOUNDS_H
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

typedef struct { int l, mul, mag, enabled; } CpuWorkspaceBounds;

/* Only PAO scratch dimensions are planned here. The actual species basis,
   integration loops, projector dimensions and global List_YOUSO stay intact. */
static int CpuWorkspaceBounds_Plan(int nspecies,const int *maxl,
                                   int *const *counts,int legacy_l,int legacy_mul,
                                   CpuWorkspaceBounds *out)
{
  int s,l,largest_l=0,largest_mul=0;
  if(!out || nspecies<=0 || !maxl || !counts || legacy_l<=0 || legacy_mul<=0) return 0;
  for(s=0;s<nspecies;s++){
    if(!counts[s] || maxl[s]<0 || maxl[s]>=legacy_l || maxl[s]>(INT_MAX-1)/2) return 0;
    if(largest_l<maxl[s]) largest_l=maxl[s];
    for(l=0;l<=maxl[s];l++){
      if(counts[s][l]<0 || counts[s][l]>legacy_mul) return 0;
      if(largest_mul<counts[s][l]) largest_mul=counts[s][l];
    }
  }
  if(largest_mul==0) return 0;
  out->l=largest_l+1;
  out->mul=largest_mul;
  out->mag=2*largest_l+1;
  out->enabled=1;
  return 1;
}

static CpuWorkspaceBounds CpuWorkspaceBounds_Select(const char *label,int *reported,
                     int nspecies,const int *maxl,int *const *counts,
                     int legacy_l,int legacy_mul,MPI_Comm comm)
{
  CpuWorkspaceBounds out;
  const char *env=getenv("OPENMX_CPU_TIGHT_WORKSPACE");
  int rank=0;
  MPI_Comm_rank(comm,&rank);
  if(env && strcmp(env,"0") && strcmp(env,"1")){
    if(rank==0) fprintf(stderr,"[CPU_TIGHT_WORKSPACE] OPENMX_CPU_TIGHT_WORKSPACE must be 0 or 1\n");
    MPI_Abort(comm,1); abort();
  }
  out.l=legacy_l; out.mul=legacy_mul; out.mag=2*legacy_l+1; out.enabled=0;
  if(!env || !strcmp(env,"0")) return out;
  if(!CpuWorkspaceBounds_Plan(nspecies,maxl,counts,legacy_l,legacy_mul,&out)){
    if(rank==0) fprintf(stderr,"[CPU_TIGHT_WORKSPACE] %s invalid selected-basis scratch bounds\n",label);
    MPI_Abort(comm,1); abort();
  }
  if(!*reported){
    if(rank==0){
      printf("[CPU_TIGHT_WORKSPACE] function=%s legacy_L_mul_m=%d,%d,%d selected_L_mul_m=%d,%d,%d basis_unchanged=1 projector_unchanged=1\n",
             label,legacy_l,legacy_mul,2*legacy_l+1,out.l,out.mul,out.mag);
      fflush(stdout);
    }
    *reported=1;
  }
  return out;
}
#endif
