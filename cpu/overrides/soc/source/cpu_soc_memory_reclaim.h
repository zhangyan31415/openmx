#ifndef CPU_SOC_MEMORY_RECLAIM_H
#define CPU_SOC_MEMORY_RECLAIM_H
#include <malloc.h>
#include "mkl_service.h"
/* Standalone legacy SOC builds do not link the optional unified experiment. */
extern int CpuSOCMemoryEligible(void) __attribute__((weak));
static void CpuSOCMemoryReclaim(const char *stage)
{
  if(!CpuSOCMemoryEligible||!CpuSOCMemoryEligible())return;
  int rank,buffers=0;
  unsigned long long local[3],total[3];
  local[0]=(unsigned long long)MKL_Mem_Stat(&buffers);
  MKL_Free_Buffers();
  local[1]=(unsigned long long)MKL_Mem_Stat(&buffers);
  local[2]=(unsigned long long)(malloc_trim(0)!=0);
  MPI_Reduce(local,total,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,Host_ID,mpi_comm_level1);
  MPI_Comm_rank(mpi_comm_level1,&rank);
  if(rank==Host_ID){
    printf("CPU_SOC_MEMORY_RECLAIM stage=%s logical_MKL_before_bytes=%llu logical_MKL_after_bytes=%llu heap_trimmed_ranks=%llu active_matrices_retained=1\n",stage,total[0],total[1],total[2]);fflush(stdout);
  }
}
#endif
