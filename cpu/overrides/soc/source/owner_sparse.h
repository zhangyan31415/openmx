#ifndef CPU_OWNER_SPARSE_H
#define CPU_OWNER_SPARSE_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
static int CpuOwnerSparseRequested(void)
{
  const char *p=getenv("OPENMX_CPU_OWNER_SPARSE");
  if(!p||!strcmp(p,"0"))return 0;
  if(!strcmp(p,"1"))return 1;
  fprintf(stderr,"OPENMX_CPU_OWNER_SPARSE must be 0 or 1\n");
  MPI_Abort(mpi_comm_level1,95);abort();return 0;
}
#endif
