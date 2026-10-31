#ifndef CPU_NSOC_GHOST_H
#define CPU_NSOC_GHOST_H
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
static int CpuNSOCGhostFlag(const char*name){const char*s=getenv(name);return !s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:-1;}
static int CpuNSOCGhostAgree(MPI_Comm comm,int valid,int*hs,int*olp){
 int local[3]={CpuNSOCGhostFlag("OPENMX_CPU_LOCAL_HS"),CpuNSOCGhostFlag("OPENMX_CPU_LOCAL_OLP"),!!valid},lo[3],hi[3];
 MPI_Allreduce(local,lo,3,MPI_INT,MPI_MIN,comm);MPI_Allreduce(local,hi,3,MPI_INT,MPI_MAX,comm);
 if(lo[0]<0||lo[1]<0||lo[0]!=hi[0]||lo[1]!=hi[1]||((hi[0]||hi[1])&&!lo[2])||(hi[1]&&!hi[0]))return -1;
 *hs=lo[0];*olp=lo[1];return 0;
}
#endif
