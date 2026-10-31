#ifndef CPU_DM_BALANCE_CONFIG_H
#define CPU_DM_BALANCE_CONFIG_H
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
/* Read and agree once at DFT's common entry, before feature branches. */
static int CpuDmBalanceFlag(const char *name)
{
  const char *v=getenv(name);
  if(!v||!strcmp(v,"0"))return 0;
  return !strcmp(v,"1")?1:-1;
}
static int CpuDmBalanceValidate(MPI_Comm comm,int compact,int owner,int flags[4])
{
  int j,lo[6],hi[6],local[6];
  const char *names[4]={"OPENMX_CPU_DM_SORT","OPENMX_CPU_DM_TRIM",
    "OPENMX_CPU_DM_CYCLIC","OPENMX_CPU_DM_BALANCE_PROFILE"};
  for(j=0;j<4;j++)local[j]=CpuDmBalanceFlag(names[j]);
  local[4]=compact;local[5]=owner;
  MPI_Allreduce(local,lo,6,MPI_INT,MPI_MIN,comm);
  MPI_Allreduce(local,hi,6,MPI_INT,MPI_MAX,comm);
  for(j=0;j<4;j++)if(lo[j]<0||lo[j]!=hi[j])return -1;
  if((hi[0]||hi[1]||hi[2]||hi[3])&&
    (lo[4]!=1||hi[4]!=1||lo[5]!=1||hi[5]!=1))return -2;
  if(hi[2]&&(!hi[0]||!hi[1]))return -3;
  for(j=0;j<4;j++)flags[j]=lo[j];
  return 0;
}
#endif
