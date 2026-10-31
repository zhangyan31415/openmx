#ifndef CPU_NSOC_HISTORY_CAPACITY_H
#define CPU_NSOC_HISTORY_CAPACITY_H
/* Include after openmx_common.h and mpi.h.
 * Compact only unused Hamiltonian-history capacity in fixed scalar RMM-DIIS.
 * The real DM/ResidualDM history length and Pulay schedule remain unchanged.
 * All consumers must retain the simple-H-mixing slot/access guard.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int CpuNSocHamiltonianHistoryCompact(void)
{
  static int option=-1;
  if(option<0){
    const char *value=getenv("OPENMX_CPU_NSOC_H_HISTORY_COMPACT");
    if(!value||!strcmp(value,"0"))option=0;
    else if(!strcmp(value,"1"))option=1;
    else {fprintf(stderr,"CPU_NSOC_HISTORY_REJECT H_HISTORY_COMPACT must be 0 or 1\n");MPI_Abort(mpi_comm_level1,96);exit(96);}
  }
  if(option && !(MD_switch==0 && Solver==2 && SpinP_switch==0 && SO_switch==0 &&
      Cnt_switch==0 && Mixing_switch==1 && Num_Mixing_pDM>0 &&
      SCF_Control_Temp==0 && List_YOUSO[39]>=2)){
    fprintf(stderr,"CPU_NSOC_HISTORY_REJECT fixed geometry scalar uncontracted NSOC Cluster RMM-DIIS with positive DM history without temperature control required\n");
    MPI_Abort(mpi_comm_level1,96);exit(96);
  }
  return option;
}
static int CpuNSocHamiltonianHistorySlots(void)
{
  return CpuNSocHamiltonianHistoryCompact()?2:List_YOUSO[39];
}
#endif
