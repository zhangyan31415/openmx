#ifndef CPU_XC_WORKSPACE_H
#define CPU_XC_WORKSPACE_H
#include <string.h>
#include <stdint.h>

static int CpuXcWorkspaceMode(void)
{
  static int cached=-1;
  if(cached>=0)return cached;
  const char *s=getenv("OPENMX_CPU_XC_WORKSPACE");
  int mode=!s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:!strcmp(s,"2")?2:-1;
  /* Set_XC_Grid is not a collective API. Allocation checks and reporting
   * must not introduce collectives into its existing call graph. */
  if(mode<0 || (mode && (SpinP_switch!=0||SO_switch!=0||MD_switch!=0))){
    fprintf(stderr,"CPU_XC_WORKSPACE_REJECT requires0/1/2 and fixed scalar NSOC\n");
    MPI_Abort(mpi_comm_level1,96);exit(96);
  }
  return cached=mode;
}

static void CpuXcWorkspaceCheckInputs(const double *up,const double *down)
{
  int bad=0;
  for(int i=0;i<My_NumGridD;i++){
    if(!isfinite(up[i])||!isfinite(down[i])||up[i]!=down[i] ||
       (PCC_switch && (!isfinite(PCCDensity_Grid_D[0][i]) ||
        !isfinite(PCCDensity_Grid_D[1][i]) ||
        PCCDensity_Grid_D[0][i]!=PCCDensity_Grid_D[1][i]))){bad=1;break;}
  }
  if(bad){fprintf(stderr,"CPU_XC_WORKSPACE_REJECT spin alias requires exactly equal finite spin densities/core corrections\n");MPI_Abort(mpi_comm_level1,96);exit(96);}
}

static void CpuXcWorkspaceReport(int mode,int points)
{
  static int reported=0;
  if(reported)return;
  int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
  unsigned long long local=(unsigned long long)points*8*(mode==2?3:6);
  printf("CPU_XC_WORKSPACE active=%d rank=%d old_fields=12 retained_fields=%d local_payload_bytes=%llu external_derivative_outputs_unchanged=1\n",mode,rank,mode==2?3:6,local);fflush(stdout);
  reported=1;
}
#endif
