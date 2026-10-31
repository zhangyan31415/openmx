#ifndef CPU_NSOC_SYMHISTORY_GUARD_H
#define CPU_NSOC_SYMHISTORY_GUARD_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
extern int flag_export_DM;
static int CpuNSocSymHistoryEnabled(void)
{
 static int enabled=-1;
 if(enabled<0){const char *v=getenv("OPENMX_CPU_NSOC_SYM_HISTORY");
  if(!v||!strcmp(v,"0"))enabled=0;
  else if(!strcmp(v,"1"))enabled=1;
  else enabled=2;
  int low,high;
  MPI_Allreduce(&enabled,&low,1,MPI_INT,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(&enabled,&high,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(low!=high||high<0||high>1){fprintf(stderr,"CPU_NSOC_SYMHISTORY_REJECT switch must be agreed0or1 acrossMPIranks\n");MPI_Abort(mpi_comm_level1,96);exit(96);}
 }
 static int validated=0;
 if(enabled&&!validated){
  const char *project=getenv("OPENMX_CPU_NSOC_FULL_PROJECT"),*realspace=getenv("OPENMX_CPU_REALSPACE_D3");
  int valid=(MD_switch==0&&Solver==2&&SpinP_switch==0&&SO_switch==0&&Cnt_switch==0&&Mixing_switch==1&&
      Num_Mixing_pDM==30&&List_YOUSO[16]>=32&&SCF_Control_Temp==0&&Scf_RestartFromFile==0&&
      LNO_flag==0&&cal_partial_charge==0&&xanes_calc==0&&CWF_Calc==0&&flag_export_DM==0&&
      MO_fileout==0&&CDDF_on==0&&LNAO_calc_flag==0&&LNBO_calc_flag==0&&NBO_switch==0&&
      Dos_fileout==0&&DosGauss_fileout==0&&Hub_U_switch==0&&core_hole_state_flag==0&&
      project&&!strcmp(project,"1")&&realspace&&(!strcmp(realspace,"scf")||!strcmp(realspace,"1")));
  int allvalid;MPI_Allreduce(&valid,&allvalid,1,MPI_INT,MPI_MIN,mpi_comm_level1);
  if(!allvalid){
   fprintf(stderr,"CPU_NSOC_SYMHISTORY_REJECT fixed fresh projected realspace scalar NSOC Cluster RMM-DIIS History30 required\n");MPI_Abort(mpi_comm_level1,96);exit(96);
  }
  validated=1;
 }
 return enabled;
}
#endif
