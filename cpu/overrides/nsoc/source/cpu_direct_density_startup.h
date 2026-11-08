#ifndef CPU_DIRECT_DENSITY_STARTUP_H
#define CPU_DIRECT_DENSITY_STARTUP_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
static int CpuNSocDirectDensityStartup(void)
{
 static int selected=-1;
 if(selected<0){
  const char *v=getenv("OPENMX_CPU_NSOC_STARTUP");
  selected=!v||!strcmp(v,"0")?0:!strcmp(v,"1")?1:2;
  int lo,hi;MPI_Allreduce(&selected,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(&selected,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(lo!=hi||hi>1){fprintf(stderr,"CPU_NSOC_STARTUP requires agreed0or1\n");MPI_Abort(mpi_comm_level1,96);exit(96);}
 }
 if(selected){
  const char *history=getenv("OPENMX_CPU_NSOC_SYM_HISTORY");
  const char *project=getenv("OPENMX_CPU_NSOC_FULL_PROJECT");
  if(MD_switch||Solver!=2||SpinP_switch||SO_switch||Cnt_switch||LNO_flag||
     Scf_RestartFromFile||SCF_Control_Temp||Mixing_switch!=1||Num_Mixing_pDM!=30||
     !history||strcmp(history,"1")||!project||strcmp(project,"1")){
   fprintf(stderr,"CPU_NSOC_STARTUP requires fixed fresh projected scalar Cluster RMM-DIIS History30\n");
   MPI_Abort(mpi_comm_level1,96);exit(96);
  }
 }
 return selected;
}
#endif
