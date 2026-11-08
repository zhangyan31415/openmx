#ifndef CPU_RHOK_SCOPE_H
#define CPU_RHOK_SCOPE_H
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* Explicit compatibility prototype. No method enum or dimensions change.
 * 1: native RMM-DIISK with scalar memory optimizations and native stopping. */
static int CpuNSocRhokMode(void)
{
 static int mode=-1;
 if(mode<0){
  const char *v=getenv("OPENMX_CPU_NSOC_RHOK_MODE");
  mode=!v||!strcmp(v,"0")?0:!strcmp(v,"1")?1:2;
  int lo,hi;MPI_Allreduce(&mode,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(&mode,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(lo!=hi||hi>1){fprintf(stderr,"CPU_NSOC_RHOK rejects inconsistent mode; require0,1\n");MPI_Abort(mpi_comm_level1,96);exit(96);}
 }
 if(mode){
  const char *p=getenv("OPENMX_CPU_NSOC_FULL_PROJECT"),*g=getenv("OPENMX_CPU_INPUT_SYMMETRY_GROUP"),*r=getenv("OPENMX_CPU_REALSPACE_D3");
  if(MD_switch||Solver!=2||SpinP_switch||SO_switch||Cnt_switch||LNO_flag||
     Mixing_switch!=4||SCF_Control_Temp||Scf_RestartFromFile||Num_Mixing_pDM!=30||
     List_YOUSO[16]!=3||List_YOUSO[38]<32||!p||strcmp(p,"1")||!g||strcmp(g,"D3")||
     !r||(strcmp(r,"scf")&&strcmp(r,"1"))){
   fprintf(stderr,"CPU_NSOC_RHOK requires fixed fresh scalar D3 Cluster RMM-DIISK History30 native3DM slots\n");MPI_Abort(mpi_comm_level1,96);exit(96);
  }
 }
 return mode;
}
static int CpuNSocMemoryMixingAccepted(void)
{return Mixing_switch==1||(Mixing_switch==4&&CpuNSocRhokMode()!=0);}
#endif
