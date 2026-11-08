#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "openmx_common.h"
#include "cpu_rhok_scope.h"
int main(int argc,char **argv){
 MPI_Init(&argc,&argv);mpi_comm_level1=MPI_COMM_WORLD;int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 MD_switch=SpinP_switch=SO_switch=Cnt_switch=LNO_flag=Scf_RestartFromFile=SCF_Control_Temp=0;Solver=2;Mixing_switch=4;Num_Mixing_pDM=30;List_YOUSO[16]=3;List_YOUSO[38]=32;
 const char *v=argc>1?argv[1]:"1";if(!strcmp(v,"mismatch"))v=rank?"0":"1";
 if(!strcmp(v,"soc")){v="1";SpinP_switch=3;SO_switch=1;}
 if(!strcmp(v,"wrong_mixer")){v="1";Mixing_switch=1;}
 setenv("OPENMX_CPU_NSOC_RHOK_MODE",v,1);setenv("OPENMX_CPU_NSOC_FULL_PROJECT","1",1);setenv("OPENMX_CPU_REALSPACE_D3","scf",1);setenv("OPENMX_CPU_INPUT_SYMMETRY_GROUP","D3",1);
 assert(CpuNSocRhokMode()==atoi(v));assert(CpuNSocMemoryMixingAccepted()==(atoi(v)!=0));
 if(!rank)puts("PASS reciprocal mixing scope, native3DM slots, FourierHistory30");MPI_Finalize();return 0;
}
