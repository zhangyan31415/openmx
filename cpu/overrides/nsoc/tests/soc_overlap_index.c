/* Exercise the production overlap store/decoder with fixed one-step SOC.
 * The fixture substitutes only the shared-window row allocator. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "openmx_common.h"
#include "cpu_overlap_index.h"
static float rows[32][3];
int flag_export_DM;
int CpuGridShared_Requested(void){return 1;}
int CpuGridShared_IsAtom(int atom){return atom==2;}
float *CpuGridShared_ReadRow(int atom,int row){assert(atom==2&&row>=0&&row<32);return rows[row];}
int main(int argc,char **argv){
 MPI_Init(&argc,&argv);mpi_comm_level1=MPI_COMM_WORLD;
 MD_switch=0;Solver=2;SpinP_switch=3;SO_switch=1;DFTSCF_loop=1;
 Mixing_switch=1;Scf_RestartFromFile=1;
 setenv("OPENMX_CPU_SOC_MEMORY","1",1);
 assert(CpuOverlapEligible());
 DFTSCF_loop=2;assert(!CpuOverlapEligible());DFTSCF_loop=1;
 scf_stress_flag=1;assert(!CpuOverlapEligible());scf_stress_flag=0;
 Band_disp_switch=1;assert(!CpuOverlapEligible());Band_disp_switch=0;
 setenv("OPENMX_CPU_SOC_MEMORY","0",1);assert(!CpuOverlapEligible());
 setenv("OPENMX_CPU_SOC_MEMORY","1",1);
 atomnum=2;Matomnum=1;M2G=calloc(2,sizeof(int));M2G[1]=1;
 FNAN=calloc(3,sizeof(int));FNAN[1]=2;natn=calloc(3,sizeof(int*));
 natn[1]=calloc(3,sizeof(int));natn[1][0]=1;natn[1][1]=natn[1][2]=2;
 GListTAtoms1=calloc(2,sizeof(*GListTAtoms1));GListTAtoms2=calloc(2,sizeof(*GListTAtoms2));
 for(int a=0;a<=1;a++){int n=a?3:1;GListTAtoms1[a]=calloc(n,sizeof(**GListTAtoms1));GListTAtoms2[a]=calloc(n,sizeof(**GListTAtoms2));}
 CpuOverlapInit();assert(CpuOverlapPairs);
 int nc[16],nh[16];for(int i=0;i<16;i++){nc[i]=i;nh[i]=8+i;}
 CpuOverlapStore(1,1,16,nc,nh);assert(CpuOverlapPairs[1][1].runs);
 /* A second periodic image of the same atom has its own exact lookup. */
 int nc2[16],nh2[16];for(int i=0;i<16;i++){nc2[i]=2*i;nh2[i]=(i*7)%31;}
 CpuOverlapStore(1,2,16,nc2,nh2);assert(!CpuOverlapPairs[1][2].runs);
 for(int i=15;i>=0;i--){
  assert(CpuOverlapIndex1(1,1,i)==nc[i]&&CpuOverlapIndex2(1,1,i)==nh[i]);
  assert(CpuOverlapIndex1(1,2,i)==nc2[i]&&CpuOverlapIndex2(1,2,i)==nh2[i]);
  assert(CpuOverlapRow(1,1,i)==rows[nh[i]]);
  assert(CpuOverlapRow(1,2,i)==rows[nh2[i]]);
 }
 CpuOverlapReport();CpuOverlapDestroy();assert(!CpuOverlapPairs);
 SpinP_switch=0;SO_switch=0;Scf_RestartFromFile=0;
 setenv("OPENMX_CPU_SOC_MEMORY","0",1);assert(CpuOverlapEligible());
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);if(!rank)puts("PASS SOC exact overlap indices, periodic images, shared rows and safe scope");
 MPI_Finalize();return 0;
}
