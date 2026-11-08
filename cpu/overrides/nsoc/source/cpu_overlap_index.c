#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "openmx_common.h"
#include "cpu_rhok_scope.h"
#include "cpu_grid_shared.h"
#include "cpu_overlap_index.h"
extern int flag_export_DM;
CpuOverlapPair **CpuOverlapPairs=NULL;
uint64_t CpuOverlapGeneration=0;
static int overlap_atoms=0;
static int *overlap_neighbors=NULL;

/* Auxiliary paths not converted to index access retain their native arrays.
 * Fixed one-step SOC opts into the same exact geometry-only representation.
 * Moving geometry, stress, orbital exports and response retain native arrays. */
int CpuOverlapEligible(void)
{
  const char *setting=getenv("OPENMX_CPU_SOC_MEMORY");
  if(setting && strcmp(setting,"0") && strcmp(setting,"1")){
    fprintf(stderr,"OPENMX_CPU_SOC_MEMORY requires 0 or 1\n");
    MPI_Abort(mpi_comm_level1,95);abort();
  }
  int scalar=SpinP_switch==0 && SO_switch==0 && Scf_RestartFromFile==0;
  int soc=setting && !strcmp(setting,"1") && SpinP_switch==3 &&
    SO_switch==1 && DFTSCF_loop==1;
  return CpuGridShared_Requested() && MD_switch==0 && Solver==2 &&
    (scalar || soc) && Cnt_switch==0 && LNO_flag==0 &&
    (Mixing_switch==1||(scalar&&CpuNSocMemoryMixingAccepted())) &&
    !Band_disp_switch && !MO_fileout && !Dos_fileout && !DosGauss_fileout &&
    !Wannier_Func_Calc && !Wannier90_fileout && !CWF_fileout_flag && !CWF_Calc &&
    !pop_anal_aow_flag && !unfold_electronic_band && !Polarization_ABC_flag &&
    !xanes_gs_fileout && !xanes_calc && !scf_stress_flag && !stress_flag &&
    !StressConsistency_flag && !Calc_ABC_flag &&
    !LNAO_calc_flag && !LNBO_calc_flag && !NBO_switch && !NHO_fileout &&
    !NBO_fileout && !CDDF_on && !Hub_U_switch && !Voronoi_OrbM_flag &&
    !core_hole_state_flag && !cal_partial_charge && !flag_export_DM;
}
int CpuSOCMemoryEligible(void)
{
  return SpinP_switch==3 && SO_switch==1 && CpuOverlapEligible();
}
static void fail(const char *why)
{
  fprintf(stderr,"CPU_OVERLAP_INDEX %s\n",why);fflush(stderr);
  MPI_Abort(mpi_comm_level1,95);abort();
}
void CpuOverlapDestroy(void)
{
  if(CpuOverlapPairs) {
    for(int a=0;a<=overlap_atoms;a++) {
      for(int h=0;h<=overlap_neighbors[a];h++)free(CpuOverlapPairs[a][h].runs);
      free(CpuOverlapPairs[a]);
    }
    free(CpuOverlapPairs);free(overlap_neighbors);
  }
  CpuOverlapPairs=NULL;overlap_neighbors=NULL;overlap_atoms=0;
  ++CpuOverlapGeneration;
}
void CpuOverlapInit(void)
{
  CpuOverlapDestroy();
  int eligible=CpuOverlapEligible(),lo,hi;
  MPI_Allreduce(&eligible,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(&eligible,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(lo!=hi)fail("inconsistent eligibility");
  if(!lo)return;
  overlap_atoms=Matomnum;
  CpuOverlapPairs=calloc((size_t)Matomnum+1,sizeof(*CpuOverlapPairs));
  overlap_neighbors=calloc((size_t)Matomnum+1,sizeof(*overlap_neighbors));
  if(!CpuOverlapPairs || !overlap_neighbors)fail("metadata allocation failed");
  for(int a=0;a<=Matomnum;a++) {
    int n=a?FNAN[M2G[a]]:0;
    overlap_neighbors[a]=n;
    CpuOverlapPairs[a]=calloc((size_t)n+1,sizeof(**CpuOverlapPairs));
    if(!CpuOverlapPairs[a])fail("pair allocation failed");
  }
}
void CpuOverlapStore(int mc,int h,int n,const int *nc,const int *nh)
{
  if(n<0)fail("negative overlap length");
  CpuOverlapPair *p=CpuOverlapPairs?&CpuOverlapPairs[mc][h]:NULL;
  int runs=p?CpuOverlapRunCount(n,nc,nh):0;
  if(p)p->count=n;
  if(p && n>0 && (size_t)runs*sizeof(CpuOverlapRun)<2*(size_t)n*sizeof(int)) {
    p->runs=malloc((size_t)runs*sizeof(*p->runs));
    if(!p->runs)fail("run allocation failed");
    p->nruns=runs;CpuOverlapRunEncode(n,nc,nh,p->runs);
    GListTAtoms1[mc][h]=NULL;GListTAtoms2[mc][h]=NULL;
  } else {
    GListTAtoms1[mc][h]=malloc((n?(size_t)n:1)*sizeof(int));
    GListTAtoms2[mc][h]=malloc((n?(size_t)n:1)*sizeof(int));
    if(!GListTAtoms1[mc][h] || !GListTAtoms2[mc][h])fail("raw pair allocation failed");
    if(n) {memcpy(GListTAtoms1[mc][h],nc,(size_t)n*sizeof(int));memcpy(GListTAtoms2[mc][h],nh,(size_t)n*sizeof(int));}
  }
}
void CpuOverlapReport(void)
{
  unsigned long long local[4]={0,0,0,0},sum[4];int rank;
  if(!CpuOverlapPairs)return;
  for(int a=1;a<=overlap_atoms;a++)for(int h=0;h<=overlap_neighbors[a];h++) {
    CpuOverlapPair *p=&CpuOverlapPairs[a][h];
    local[0]+=2*(unsigned long long)p->count*sizeof(int);
    local[1]+=p->runs?(unsigned long long)p->nruns*sizeof(CpuOverlapRun):2*(unsigned long long)p->count*sizeof(int);
    local[2]+=p->runs?1:0;local[3]+=p->runs?p->nruns:0;
  }
  MPI_Comm_rank(mpi_comm_level1,&rank);
  MPI_Reduce(local,sum,4,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,mpi_comm_level1);
  if(rank==0) {printf("CPU_OVERLAP_INDEX active=1 raw_pair_bytes=%llu stored_pair_bytes=%llu compressed_pairs=%llu runs=%llu precision=exact order=native\n",sum[0],sum[1],sum[2],sum[3]);fflush(stdout);}
}
