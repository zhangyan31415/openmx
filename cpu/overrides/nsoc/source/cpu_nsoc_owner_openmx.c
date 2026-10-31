#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "openmx_common.h"
#include "flpq_dm.h"
#include "cpu_nsoc_owner.h"
#include "cpu_nsoc_owner_openmx.h"
#include "cpu_hassembly_workspace.h"
static CpuNOPlan*plan=NULL;static int cached_n=0,reported=0;
static void die(const char*why){fprintf(stderr,"CPU_NSOC_OWNER_DM_REJECT %s\n",why);MPI_Abort(mpi_comm_level1,98);abort();}
int CpuNSOCOwnerMode(void){const char*s=getenv("OPENMX_CPU_NSOC_OWNER_DM");if(!s||!strcmp(s,"0"))return 0;if(!strcmp(s,"1"))return 1;die("flag must be0 or1");return 0;}
void CpuNSOCOwnerValidate(int md_iter){
 int local[2]={CpuNSOCOwnerMode(),MD_switch==0&&md_iter==1&&Solver==2&&SpinP_switch==0&&SO_switch==0&&Cnt_switch==0&&Mixing_switch==1&&SCF_Control_Temp==0&&LNO_flag==0&&cal_partial_charge==0&&xanes_calc==0&&CWF_Calc==0&&flag_export_DM==0&&MO_fileout==0&&CDDF_on==0&&LNAO_calc_flag==0&&LNBO_calc_flag==0&&NBO_switch==0&&Dos_fileout==0&&DosGauss_fileout==0&&Hub_U_switch==0&&core_hole_state_flag==0},lo[2],hi[2];
 MPI_Allreduce(local,lo,2,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(local,hi,2,MPI_INT,MPI_MAX,mpi_comm_level1);
 if(lo[0]!=hi[0]||(hi[0]&&!lo[1]))die("requires agreed flag, fixed spin-off NSOC Cluster DM-mixing, no contraction/LNO/partial-charge/CWF/MO/DOS/DMmu/XANES/Hubbard/corehole");
}
static CpuNOPlan*representative_plan=NULL;
void CpuNSOCOwnerDestroy(void){CpuNODestroy(representative_plan);representative_plan=NULL;CpuNODestroy(plan);plan=NULL;cached_n=0;reported=0;}
static double CpuNSOCOwnerDensityLegacy(int n,int first,int last,int global_entries,const double*coeff,const double*eigenvalues,double*****dm,double*****edm,MPI_Comm world){
 double start=MPI_Wtime();int rank,np,cmp;MPI_Comm_rank(mpi_comm_level1,&rank);MPI_Comm_size(mpi_comm_level1,&np);MPI_Comm_compare(world,mpi_comm_level1,&cmp);if(cmp!=MPI_IDENT&&cmp!=MPI_CONGRUENT)die("owner rank ordering differs from eigenvector communicator");
 if(!plan){
  size_t count=0,entries=0;int ao=0;
  for(int g=1;g<=atomnum;g++){count+=(size_t)FNAN[g]+1;ao+=Spe_Total_CNO[WhatSpecies[g]];}
  if(count>INT_MAX||ao!=n||global_entries<0)die("invalid fixed basis/block count");
  CpuNOBlock*b=malloc((count?count:1)*sizeof(*b));int*offset=malloc(((size_t)atomnum+1)*sizeof(int));if(!b||!offset)die("block plan allocation");
  ao=0;for(int g=1;g<=atomnum;g++){offset[g]=ao;ao+=Spe_Total_CNO[WhatSpecies[g]];}
  size_t at=0;
  for(int g=1;g<=atomnum;g++){int nr=Spe_Total_CNO[WhatSpecies[g]];for(int h=0;h<=FNAN[g];h++){int q=natn[g][h],nc=Spe_Total_CNO[WhatSpecies[q]];b[at++]=(CpuNOBlock){offset[g],nr,offset[q],nc,G2ID[g]};entries+=(size_t)nr*nc;}}
  if(entries!=(size_t)global_entries)die("block payload differs from native global sparse size");
  plan=CpuNOCreate(world,n,(int)count,b,1024*1024);free(offset);free(b);if(!plan)die("invalid owner tile plan");cached_n=n;
 }
 if(n!=cached_n)die("basis changed without owner plan reset");
 int nk=last-first+1;if(nk<0)nk=0;double*f=malloc((size_t)(nk?nk:1)*sizeof(double));if(!f)die("occupation allocation");
 int po=0;for(int kk=0;kk<nk;kk++){int k=first+kk;double x=(eigenvalues[k]-ChemP)*Beta;if(x<=-60.)x=-60.;if(x>=60.)x=60.;f[kk]=FermiFunc(x,0,k,&po,&x);}
 if(CpuNOCompute(plan,nk,coeff,f,eigenvalues+first))die("owner density evaluation failed");free(f);
 const double*d=CpuNODensity(plan),*e=CpuNOEnergyDensity(plan);size_t at=0;
 for(int g=1;g<=atomnum;g++)if(G2ID[g]==rank){int m=F_G2M[g],nr=Spe_Total_CNO[WhatSpecies[g]];if(m<1||m>Matomnum||M2G[m]!=g)die("owned atom map mismatch");for(int h=0;h<=FNAN[g];h++){int nc=Spe_Total_CNO[WhatSpecies[natn[g][h]]];for(int i=0;i<nr;i++){memcpy(dm[0][m][h][i],d+at,(size_t)nc*sizeof(double));memcpy(edm[0][m][h][i],e+at,(size_t)nc*sizeof(double));at+=(size_t)nc;}}}
 if(at!=CpuNOEntries(plan))die("owned output count mismatch");
 if(!reported){unsigned long long bytes=CpuNOEntries(plan)*2*sizeof(double),lo,hi,work=CpuNOTemporaryBytes(plan),workmax;MPI_Reduce(&bytes,&lo,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,0,world);MPI_Reduce(&bytes,&hi,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,0,world);MPI_Reduce(&work,&workmax,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,0,world);if(!rank){printf("CPU_NSOC_OWNER_DM active=1 global_DM_EDM_replicas=0 retained_HS_Work1_elements=%d owned_DM_EDM_bytes_min=%llu max=%llu tile_work_bytes_max=%llu signed_EDM=1 occupations_unchanged=1\n",(CpuHAssembly_Compact()?1:global_entries),lo,hi,workmax);fflush(stdout);}reported=1;}
 return MPI_Wtime()-start;
}

#include <math.h>
#include "cpu_nsoc_representatives.h"
static int representative_calls=0,representative_cached_mode=-1;
static int representative_mode(void){
 if(representative_cached_mode>=0)return representative_cached_mode;
 const char*s=getenv("OPENMX_CPU_NSOC_REPRESENTATIVE_DM");int mode=0;
 if(s&&strcmp(s,"0")){if(!strcmp(s,"1"))mode=1;else if(!strcmp(s,"2"))mode=2;else die("REPRESENTATIVE_DM must be 0(off),1(on),2(full-owner audit)");}
 int lo,hi;MPI_Allreduce(&mode,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(&mode,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);if(lo!=hi)die("representative flag differs across ranks");representative_cached_mode=mode;return mode;
}
#include "cpu_nsoc_rep_tail.h"
static int representative_tail_cached=-1;
static int representative_tail_mode(void){
 if(representative_tail_cached>=0)return representative_tail_cached;
 const char*s=getenv("OPENMX_CPU_NSOC_REP_TAIL_BOUND");int v=0,bad=0,lo,hi;
 if(s&&strcmp(s,"0")){if(!strcmp(s,"1"))v=1;else bad=1;}
 int q[2]={v,bad},low[2],high[2];MPI_Allreduce(q,low,2,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(q,high,2,MPI_INT,MPI_MAX,mpi_comm_level1);
 if(high[1]||low[0]!=high[0])die("REP_TAIL_BOUND requires collectively agreed 0 or1");representative_tail_cached=v;return v;
}
double CpuNSOCOwnerDensity(int n,int first,int last,int global_entries,const double*coeff,const double*eigenvalues,double*****dm,double*****edm,MPI_Comm world){
 int mode=representative_mode(),tail=representative_tail_mode();
 if(!mode){if(tail)die("REP_TAIL_BOUND requires REPRESENTATIVE_DM");return CpuNSOCOwnerDensityLegacy(n,first,last,global_entries,coeff,eigenvalues,dm,edm,world);}
 if(empty_occupation_flag||empty_states_flag)die("representative density forbids state-index-dependent occupation overrides");
 representative_calls++;double start=MPI_Wtime();int rank,np,cmp,maxn;MPI_Comm_rank(world,&rank);MPI_Comm_size(world,&np);MPI_Comm_compare(world,mpi_comm_level1,&cmp);if(cmp!=MPI_IDENT&&cmp!=MPI_CONGRUENT)die("representative communicator mismatch");
 MPI_Allreduce(&last,&maxn,1,MPI_INT,MPI_MAX,world);
 int exact=CpuRSym_RepresentativeCertificate(n,maxn,eigenvalues),lo,hi;MPI_Allreduce(&exact,&lo,1,MPI_INT,MPI_MIN,world);MPI_Allreduce(&exact,&hi,1,MPI_INT,MPI_MAX,world);if(lo!=hi)die("representative solver certificate differs across ranks");
 int admitted=exact,candidate=0,cut_count=0;const int*cut_ids=NULL;double amplification=0.;CpuRepTailEvidence evidence={0};
 if(!exact&&tail){candidate=CpuRSym_RepresentativeTailCandidate(n,maxn,eigenvalues,&cut_count,&cut_ids,&amplification);MPI_Allreduce(&candidate,&lo,1,MPI_INT,MPI_MIN,world);MPI_Allreduce(&candidate,&hi,1,MPI_INT,MPI_MAX,world);if(lo!=hi)die("representative tail candidate differs across ranks");}
 if(!exact&&!candidate){if(!rank){printf("CPU_NSOC_REP_FALLBACK step=%d reason=incomplete_paired_weighted_subspace original_owner_kernel=1\n",representative_calls);fflush(stdout);}return CpuNSOCOwnerDensityLegacy(n,first,last,global_entries,coeff,eigenvalues,dm,edm,world);}
 int nk=last-first+1;if(nk<0)nk=0;double*f=malloc((size_t)(nk?nk:1)*sizeof(double));if(!f)die("representative occupations allocation");
 int po=0;for(int k=0;k<nk;k++){double x=(eigenvalues[first+k]-ChemP)*Beta;if(x < -60.)x=-60.;if(x > 60.)x=60.;f[k]=FermiFunc(x,0,first+k,&po,&x);}
 if(candidate){admitted=CpuRepTailAdmit(world,n,first,last,coeff,f,eigenvalues+first,cut_count,cut_ids,amplification,1e-13,&evidence);if(!rank){printf("CPU_NSOC_REP_TAIL_BOUND step=%d cut_states=%d covered_states=%d D_max_entry_bound=%.12e EDM_max_entry_bound=%.12e amplification=%.9g tolerance=1e-13 invalid=%d admitted=%d\n",representative_calls,evidence.cut_states,evidence.covered_states,evidence.density_bound,evidence.energy_density_bound,evidence.amplification,evidence.invalid,evidence.admitted);fflush(stdout);}}
 if(!admitted){free(f);if(!rank){printf("CPU_NSOC_REP_FALLBACK step=%d reason=weighted_tail_bound_not_satisfied original_owner_kernel=1\n",representative_calls);fflush(stdout);}return CpuNSOCOwnerDensityLegacy(n,first,last,global_entries,coeff,eigenvalues,dm,edm,world);}
 if(mode==2)CpuNSOCOwnerDensityLegacy(n,first,last,global_entries,coeff,eigenvalues,dm,edm,world);
 if(!representative_plan){int count;CpuNOBlock*b=CpuRSym_RepresentativeBlocks(&count);representative_plan=CpuNOCreate(world,n,count,b,1024*1024);free(b);if(!representative_plan)die("representative owner plan allocation");}
 if(CpuNOCompute(representative_plan,nk,coeff,f,eigenvalues+first))die("representative D/EDM kernel failed");free(f);
 CpuRSym_RepresentativeExpand(CpuNODensity(representative_plan),CpuNOEntries(representative_plan),dm[0]);
 CpuRSym_RepresentativeExpand(CpuNOEnergyDensity(representative_plan),CpuNOEntries(representative_plan),edm[0]);
 if(mode==2){
  double error=0,maximum;size_t at=0;const double*d=CpuNODensity(plan),*e=CpuNOEnergyDensity(plan);
  for(int g=1;g<=atomnum;g++)if(G2ID[g]==rank){int m=F_G2M[g],nr=Spe_Total_CNO[WhatSpecies[g]];for(int h=0;h<=FNAN[g];h++){int nc=Spe_Total_CNO[WhatSpecies[natn[g][h]]];for(int i=0;i<nr;i++)for(int j=0;j<nc;j++,at++){if(!isfinite(d[at])||!isfinite(e[at])||!isfinite(dm[0][m][h][i][j])||!isfinite(edm[0][m][h][i][j]))die("nonfinite reference or reconstructed density during audit");error=fmax(error,fabs(d[at]-dm[0][m][h][i][j]));error=fmax(error,fabs(e[at]-edm[0][m][h][i][j]));}}}
  MPI_Allreduce(&error,&maximum,1,MPI_DOUBLE,MPI_MAX,world);if(!rank){printf("CPU_NSOC_REP_AUDIT step=%d DM_EDM_max_abs=%.12e tolerance=5e-11\n",representative_calls,maximum);fflush(stdout);}if(!isfinite(maximum)||maximum>5e-11)die("representative density disagrees with independent full-owner kernel");
 }
 if(!rank){printf("CPU_NSOC_REP_ACTIVE step=%d certificate=%s signed_EDM=1 audit=%d elapsed=%.9f\n",representative_calls,exact?"complete_paired_spectral_subspace":"weighted_tail_bound",mode==2,MPI_Wtime()-start);fflush(stdout);}return MPI_Wtime()-start;
}
