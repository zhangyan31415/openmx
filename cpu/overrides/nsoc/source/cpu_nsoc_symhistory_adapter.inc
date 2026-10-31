#include "cpu_nsoc_symhistory_rs.h"
#include "cpu_symhistory_diis.h"
static CpuSymHistoryDIIS *cpu_history_plan=NULL;
static int cpu_history_last=0;
static void CpuNSocSymHistoryInverse(int n,double **a,double **ia,void *context)
{
 (void)context;Inverse(n,a,ia);
}
static void CpuNSocSymHistoryDie(const char *msg)
{
 fprintf(stderr,"CPU_NSOC_SYMHISTORY_REJECT %s\n",msg);MPI_Abort(mpi_comm_level1,96);exit(96);
}
void CpuNSocSymHistoryDestroy(void)
{
 CpuSymHistoryDIIS_Free(cpu_history_plan);cpu_history_plan=NULL;cpu_history_last=0;
}
static void CpuNSocSymHistoryAdmit(int step,double ******residual)
{
 if(cpu_history_plan)CpuNSocSymHistoryDie("unexpected DIIS reset with live History30");
 size_t n=CpuRSym_HistorySize();
 cpu_history_plan=CpuSymHistoryDIIS_Create(n,CpuRSym_HistoryWeights(),30,mpi_comm_level1);
 if(!cpu_history_plan)CpuNSocSymHistoryDie("representative history allocation failed");
 for(int slot=0;slot<=2;slot++){
  CpuRSym_HistoryPackChecked(DM[slot][0],CpuSymHistoryDIIS_Density(cpu_history_plan,slot),step,"admit_DM");
  CpuRSym_HistoryPackChecked(residual[slot][0],CpuSymHistoryDIIS_Residual(cpu_history_plan,slot),step,"admit_R");
 }
 cpu_history_last=step;
}
static void CpuNSocSymHistoryStep(int step,double ******residual)
{
 if(!cpu_history_plan||step!=cpu_history_last+1)CpuNSocSymHistoryDie("nonsequential DIIS step or missing native admission");
 double start=MPI_Wtime();
 /* Preserve native current residual and current two input vectors at full size.
  * Full native higher slots are NULL, owned by the compact engine instead. */
 for(int a=1;a<=Matomnum;a++){
  int ga=M2G[a],nr=Spe_Total_CNO[WhatSpecies[ga]];
  for(int b=0;b<=FNAN[ga];b++){
   int gb=natn[ga][b],nc=Spe_Total_CNO[WhatSpecies[gb]];
   for(int i=0;i<nr;i++)for(int j=0;j<nc;j++){
    double r=DM[0][0][a][b][i][j]-DM[1][0][a][b][i][j];
    residual[1][0][a][b][i][j]=r;
    residual[2][0][a][b][i][j]=r;
    DM[2][0][a][b][i][j]=DM[1][0][a][b][i][j];
   }
  }
 }
 /* Admission is checked at every new raw solver input. The tolerance admits
  * machine precision symmetry noise and rejects finite complement. */
 CpuRSym_HistoryPackChecked(DM[0][0],CpuSymHistoryDIIS_Density(cpu_history_plan,0),step,"solver_DM");
 CpuRSym_HistoryPack(DM[1][0],CpuSymHistoryDIIS_Density(cpu_history_plan,1));
 /* Check native residual admission independently. Kernel subtracts these same
  * native representative entries in the original per-element order. */
 CpuRSym_HistoryPackChecked(residual[1][0],CpuSymHistoryDIIS_Residual(cpu_history_plan,1),step,"native_R1");
 CpuSymHistoryDIISInfo info;
 int status=CpuSymHistoryDIIS_Step(cpu_history_plan,step,Mixing_weight,CpuNSocSymHistoryInverse,NULL,&info);
 if(status)CpuNSocSymHistoryDie("compressed native-equivalent DIIS kernel failed");
 NormRD[0]=info.norm_current;
 CpuRSym_HistoryExpand(CpuSymHistoryDIIS_Density(cpu_history_plan,0),DM[0][0]);
 for(int a=1;a<=Matomnum;a++){
  int ga=M2G[a],nr=Spe_Total_CNO[WhatSpecies[ga]];
  for(int b=0;b<=FNAN[ga];b++){
   int gb=natn[ga][b],nc=Spe_Total_CNO[WhatSpecies[gb]];
   for(int i=0;i<nr;i++)for(int j=0;j<nc;j++)DM[1][0][a][b][i][j]=DM[0][0][a][b][i][j];
  }
 }
 cpu_history_last=step;
 double elapsed=MPI_Wtime()-start,wallmax;MPI_Reduce(&elapsed,&wallmax,1,MPI_DOUBLE,MPI_MAX,Host_ID,mpi_comm_level1);
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 if(rank==Host_ID){printf("CPU_NSOC_SYMHISTORY_STEP effective=%d history=%d mix=%d slide=%d NormRD=%.12e optimized_norm=%.12e fallback=%d wall_max=%.9f\n",step,Num_Mixing_pDM,info.num_mix,info.num_slide,info.norm_current,info.norm_optimized,info.fallback,wallmax);fflush(stdout);}
}
