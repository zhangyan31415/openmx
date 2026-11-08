#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "openmx_common.h"
#include "cpu_rhok_scope.h"
#include "flpq_dm.h"
#include "cpu_nsoc_owner_openmx.h"
#include "cpu_nsoc_imag.h"
static int stored=-1;
static void reject(const char *why)
{
  fprintf(stderr,"CPU_NSOC_IMAG_HISTORY_REJECT %s\n",why);
  MPI_Abort(mpi_comm_level1,99); abort();
}
int CpuNSOCImagSlots(void)
{
  if (stored>=0) return stored==0 ? List_YOUSO[16] : stored;
  const char *s=getenv("OPENMX_CPU_NSOC_IMAG_COMPACT");
  int requested=0;
  if(s && strcmp(s,"0")) {if(strcmp(s,"1")) reject("flag must be 0 or 1"); requested=1;}
  int owner=CpuNSOCOwnerMode();
  int eligible=MD_switch==0 && Solver==2 && SpinP_switch==0 && SO_switch==0 &&
    Cnt_switch==0 && CpuNSocMemoryMixingAccepted() && SCF_Control_Temp==0 && LNO_flag==0 &&
    cal_partial_charge==0 && xanes_calc==0 && CWF_Calc==0 && flag_export_DM==0 &&
    MO_fileout==0 && CDDF_on==0 && LNAO_calc_flag==0 && LNBO_calc_flag==0 &&
    NBO_switch==0 && Dos_fileout==0 && DosGauss_fileout==0 && Hub_U_switch==0 &&
    core_hole_state_flag==0 && Constraint_NCS_switch==0 && Zeeman_NCS_switch==0 &&
    Zeeman_NCO_switch==0;
  int local[3]={requested,CpuNSOCImagSelect(requested,owner,eligible,List_YOUSO[16]),List_YOUSO[16]},lo[3],hi[3];
  MPI_Allreduce(local,lo,3,MPI_INT,MPI_MIN,mpi_comm_level1);
  MPI_Allreduce(local,hi,3,MPI_INT,MPI_MAX,mpi_comm_level1);
  if(lo[0]!=hi[0] || lo[2]!=hi[2] || lo[1]<0) reject("requires agreed flag/capacity and fixed spin-off NSOC Cluster owner DM-mixing without auxiliary spin/DM features");
  stored=requested?1:0;
  return stored==0 ? List_YOUSO[16] : stored;
}
void CpuNSOCImagReport(unsigned long long local_payload_bytes)
{
  int slots=CpuNSOCImagSlots();
  if(slots==List_YOUSO[16]) return;
  unsigned long long actual=0;
  MPI_Reduce(&local_payload_bytes,&actual,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,Host_ID,mpi_comm_level1);
  int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
  if(rank==Host_ID){
    unsigned long long legacy=actual*(unsigned long long)List_YOUSO[16];
    printf("CPU_NSOC_IMAG_HISTORY active=1 requested_slots=%d stored_slots=%d real_history_unchanged=1 full_zero_current_retained=1 payload_bytes=%llu eliminated_payload_bytes=%llu\n",List_YOUSO[16],slots,actual,legacy-actual);fflush(stdout);
  }
}
