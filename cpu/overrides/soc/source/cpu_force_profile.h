#ifndef CPU_FORCE_PROFILE_H
#define CPU_FORCE_PROFILE_H
/* Profiling-only bounded counters: no heap/scientific arrays or math changes.
 * Coarse and nested labels are explicitly distinguished, never summed across
 * rank maxima. End-of-Force is the only profiling reduction point. */
enum {
 FP_TOTAL,FP_EARLY,FP_CACHE,FP_F1,FP_PUBLISH,FP_PRE2,FP_F2,FP_F3,FP_F4,FP_F5,
 FP_HNL,FP_HCH,FP_CLEANUP,FP_ENDWAIT,FP_F3_MASTER,
 FP_DER_THREAD_MIN,FP_DER_THREAD_MAX,FP_DER_WAIT_MIN,FP_DER_WAIT_MAX,
 FP_MASTER_WAIT_MIN,FP_MASTER_WAIT_MAX,
 FP_VNA1,FP_VNA2,FP_VNA3,FP_VNA4,FP_VNA5,FP_VNA5A,FP_VNA5B,
 FP_NL1,FP_NL2,FP_NL3,FP_NL_REMOTE_COMPUTE,FP_NL_SIZE_WAIT,FP_NL_PAYLOAD_WAIT,FP_NL_SEND_WAIT,FP_NL_PACK,FP_NL_UNPACK,FP_COUNT
};
static struct {int active,nt;double start;double value[FP_COUNT],thread[4][3];} CpuFP;
static const char *CpuFPName[FP_COUNT]={
 "force_total","early_release","cache_setup","force1","force1_publication_nested",
 "pre_force2_sync","force2","force3","force4","force5","nonlocal_force",
 "corehole_force","cleanup","end_barrier","force3_master_contract_nested",
 "force3_derivative_thread_min","force3_derivative_thread_max",
 "force3_derivative_wait_thread_min","force3_derivative_wait_thread_max",
 "force3_master_wait_thread_min","force3_master_wait_thread_max",
 "vna_part1_setup_nested","vna_part2_weight_nested","vna_part3_remote_nested",
 "vna_part4_local_nested","vna_part5_center_nested","vna_part5A_route_nested",
 "vna_part5B_compute_nested","nonlocal_part1_setup_nested",
 "nonlocal_part2_basis_nested","nonlocal_part3_center_nested",
 "nonlocal_remote_compute_nested","nonlocal_size_MPI_wait_nested",
 "nonlocal_payload_MPI_recv_wait_nested","nonlocal_payload_MPI_send_wait_nested",
 "nonlocal_payload_pack_nested","nonlocal_payload_unpack_nested"
};
static void CpuFPBegin(void){
 const char*s=getenv("OPENMX_CPU_FORCE_PROFILE");int mode=!s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:-1;
 int nt=omp_get_max_threads(),code=mode<0?-1:mode?(!omp_in_parallel()&&(nt==1||nt==4)?1+2*nt:-1):0;
 int in[2]={code,-code},out[2];MPI_Allreduce(in,out,2,MPI_INT,MPI_MIN,mpi_comm_level1);
 if(out[0]<0||out[0]!=-out[1])CpuF3TileFail("FORCE_PROFILE must agree and be exactly0/1; ON requires serial entry withOMP1or4");
 CpuFP.active=mode;if(!mode)return;
 memset(CpuFP.value,0,sizeof(CpuFP.value));memset(CpuFP.thread,0,sizeof(CpuFP.thread));CpuFP.nt=0;CpuFP.start=omp_get_wtime();
}
static double CpuFPTick(void){return CpuFP.active?omp_get_wtime():0.;}
static void CpuFPAdd(int field,double start){if(CpuFP.active)CpuFP.value[field]+=omp_get_wtime()-start;}
static void CpuFPValue(int field,double seconds){if(CpuFP.active)CpuFP.value[field]+=seconds;}
static void CpuFPThread(int thread,int field,double seconds){if(CpuFP.active)CpuFP.thread[thread][field]+=seconds;}
static void CpuFPFinish(void){
 if(!CpuFP.active)return;
 CpuFP.value[FP_TOTAL]=omp_get_wtime()-CpuFP.start;
 for(int k=0;k<3;k++){
  double lo=0.,hi=0.;if(CpuFP.nt){lo=hi=CpuFP.thread[0][k];for(int t=1;t<CpuFP.nt;t++){if(CpuFP.thread[t][k]<lo)lo=CpuFP.thread[t][k];if(CpuFP.thread[t][k]>hi)hi=CpuFP.thread[t][k];}}
  int at=k==0?FP_DER_THREAD_MIN:k==1?FP_DER_WAIT_MIN:FP_MASTER_WAIT_MIN;CpuFP.value[at]=lo;CpuFP.value[at+1]=hi;
 }
 struct {double value;int rank;} local[2*FP_COUNT],global[2*FP_COUNT];int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 for(int k=0;k<FP_COUNT;k++){local[2*k].value=CpuFP.value[k];local[2*k].rank=rank;local[2*k+1].value=-CpuFP.value[k];local[2*k+1].rank=rank;}
 /* Standard C double/int pair matches MPI_DOUBLE_INT. No collectives inside
  * worker/atom loops; positive and negative MAXLOC give max/min in one call. */
 MPI_Reduce(local,global,2*FP_COUNT,MPI_DOUBLE_INT,MPI_MAXLOC,Host_ID,mpi_comm_level1);
 if(rank==Host_ID){for(int k=0;k<FP_COUNT;k++)printf("CPU_FORCE_PROFILE stage=%s min=%.9f rank_min=%d max=%.9f rank_max=%d spread=%.9f nested=%d\n",CpuFPName[k],-global[2*k+1].value,global[2*k+1].rank,global[2*k].value,global[2*k].rank,global[2*k].value+global[2*k+1].value,k==FP_PUBLISH||k>=FP_F3_MASTER);fflush(stdout);}
 printf("CPU_FORCE_PROFILE_THREADS rank=%d threads=%d",rank,CpuFP.nt);
 for(int t=0;t<CpuFP.nt;t++)printf(" t%d_derivative=%.9f t%d_derivative_wait=%.9f t%d_master_wait=%.9f",t,CpuFP.thread[t][0],t,CpuFP.thread[t][1],t,CpuFP.thread[t][2]);
 printf("\n");fflush(stdout);CpuFP.active=0;
}
#endif
