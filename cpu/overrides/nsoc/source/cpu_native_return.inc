#include "sector_native_pack.h"
#include "state_coefficients.h"
static CpuStateCoefficients native_states;
static int native_states_retired_step=0;
static void native_state_release(void){CpuStateDestroy(&native_states);native_states_retired_step=0;}
int CpuNSym_StateProvider(int n,int first,int last,CpuNORowProvider*provider){
 if(native_return_mode()<3)return 0;
 if(!provider||native_return.returned!=native_return.step||native_states.n!=n||native_states.first!=first||native_states.last!=last)die("stale native state coefficient provider");
 *provider=(CpuNORowProvider){&native_states,CpuStateRows};return 1;
}
int CpuNSym_StateCutNorms(int first,int last,int count,const int*ids,long double*norms){
 if(native_return_mode()<3||native_return.returned!=native_return.step||!native_states.data||
    first!=native_states.first||last!=native_states.last)return -1;
 return CpuStateCutNorms(&native_states,count,ids,norms);
}
int CpuNSym_StateProviderAudit(double*error){if(native_return_mode()<3)return 0;if(error)*error=native_states.audit_error;return native_return_mode()==4;}
void CpuNSym_RetireStateCoefficients(int step){
 if(native_return_mode()<3)return;
 if(step==native_return.step&&step==native_states_retired_step)return;
 if(step<1||step!=native_return.step||native_return.returned!=step||
    native_states.data!=B.U||!borrowed_eigen_scratch||B.U!=eigen_scratch.ptr)
  die("state coefficient retirement requires the current completed borrowed-Q solve");
 size_t live=native_states.used*sizeof(double);
 size_t reference=native_states.reference?(size_t)native_states.n*native_states.nk*sizeof(double):0;
 /* Q is dead, including its boundary page inside the larger reserved EVec1
  * payload. The shared reclaim helper protects allocation boundary pages. */
 size_t discarded=CpuReclaim_Pages(B.U,0,eigen_scratch.capacity*sizeof(double));
 free(native_states.reference);native_states.reference=NULL;
 native_states.data=NULL;native_return.returned=0;native_states_retired_step=step;
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 printf("CPU_NSOC_STATE_COEFFICIENTS_RETIRED rank=%d step=%d live_Q_bytes=%zu discarded_allocation_page_bytes=%zu audit_reference_bytes_freed=%zu allocation_address_retained=1 sparse_U_row_metadata_retained=1 late_provider_reads_rejected=1\n",
  rank,step,live,discarded,reference);fflush(stdout);
}
static double*CpuNativeStateReference(int n,int maxn,double*H,double*W,const int*inv,const CpuNativeSector*sectors,const int*ends,int phase){
 size_t offsets[48],used=0,temporary=0;
 if(CpuNativePack(mpi_comm_level1,n,nblk,na_rows,na_cols,np_rows,np_cols,my_prow,my_pcol,B.rankmap,sectors,B.nsec,H,W,B.cap,offsets,&used,1024,&temporary))die("state audit pack failed");
 double*full=mem(B.active,sizeof(double)),*ao=mem(B.cap,sizeof(double));
 int rs[48],cs[48],rn[48],r=0,c=0,previous=0;
 for(int s=0;s<B.nsec;s++){
  rs[s]=r;cs[s]=c;for(int q=previous;q<ends[s];q++){r+=(q/nblk)%np_rows==my_prow;c+=(q/nblk)%np_cols==my_pcol;}rn[s]=r-rs[s];previous=ends[s];
 }
 for(int j=0;j<na_cols;j++){
  int gj=global(j,nblk,my_pcol,np_cols);if(gj>=n)continue;int s=sector(gj);
  for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi<n&&sector(gi)==s)full[(size_t)j*na_rows+i]=W[offsets[s]+(size_t)(j-cs[s])*rn[s]+i-rs[s]];}
 }
 int*columns=mem(na_cols?na_cols:1,sizeof(int)),count=0;for(int j=0;j<na_cols;j++){int q=global(j,nblk,my_pcol,np_cols);if(q<n&&inv[q]<maxn)columns[count++]=j;}
 CfsSelectedTimes timing={0};char error[256];double audit=0.;
 if(CfsSelectedLeft(CfsNativeSelectedPlan(),full,na_rows,ao,na_rows,columns,count,phase,&timing,error,sizeof(error)))die(error);
 int nk=native_return.last>=native_return.first?native_return.last-native_return.first+1:0;
 double*reference=mem((size_t)n*(nk?nk:1),sizeof(double));
 if(CfsSelectedNativeLeft(CfsNativeSelectedPlan(),mpi_comm_level1,W,na_rows,inv,maxn,native_return.first,native_return.last,reference,(size_t)n*nk,ends,B.nsec,offsets,1024,ao,phase,&timing,&audit,error,sizeof(error)))die(error);
 if(audit!=0.)die("state provider reference differs from original sparse product");
 free(columns);free(full);free(ao);return reference;
}
static void CpuNativeStateReturnTail(int step,int n,int maxn,double*eval,double*ko,double*H,double*W,int phase){
 double start=MPI_Wtime();int*inv=CfsNativePermutation(eval,ko),ends[48];CpuNativeSector sectors[48];
 for(int k=0;k<B.nsec;k++){Sector*t=B.sec+k;ends[k]=t->ofs+t->n;sectors[k]=(CpuNativeSector){t->n,t->ofs,t->nr,t->nc,t->nb,B.hoff[k]};}
 if(CpuStateInit(&native_states,n,B.cp,B.ri,B.cv)||
    CpuStateBind(&native_states,B.nsec,ends,inv,maxn,native_return.first,native_return.last,B.U,eigen_scratch.capacity))die("state coefficient descriptor initialization failed");
 if(native_return_mode()==4)native_states.reference=CpuNativeStateReference(n,maxn,H,W,inv,sectors,ends,phase);
 size_t old_w=CpuReclaim_Pages(W,0,B.cap*sizeof(double));
 size_t old_u=CpuReclaim_Pages(B.U,0,eigen_scratch.capacity*sizeof(double));size_t temporary=0;
 if(CpuStatePack(mpi_comm_level1,np_rows,np_cols,my_prow,my_pcol,sectors,H,inv,maxn,&native_states,1024,&temporary))die("state-owned compact coefficient redistribution failed");
 size_t returned_h=CpuReclaim_Pages(H,0,B.active*sizeof(double));
 CpuReclaim_Pages(B.U,native_states.used*sizeof(double),eigen_scratch.capacity*sizeof(double));
 native_return.returned=step;free(inv);int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 printf("CPU_NSOC_STATE_COEFFICIENTS rank=%d step=%d mode=%d n=%d maxn=%d local_states=%d Q_bytes=%zu U_row_metadata_bytes=%zu native_full_AO_coefficients=0 H_pages_returned=%zu old_W_pages_returned=%zu old_U_pages_returned=%zu conversion_temporary_bytes=%zu reference_bytes=%zu states_occupations_unchanged=1 elapsed=%.9f\n",
  rank,step,native_return_mode(),n,maxn,native_states.nk,native_states.used*sizeof(double),
  (size_t)(n+1)*sizeof(int)+(size_t)native_states.rowptr[n]*(sizeof(int)+sizeof(double)),returned_h,old_w,old_u,temporary,
  native_states.reference?(size_t)n*native_states.nk*sizeof(double):0,MPI_Wtime()-start);fflush(stdout);
}
/* H contains only reset-origin packed sectors; W's sector solve scratch is
 * dead. B.U still aliases dead EVec1 scratch. Pack H into W, release H and
 * B.U pages, then write native coefficients directly. No input aliases output. */
static void CpuNativeReturnTail(int step,int n,int maxn,double*eval,double*ko,double*H,double*W,int phase)
{
 if(native_return_mode()>=3){CpuNativeStateReturnTail(step,n,maxn,eval,ko,H,W,phase);return;}
 int*inv=CfsNativePermutation(eval,ko),ends[48];size_t offsets[48],used=0,temporary=0;
 CpuNativeSector sectors[48];double start=MPI_Wtime();
 for(int k=0;k<B.nsec;k++){
  Sector*t=B.sec+k;sectors[k]=(CpuNativeSector){t->n,t->ofs,t->nr,t->nc,t->nb,B.hoff[k]};ends[k]=t->ofs+t->n;
 }
 if(CpuNativePack(mpi_comm_level1,n,nblk,na_rows,na_cols,np_rows,np_cols,my_prow,my_pcol,
    B.rankmap,sectors,B.nsec,H,W,B.cap,offsets,&used,1024,&temporary))die("packed sector/native-coordinate redistribution failed");
 double pack_end=MPI_Wtime();size_t h_returned=CpuReclaim_Pages(H,0,B.active*sizeof(double));
 size_t w_returned=CpuReclaim_Pages(W,used*sizeof(double),B.cap*sizeof(double));
 size_t u_returned=CpuReclaim_Pages(B.U,0,eigen_scratch.capacity*sizeof(double));
 double*reference=NULL;int mode=native_return_mode();char error[256];
 if(mode==2){
  /* Independent previous sparse multiplication oracle. Audit-only full H/Q
   * arrays intentionally exclude this run from memory/timing measurements. */
  memset(H,0,B.active*sizeof(double));
  int row_start[48],col_start[48],row_count[48],r=0,c=0,previous=0;
  for(int s=0;s<B.nsec;s++){
   row_start[s]=r;col_start[s]=c;for(int q=previous;q<ends[s];q++){r+=(q/nblk)%np_rows==my_prow;c+=(q/nblk)%np_cols==my_pcol;}
   row_count[s]=r-row_start[s];previous=ends[s];
  }
  for(int j=0;j<na_cols;j++){
   int gj=global(j,nblk,my_pcol,np_cols);if(gj>=n)continue;int s=sector(gj);
   for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi<n&&sector(gi)==s)H[(size_t)j*na_rows+i]=W[offsets[s]+(size_t)(j-col_start[s])*row_count[s]+i-row_start[s]];}
  }
  reference=mem(B.cap,sizeof(double));int*columns=mem(na_cols?na_cols:1,sizeof(int)),count=0;
  for(int j=0;j<na_cols;j++){int q=global(j,nblk,my_pcol,np_cols);if(q<n&&inv[q]<maxn)columns[count++]=j;}
  CfsSelectedTimes old={0};if(CfsSelectedLeft(CfsNativeSelectedPlan(),H,na_rows,reference,na_rows,columns,count,phase,&old,error,sizeof(error)))die(error);
  free(columns);
 }
 int tile=(int)option("OPENMX_CPU_NSOC_REORDER_TILE",1024);if(tile<64)tile=64;if(tile>1024)tile=1024;
 CfsSelectedTimes timing={0};double audit_error=0.;
 if(CfsSelectedNativeLeft(CfsNativeSelectedPlan(),mpi_comm_level1,W,na_rows,inv,maxn,
     native_return.first,native_return.last,native_return.coeff,native_return.capacity,
     ends,B.nsec,offsets,tile,reference,phase,&timing,&audit_error,error,sizeof(error)))die(error);
 if(mode==2&&audit_error!=0.)die("native return sparse-product audit differs from previous exact path");
 free(reference);free(inv);
 if(mode==2)h_returned+=CpuReclaim_Pages(H,0,B.active*sizeof(double));
 native_return.returned=step;
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 printf("CPU_NSOC_NATIVE_RETURN rank=%d step=%d mode=%d n=%d maxn=%d packed_H_bytes=%zu packed_back_input_bytes=%zu H_pages_returned=%zu eigen_scratch_pages_returned=%zu W_tail_pages_returned=%zu native_coefficient_bytes=%zu conversion_temporary_bytes=%zu audit_max_abs=%.17g first_full_states_retained=1\n",
  rank,step,mode,n,maxn,B.hcap*sizeof(double),used*sizeof(double),h_returned,u_returned,w_returned,
  (size_t)n*(native_return.last>=native_return.first?native_return.last-native_return.first+1:0)*sizeof(double),temporary,audit_error);fflush(stdout);
 if(phase){printf("CPU_NSOC_NATIVE_RETURN_PHASE rank=%d step=%d pack_sector=%.9f sparse_pack=%.9f collective=%.9f multiply=%.9f total=%.9f temporary_bytes=%zu\n",
  rank,step,pack_end-start,timing.pack,timing.collective,timing.multiply,MPI_Wtime()-start,timing.plan_bytes);fflush(stdout);}
}
