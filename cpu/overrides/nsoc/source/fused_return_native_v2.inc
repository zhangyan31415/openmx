#include "fused_return.h"
#include "selected_left.h"
static int CfsNativeRequest(void){
 static int cached=-1;if(cached>=0)return cached;const char*s=getenv("OPENMX_CPU_NSOC_FUSED_RETURN");int mode=!s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:!strcmp(s,"2")?2:-1,in[2]={mode,-mode},out[2];MPI_Allreduce(in,out,2,MPI_INT,MPI_MIN,mpi_comm_level1);
 if(out[0]<0||out[0]!=-out[1])die("FUSED_RETURN must agree and be0,1,2");
 if(mode){int ok=d3_basis()&&SpinP_switch==0&&Cnt_switch==0&&CpuNSym_FullProjected()&&option("OPENMX_CPU_NSOC_PARTIAL_EVP",0)==0&&option("OPENMX_CPU_NSOC_SPARSE_BACK",0)==1&&option("OPENMX_CPU_NSOC_CHECK",0)==0,all;MPI_Allreduce(&ok,&all,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!all)die("FUSED_RETURN requires full-projected real nonmagnetic uncontracted D3/fullEVP/SPARSE_BACK1/CHECK0");}
 cached=mode;if(rank0()){printf("CPU_NSOC_FUSED_RETURN_CONFIG mode=%d value_only=1 returned_state_rule_unchanged=1\n",mode);fflush(stdout);}return mode;
}
static int CfsNativeSelectedRequest(int fused){
 static int cached=-1;if(cached>=0)return cached;const char*s=getenv("OPENMX_CPU_NSOC_SELECTED_BACK");int mode=!s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:!strcmp(s,"2")?2:-1,in[2]={mode,-mode},out[2];MPI_Allreduce(in,out,2,MPI_INT,MPI_MIN,mpi_comm_level1);
 if(out[0]<0||out[0]!=-out[1])die("SELECTED_BACK must agree and be0,1,2");if(mode&&(fused!=1||option("OPENMX_CPU_NSOC_BACK_BATCH_ROWS",0)!=0))die("SELECTED_BACK requires FUSED_RETURN1 and BACK_BATCH_ROWS0");cached=mode;return mode;
}
static int*CfsNativePermutation(const double*ev,double*ko){Eigen*s=mem(B.N,sizeof(Eigen));int*inv=mem(B.N,sizeof(int));for(int i=0;i<B.N;i++){if(!isfinite(ev[i]))die("nonfinite eigenvalue");s[i].value=ev[i];s[i].source=i;}qsort(s,B.N,sizeof(Eigen),eigcmp);for(int i=0;i<B.N;i++){ko[i+1]=s[i].value;inv[s[i].source]=i;}free(s);return inv;}
static CfsSelectedPlan*CfsNativeSelectedPlan(void){
 static CfsSelectedPlan*p=NULL;static int size,block,rows,cols;static const int64_t*columns;
 if(p){if(size!=B.N||block!=nblk||rows!=na_rows||cols!=na_cols||columns!=B.cp)die("selected sparse plan invalidation");return p;}
 char error[256];p=CfsSelectedCreate(B.N,nblk,na_rows,na_cols,np_rows,np_cols,my_prow,my_pcol,B.cp,B.ri,B.cv,mpi_comm_level1,error,sizeof(error));if(!p)die(error);size=B.N;block=nblk;rows=na_rows;cols=na_cols;columns=B.cp;return p;
}
static void CfsNativeTail(int step,int n,int maxn,double*eval,double*ko,double*H,double*W,int phase,int mode){
 int one=1,selected=CfsNativeSelectedRequest(mode);double a=1.,b=0.,start=MPI_Wtime(),left=0.,reference=0.;CfsStats stats={0};int*inv=CfsNativePermutation(eval,ko);
 if(mode==2){backtransform(step,n,H,W,phase);reorder(eval,ko,maxn,W,H);pdtran_(&n,&n,&a,W,&one,&one,descC,&b,B.U,&one,&one,descC);reference=MPI_Wtime()-start;}
 else{
  for(int j=0;j<na_cols;j++){int gj=global(j,nblk,my_pcol,np_cols);if(gj>=n)continue;for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi<n&&sector(gi)!=sector(gj))H[(size_t)j*na_rows+i]=0.;}}
  int batch=(int)option("OPENMX_CPU_NSOC_BACK_BATCH_ROWS",0);if(batch<0||batch>1024)die("invalid BACK_BATCH_ROWS");char error[256];CpuNSparseTimings oldtimes={0};double plan_time=0.;
  if(!selected){CpuNSparsePlan*p=sparse_plan_get(&plan_time);if(CpuNSparse_Left(p,H,na_rows,B.U,na_rows,batch,phase,&oldtimes,error,sizeof(error)))die(error);}
  else{
   int*columns=mem(na_cols?na_cols:1,sizeof(int)),count=0;for(int j=0;j<na_cols;j++){int q=global(j,nblk,my_pcol,np_cols);if(q<n&&inv[q]<maxn)columns[count++]=j;}
   if(selected==2){CpuNSparsePlan*old=sparse_plan_get(&plan_time);if(CpuNSparse_Left(old,H,na_rows,W,na_rows,0,phase,&oldtimes,error,sizeof(error)))die(error);}
   CfsSelectedTimes timing={0};if(CfsSelectedLeft(CfsNativeSelectedPlan(),H,na_rows,B.U,na_rows,columns,count,phase,&timing,error,sizeof(error)))die(error);
   double audit_error=0.;
   if(selected==2){double v[2]={0.,0.},g[2];int bad=0,gbad;for(int c=0;c<count;c++)for(int i=0;i<na_rows;i++){size_t at=(size_t)columns[c]*na_rows+i;double x=W[at],y=B.U[at];if(!isfinite(x)||!isfinite(y)){bad=1;continue;}if(fabs(x-y)>v[0])v[0]=fabs(x-y);if(fabs(x)>v[1])v[1]=fabs(x);}MPI_Allreduce(v,g,2,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);MPI_Allreduce(&bad,&gbad,1,MPI_INT,MPI_MAX,mpi_comm_level1);audit_error=g[0];if(gbad||g[0]>3e-13*fmax(1.,g[1]))die("selected sparse product all-selected-entry audit failed");}
   if(rank0()){printf("CPU_NSOC_SELECTED_BACK step=%d mode=%d n=%d maxn=%d original_local_cols=%d selected_local_cols=%d rank0_plan_bytes=%zu audit_max_abs=%.17g\n",step,selected,n,maxn,na_cols,count,timing.plan_bytes,audit_error);fflush(stdout);}
   if(phase){int rank;MPI_Comm_rank(mpi_comm_level1,&rank);printf("CPU_NSOC_SELECTED_BACK_PHASE rank=%d step=%d columns=%d pack=%.9f collective=%.9f multiply=%.9f logical_values=%llu\n",rank,step,count,timing.pack,timing.collective,timing.multiply,timing.logical_values);fflush(stdout);}free(columns);
  }
  left=MPI_Wtime()-start;
 }
 int tile=(int)option("OPENMX_CPU_NSOC_REORDER_TILE",1024);if(tile<1||tile>n)tile=n;double*output=mode==2?W:H;
 if(CfsTransposeSelected(mpi_comm_level1,n,maxn,nblk,np_rows,np_cols,my_prow,my_pcol,B.rankmap,inv,B.U,na_rows,output,na_rows,tile,&stats))die("fused return validation failed");free(inv);double maxerror=0.;
 if(mode==2){double v[2]={0.,0.},g[2];int bad=0,gbad;for(size_t i=0;i<B.active;i++){double x=H[i],y=W[i];if(!isfinite(x)||!isfinite(y)){bad=1;continue;}if(fabs(x-y)>v[0])v[0]=fabs(x-y);if(fabs(x)>v[1])v[1]=fabs(x);}MPI_Allreduce(v,g,2,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);MPI_Allreduce(&bad,&gbad,1,MPI_INT,MPI_MAX,mpi_comm_level1);maxerror=g[0];if(gbad||maxerror!=0.)die("fused return all-entry audit failed");}
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);if(!rank){printf("CPU_NSOC_FUSED_RETURN step=%d mode=%d n=%d maxn=%d logical_global_values=%llu batches=%d rank0_panel_bytes=%zu audit_max_abs=%.17g\n",step,mode,n,maxn,(unsigned long long)n*(unsigned)maxn,stats.batches,stats.peak_bytes,maxerror);fflush(stdout);}
 if(phase){printf("CPU_NSOC_FUSED_RETURN_PHASE rank=%d step=%d mode=%d left=%.9f reference=%.9f pack=%.9f collective=%.9f unpack=%.9f total=%.9f sent_values=%llu received_values=%llu panel_bytes=%zu\n",rank,step,mode,left,reference,stats.pack,stats.collective,stats.unpack,MPI_Wtime()-start,stats.sent_values,stats.received_values,stats.peak_bytes);fflush(stdout);}
}
