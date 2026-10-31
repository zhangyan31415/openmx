/* Independent namespaced sparse plan: reuses audited component construction,
 * does NOT assume the ABI/layout of the separately linked legacy plan object. */
#define CpuNSparsePlan CfsSelectedPlan
#define CpuNSparse_Create CfsRefCreate
#define CpuNSparse_Congruence CfsRefCongruence
#define CpuNSparse_Left CfsRefLeft
#define CpuNSparse_BackPanelBytes CfsRefBackPanelBytes
#define CpuNSparse_BackBatches CfsRefBackBatches
#define CpuNSparse_PanelBytes CfsRefPanelBytes
#define CpuNSparse_MetadataBytes CfsRefMetadataBytes
#define CpuNSparse_Components CfsRefComponents
#define CpuNSparse_MaxBlock CfsRefMaxBlock
#define CpuNSparse_Destroy CfsRefDestroy
#include "cpu_nsoc_sparse.frozen.c"
#undef CpuNSparsePlan
#include "selected_left.h"
CfsSelectedPlan*CfsSelectedCreate(int n,int nb,int nr,int nc,int pr,int pc,int row,int col,const int64_t*cp,const int*ri,const double*cv,MPI_Comm comm,char*e,size_t z){return CfsRefCreate(n,nb,nr,nc,pr,pc,row,col,cp,ri,cv,comm,e,z);}
void CfsSelectedDestroy(CfsSelectedPlan*p){CfsRefDestroy(p);}
int CfsSelectedLeft(CfsSelectedPlan*p,const double*A,int lda,double*W,int ldw,const int*columns,int count,int profile,CfsSelectedTimes*t,char*error,size_t ez){
 if(!p)return fail(error,ez,"selected sparse plan missing");
 int okay=A&&W&&A!=W&&lda>=(p->nr?p->nr:1)&&ldw>=(p->nr?p->nr:1)&&count>=0&&count<=p->nc&&(!count||columns),all;
 int*seen=calloc(p->nc?p->nc:1,sizeof(int));if(!seen){MPI_Abort(p->colcomm,96);abort();}
 unsigned long long hash=1469598103934665603ULL;hash^=(unsigned)count;hash*=1099511628211ULL;
 if(okay)for(int j=0;j<count;j++){int q=columns[j];if(q<0||q>=p->nc)okay=0;else if(seen[q]++)okay=0;hash^=(unsigned)q;hash*=1099511628211ULL;}
 free(seen);unsigned long long in[2]={hash,~hash},out[2];MPI_Allreduce(in,out,2,MPI_UNSIGNED_LONG_LONG,MPI_MIN,p->colcomm);if(out[0]!=hash||out[1]!=~hash)okay=0;
 MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,p->colcomm);if(!all)return fail(error,ez,"selected column list/buffers differ within column communicator");
 CfsSelectedTimes tt={0};tt.columns=count;tt.plan_bytes=p->panel_bytes+p->metadata_bytes;double start=0.;
 double coeff[NSOC_SPARSE_MAX_BLOCK*NSOC_SPARSE_MAX_BLOCK];
 for(int component=0;component<p->components;component++){
  int begin=p->ptr[component],end=p->ptr[component+1],m=end-begin;tt.batches++;tt.logical_values+=(unsigned long long)m*(unsigned)count;
  if(!count)continue;
  if(profile)start=MPI_Wtime();memset(p->panel,0,(size_t)m*count*sizeof(double));
  for(int a=0;a<m;a++){int g=p->cols[begin+a];if(owner(g,p->nb,p->pr)==p->row){int i=local_index(g,p->nb,p->pr);for(int j=0;j<count;j++)p->panel[(size_t)a*count+j]=A[(size_t)columns[j]*lda+i];}}
  if(profile)tt.pack+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  if(p->pr>1)MPI_Allreduce(MPI_IN_PLACE,p->panel,m*count,MPI_DOUBLE,MPI_SUM,p->colcomm);
  if(profile)tt.collective+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  memset(coeff,0,(size_t)m*m*sizeof(double));
  for(int a=begin;a<end;a++){int g=p->cols[a];for(int64_t z=p->cp[g];z<p->cp[g+1];z++){double v=p->cv[z];if(v!=0.)coeff[(size_t)p->row_slot[p->ri[z]]*m+a-begin]+=v;}}
  for(int a=begin;a<end;a++){int g=p->rows[a];if(owner(g,p->nb,p->pr)==p->row){
   int i=local_index(g,p->nb,p->pr);memset(p->scratch,0,(size_t)count*sizeof(double));
   for(int b=0;b<m;b++){double v=coeff[(size_t)(a-begin)*m+b];if(v==0.)continue;const double*input=p->panel+(size_t)b*count;for(int j=0;j<count;j++)p->scratch[j]+=v*input[j];}
   for(int j=0;j<count;j++)W[(size_t)columns[j]*ldw+i]=p->scratch[j];
  }}
  if(profile)tt.multiply+=MPI_Wtime()-start;
 }
 if(t)*t=tt;return 0;
}
