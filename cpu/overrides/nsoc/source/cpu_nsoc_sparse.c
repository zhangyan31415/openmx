/* Distributed sparse orthogonal-basis congruence. No replicated dense N*N
 * workspace: one <=64-row/column panel plus one row scratch per rank. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "cpu_nsoc_sparse.h"
#define NSOC_SPARSE_MAX_BLOCK 64
struct CpuNSparsePlan {
 int n,nb,nr,nc,pr,pc,row,col,components,maxblock;
 int *ptr,*rows,*cols,*row_slot;
 const int64_t *cp;const int*ri;const double*cv;
 MPI_Comm rowcomm,colcomm;
 double *panel,*scratch;size_t panel_bytes,metadata_bytes;
 double *back_panel;size_t back_panel_bytes;int back_panel_rows,back_batches;
};
static int fail(char*e,size_t z,const char*s){if(e&&z)snprintf(e,z,"%s",s);return -1;}
static int root(int*p,int a){while(p[a]!=a){p[a]=p[p[a]];a=p[a];}return a;}
static void join(int*p,int a,int b){a=root(p,a);b=root(p,b);if(a!=b){if(a<b)p[b]=a;else p[a]=b;}}
static int owner(int g,int b,int np){return(g/b)%np;}
static int local_index(int g,int b,int np){return(g/(b*np))*b+g%b;}
static int local_count(int n,int b,int c,int np){int k=0;for(int g=0;g<n;g++)if(owner(g,b,np)==c)k++;return k;}
CpuNSparsePlan *CpuNSparse_Create(int n,int nb,int nr,int nc,int pr,int pc,int row,int col,const int64_t*cp,const int*ri,const double*cv,MPI_Comm world,char*error,size_t ez){
 int size;MPI_Comm_size(world,&size);
 if(n<1||n>INT_MAX/2||nb<1||pr<1||pc<1||(int64_t)pr*pc!=size||row<0||row>=pr||col<0||col>=pc||!cp||!ri||!cv){fail(error,ez,"invalid sparse congruence grid/basis");return NULL;}
 if(nr!=local_count(n,nb,row,pr)||nc!=local_count(n,nb,col,pc)){fail(error,ez,"actual local dimensions do not match block-cyclic map");return NULL;}
 if(cp[0]!=0||cp[n]<1||cp[n]>INT_MAX){fail(error,ez,"invalid sparse CSC extent");return NULL;}
 CpuNSparsePlan*p=calloc(1,sizeof(*p));int*uf=malloc((size_t)2*n*4),*map=malloc((size_t)2*n*4),*counts=calloc((size_t)2*n,4),*pos=calloc((size_t)2*n,4),*coords=malloc((size_t)2*size*4);
 if(!p||!uf||!map||!counts||!pos||!coords){fail(error,ez,"sparse plan allocation failed");free(p);free(uf);free(map);free(counts);free(pos);free(coords);return NULL;}
 p->rowcomm=p->colcomm=MPI_COMM_NULL;p->n=n;p->nb=nb;p->nr=nr;p->nc=nc;p->pr=pr;p->pc=pc;p->row=row;p->col=col;p->cp=cp;p->ri=ri;p->cv=cv;
 for(int k=0;k<2*n;k++){uf[k]=k;map[k]=-1;}
 int me[2]={row,col};MPI_Allgather(me,2,MPI_INT,coords,2,MPI_INT,world);
 for(int k=0;k<size;k++){int r=coords[2*k],c=coords[2*k+1];if(r<0||r>=pr||c<0||c>=pc){fail(error,ez,"invalid actual BLACS coordinate");goto bad;}for(int q=0;q<k;q++)if(coords[2*q]==r&&coords[2*q+1]==c){fail(error,ez,"duplicate actual BLACS coordinate");goto bad;}}
 memset(pos,0,(size_t)2*n*4);
 for(int j=0;j<n;j++){
  if(cp[j]<0||cp[j]>cp[j+1]||cp[j+1]>cp[n]){fail(error,ez,"invalid sparse CSC pointers");goto bad;}
  for(int64_t z=cp[j];z<cp[j+1];z++){if(ri[z]<0||ri[z]>=n||!isfinite(cv[z])){fail(error,ez,"sparse CSC row/value invalid");goto bad;}if(cv[z]!=0)join(uf,ri[z],n+j);}
 }
 for(int i=0;i<n;i++){int q=root(uf,i);if(map[q]<0)map[q]=p->components++;counts[map[q]]++;}
 for(int j=0;j<n;j++){int q=map[root(uf,n+j)];if(q<0){fail(error,ez,"empty sparse output column");goto bad;}counts[n+q]++;}
 p->ptr=calloc((size_t)p->components+1,4);p->rows=malloc((size_t)n*4);p->cols=malloc((size_t)n*4);p->row_slot=malloc((size_t)n*4);
 if(!p->ptr||!p->rows||!p->cols||!p->row_slot){fail(error,ez,"sparse component allocation failed");goto bad;}
 for(int k=0;k<p->components;k++){
  if(counts[k]!=counts[n+k]||counts[k]>NSOC_SPARSE_MAX_BLOCK){fail(error,ez,"sparse-U component must be square and at most64 (use dense backend for wider components)");goto bad;}
  if(counts[k]>p->maxblock)p->maxblock=counts[k];p->ptr[k+1]=p->ptr[k]+counts[k];pos[k]=pos[n+k]=p->ptr[k];
 }
 for(int i=0;i<n;i++){int q=map[root(uf,i)],at=pos[q]++;p->rows[at]=i;p->row_slot[i]=at-p->ptr[q];}
 for(int j=0;j<n;j++){int q=map[root(uf,n+j)];p->cols[pos[n+q]++]=j;}
 int dim=nr>nc?nr:nc;if(dim<1)dim=1;
 if((size_t)(p->maxblock+1)>SIZE_MAX/8/(size_t)dim||(int64_t)p->maxblock*dim>INT_MAX){fail(error,ez,"sparse panel/MPI count overflow");goto bad;}
 p->panel_bytes=(size_t)(p->maxblock+1)*dim*8;p->panel=malloc(p->panel_bytes);if(!p->panel){fail(error,ez,"bounded sparse panel allocation failed");goto bad;}p->scratch=p->panel+(size_t)p->maxblock*dim;
 p->metadata_bytes=sizeof(*p)+(size_t)(p->components+1)*4+(size_t)3*n*4;
 MPI_Comm_split(world,row,col,&p->rowcomm);MPI_Comm_split(world,col,row,&p->colcomm);
 free(uf);free(map);free(counts);free(pos);free(coords);return p;
 bad:free(uf);free(map);free(counts);free(pos);free(coords);CpuNSparse_Destroy(p);return NULL;
}
int CpuNSparse_Congruence(CpuNSparsePlan*p,double*A,int lda,double*W,int ldw,int profile,CpuNSparseTimings*t,char*error,size_t ez){
 if(!p||!A||!W||A==W||lda<(p->nr?p->nr:1)||ldw<(p->nr?p->nr:1))return fail(error,ez,"invalid sparse congruence buffers/leading dimensions");
 CpuNSparseTimings tt={0};double total=profile?MPI_Wtime():0,start=0;
 for(int component=0;component<p->components;component++){
  int begin=p->ptr[component],end=p->ptr[component+1],m=end-begin,nr=p->nr;
  if(profile)start=MPI_Wtime();memset(p->panel,0,(size_t)m*nr*8);
  for(int a=0;a<m;a++){int g=p->rows[begin+a];if(owner(g,p->nb,p->pc)==p->col){int j=local_index(g,p->nb,p->pc);for(int i=0;i<nr;i++)p->panel[(size_t)a*nr+i]=A[(size_t)j*lda+i];}}
  if(profile)tt.right_pack+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  if(nr&&p->pc>1)MPI_Allreduce(MPI_IN_PLACE,p->panel,m*nr,MPI_DOUBLE,MPI_SUM,p->rowcomm);
  if(profile)tt.right_collective+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  for(int a=begin;a<end;a++){int g=p->cols[a];if(owner(g,p->nb,p->pc)==p->col){double*out=W+(size_t)local_index(g,p->nb,p->pc)*ldw;for(int i=0;i<nr;i++)out[i]=0;for(int64_t z=p->cp[g];z<p->cp[g+1];z++){double v=p->cv[z];if(v==0)continue;const double*in=p->panel+(size_t)p->row_slot[p->ri[z]]*nr;for(int i=0;i<nr;i++)out[i]+=v*in[i];}}}
  if(profile)tt.right_multiply+=MPI_Wtime()-start;
 }
 for(int component=0;component<p->components;component++){
  int begin=p->ptr[component],end=p->ptr[component+1],m=end-begin,nc=p->nc;
  if(profile)start=MPI_Wtime();memset(p->panel,0,(size_t)m*nc*8);
  for(int a=0;a<m;a++){int g=p->rows[begin+a];if(owner(g,p->nb,p->pr)==p->row){int i=local_index(g,p->nb,p->pr);for(int j=0;j<nc;j++)p->panel[(size_t)a*nc+j]=W[(size_t)j*ldw+i];}}
  if(profile)tt.left_pack+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  if(nc&&p->pr>1)MPI_Allreduce(MPI_IN_PLACE,p->panel,m*nc,MPI_DOUBLE,MPI_SUM,p->colcomm);
  if(profile)tt.left_collective+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  for(int a=begin;a<end;a++){int g=p->cols[a];if(owner(g,p->nb,p->pr)==p->row){int i=local_index(g,p->nb,p->pr);for(int j=0;j<nc;j++)p->scratch[j]=0;for(int64_t z=p->cp[g];z<p->cp[g+1];z++){double v=p->cv[z];if(v==0)continue;const double*in=p->panel+(size_t)p->row_slot[p->ri[z]]*nc;for(int j=0;j<nc;j++)p->scratch[j]+=v*in[j];}for(int j=0;j<nc;j++)A[(size_t)j*lda+i]=p->scratch[j];}}
  if(profile)tt.left_multiply+=MPI_Wtime()-start;
 }
 if(profile)tt.total=MPI_Wtime()-total;if(t)*t=tt;return 0;
}
int CpuNSparse_Left(CpuNSparsePlan*p,const double*A,int lda,double*W,int ldw,int batch_rows,int profile,CpuNSparseTimings*t,char*error,size_t ez){
 if(!p||!A||!W||A==W||lda<(p->nr?p->nr:1)||ldw<(p->nr?p->nr:1))return fail(error,ez,"invalid sparse left-product buffers/leading dimensions");
 if(batch_rows<0||batch_rows>1024)return fail(error,ez,"backtransform batch rows must be0..1024");
 CpuNSparseTimings tt={0};double total=profile?MPI_Wtime():0,start=0;
 double coeff[NSOC_SPARSE_MAX_BLOCK*NSOC_SPARSE_MAX_BLOCK];
 int limit=batch_rows>p->maxblock?batch_rows:p->maxblock,nc=p->nc;
 double*panel=p->panel;
 if(batch_rows){
  int dim=nc?nc:1;
  if((size_t)limit>SIZE_MAX/sizeof(double)/(size_t)dim||(int64_t)limit*nc>INT_MAX)return fail(error,ez,"backtransform batch panel/MPI count overflow");
  if(!p->back_panel){
   size_t bytes=(size_t)limit*dim*sizeof(double);
   p->back_panel=malloc(bytes);if(!p->back_panel)return fail(error,ez,"backtransform batch panel allocation failed");
   p->back_panel_bytes=bytes;p->back_panel_rows=limit;
  }
  if(limit>p->back_panel_rows)return fail(error,ez,"backtransform batch capacity changed; use a fresh process");
  panel=p->back_panel;
 }
 p->back_batches=0;
 for(int first=0,last;first<p->components;first=last){
  int batch_begin=p->ptr[first];last=first+1;
  if(batch_rows)while(last<p->components&&p->ptr[last+1]-batch_begin<=limit)last++;
  int batch_end=p->ptr[last],rows=batch_end-batch_begin;p->back_batches++;
  if(profile)start=MPI_Wtime();memset(panel,0,(size_t)rows*nc*sizeof(double));
  /* U's columns are symmetry-basis indices, hence source rows of E.
   * Group adjacent whole components; only communication is batched. */
  for(int a=batch_begin;a<batch_end;a++){int g=p->cols[a];if(owner(g,p->nb,p->pr)==p->row){int i=local_index(g,p->nb,p->pr);for(int j=0;j<nc;j++)panel[(size_t)(a-batch_begin)*nc+j]=A[(size_t)j*lda+i];}}
  if(profile)tt.left_pack+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  if(nc&&p->pr>1)MPI_Allreduce(MPI_IN_PLACE,panel,rows*nc,MPI_DOUBLE,MPI_SUM,p->colcomm);
  if(profile)tt.left_collective+=MPI_Wtime()-start;if(profile)start=MPI_Wtime();
  for(int component=first;component<last;component++){
   int begin=p->ptr[component],end=p->ptr[component+1],m=end-begin;
  /* A fixed <=64-by-64 stack block avoids repeated sparse searches. The
   * existing contiguous scratch row keeps arithmetic off strided W storage. */
  memset(coeff,0,(size_t)m*m*sizeof(double));
  for(int a=begin;a<end;a++){
   int g=p->cols[a];
   for(int64_t z=p->cp[g];z<p->cp[g+1];z++){
    double v=p->cv[z];if(v!=0)coeff[(size_t)p->row_slot[p->ri[z]]*m+a-begin]+=v;
   }
  }
  for(int a=begin;a<end;a++){int g=p->rows[a];if(owner(g,p->nb,p->pr)==p->row){
   int i=local_index(g,p->nb,p->pr);memset(p->scratch,0,(size_t)nc*sizeof(double));
   for(int b=0;b<m;b++){double v=coeff[(size_t)(a-begin)*m+b];if(v==0)continue;const double*in=panel+(size_t)(begin-batch_begin+b)*nc;for(int j=0;j<nc;j++)p->scratch[j]+=v*in[j];}
   for(int j=0;j<nc;j++)W[(size_t)j*ldw+i]=p->scratch[j];
  }}
  }
  if(profile)tt.left_multiply+=MPI_Wtime()-start;
 }
 if(profile)tt.total=MPI_Wtime()-total;if(t)*t=tt;return 0;
}
size_t CpuNSparse_BackPanelBytes(const CpuNSparsePlan*p){return p?p->back_panel_bytes:0;}
int CpuNSparse_BackBatches(const CpuNSparsePlan*p){return p?p->back_batches:0;}
size_t CpuNSparse_PanelBytes(const CpuNSparsePlan*p){return p?p->panel_bytes:0;}
size_t CpuNSparse_MetadataBytes(const CpuNSparsePlan*p){return p?p->metadata_bytes:0;}
int CpuNSparse_Components(const CpuNSparsePlan*p){return p?p->components:0;}
int CpuNSparse_MaxBlock(const CpuNSparsePlan*p){return p?p->maxblock:0;}
void CpuNSparse_Destroy(CpuNSparsePlan*p){if(!p)return;if(p->rowcomm!=MPI_COMM_NULL)MPI_Comm_free(&p->rowcomm);if(p->colcomm!=MPI_COMM_NULL)MPI_Comm_free(&p->colcomm);free(p->ptr);free(p->rows);free(p->cols);free(p->row_slot);free(p->panel);free(p->back_panel);free(p);}
