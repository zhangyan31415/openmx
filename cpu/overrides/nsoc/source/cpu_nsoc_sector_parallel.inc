/* ELPA_SCOPE_PATCH_BEGIN */
#include "cpu_nsoc_elpa_omp_scope.h"
/* ELPA_SCOPE_PATCH_END */
/* Experimental real-D3 independent ELPA groups. No new dense allocations.
 * Source matrices stay in the caller world descriptors. Only W, borrowed B.U
 * and cached Ss are reinterpreted inside this explicit solver lifetime. */
#include "cpu_nsoc_sector_layout.h"
extern int Csys2blacs_handle(MPI_Comm);
extern void Cblacs_gridmap(int*,int*,int,int,int);
extern void Cblacs_gridinfo(int,int*,int*,int*,int*);
extern void Cblacs_gridexit(int);
typedef struct {
 int sector,p,pr,pc,row,col,member,root,ctxt,nb,nr,nc,ld,desc[9];
 size_t count;MPI_Comm comm,rows,cols;
} CpuNSectorGroup;
static struct {
 int loaded,migrated,n,assigned,rank,size,desc[9];
 int decision_ready,enabled,decision_n,decision_version,decision_natoms,counts[3],threads,mkl_threads;MPI_Comm decision_world;
 /* ELPA_SCOPE_PATCH_BEGIN */
 int elpa_omp_threads;unsigned long long elpa_omp_calls;
/* ELPA_SCOPE_PATCH_END */
 double *S,*H,*W,*U;CpuNSectorGroup group[3];
} sector_parallel;
static int sector_parallel_local_configuration(int*enabled,int*p){
 const char*flag=getenv("OPENMX_CPU_NSOC_SECTOR_PARALLEL");char*end;double q=0;
 if(flag&&*flag){q=strtod(flag,&end);if(end==flag||*end||!isfinite(q)||(q!=0&&q!=1))return 0;}
 *enabled=(int)q;const char*text=getenv("OPENMX_CPU_NSOC_SECTOR_RANKS");if(!text||!*text)text=sector_parallel.size==64?"32,16,16":"224,16,16";
 for(int k=0;k<3;k++){
  long n=strtol(text,&end,10);if(end==text||n<1||n>256)return 0;p[k]=(int)n;
  if(k<2){if(*end!=',')return 0;text=end+1;}else if(*end)return 0;
 }
 return 1;
}
static int sector_parallel_request(int step){
 if(!sector_parallel.decision_ready){
  int size;MPI_Comm_size(mpi_comm_level1,&size);sector_parallel.size=size;
  int enabled=0,counts[3]={0,0,0},threads=omp_get_max_threads(),mkl_threads=mkl_get_max_threads();int valid=sector_parallel_local_configuration(&enabled,counts);
  if(threads<1||threads>256||mkl_threads<1||mkl_threads>256)valid=0;
/* ELPA_SCOPE_PATCH_BEGIN */
  int elpa_omp_threads=0;if(!CpuNSocElpaOmpParse(getenv("OPENMX_CPU_NSOC_ELPA_OMP_THREADS"),&elpa_omp_threads))valid=0;
/* ELPA_SCOPE_PATCH_END */

  /* Uniquely encode the bounded flag/ranks/thread tuple. MIN(x,-x) simultaneously
   * gives the global min/max in ONE tiny collective, including the OFF path.
   * Parsing never aborts a subset before this common agreement point. */
  long long code=valid?((long long)enabled<<27)|((long long)counts[0]<<18)|((long long)counts[1]<<9)|counts[2]|((long long)threads<<28)|((long long)mkl_threads<<37)/* ELPA_SCOPE_PATCH_BEGIN */
|((long long)elpa_omp_threads<<46)
/* ELPA_SCOPE_PATCH_END */
:(1LL<<60);
  long long pair[2]={code,-code},extrema[2];MPI_Allreduce(pair,extrema,2,MPI_LONG_LONG_INT,MPI_MIN,mpi_comm_level1);
  if(extrema[0]!=-extrema[1])die("SECTOR_PARALLEL configuration differs across MPI ranks");
  if(!valid)die("invalid collective SECTOR_PARALLEL flag/ranks configuration");
  /* ELPA_SCOPE_PATCH_BEGIN */
  sector_parallel.elpa_omp_threads=elpa_omp_threads;
  if(elpa_omp_threads&&(!enabled||size!=64||threads!=4||mkl_threads!=4))die("ELPA_OMP_THREADS=1 requires enabled64MPIx4 private real-D3 path with MKL4");
  if(rank0()){printf("CPU_NSOC_ELPA_OMP_CONFIG threads=%d outer_omp=%d mkl_threads=%d uniform=1 cached=1\n",elpa_omp_threads,threads,mkl_threads);fflush(stdout);}
/* ELPA_SCOPE_PATCH_END */
  sector_parallel.enabled=enabled;sector_parallel.threads=threads;sector_parallel.mkl_threads=mkl_threads;memcpy(sector_parallel.counts,counts,sizeof(counts));sector_parallel.decision_n=B.N;sector_parallel.decision_version=B.version;sector_parallel.decision_natoms=B.natoms;sector_parallel.decision_world=mpi_comm_level1;sector_parallel.decision_ready=1;
  if(rank0()){printf("CPU_NSOC_SECTOR_PARALLEL_CONFIG enabled=%d ranks=%d,%d,%d omp_threads=%d mkl_threads=%d agreement=1 cached=1\n",enabled,counts[0],counts[1],counts[2],threads,mkl_threads);fflush(stdout);}
 }
 if(sector_parallel.decision_n!=B.N||sector_parallel.decision_version!=B.version||sector_parallel.decision_natoms!=B.natoms||sector_parallel.decision_world!=mpi_comm_level1)die("SECTOR_PARALLEL cached configuration geometry/basis/world changed");
 if(!sector_parallel.enabled)return 0;
 if(!d3_basis()||!CpuNSym_FullProjected()||option("OPENMX_CPU_NSOC_PARTIAL_EVP",0)!=0||raw_reset_at()!=0)
  die("SECTOR_PARALLEL requires real four-sector full-projected D3, PARTIAL_EVP0 and no raw reset");
 if(step!=1&&!sector_parallel.loaded)die("SECTOR_PARALLEL must start on the first SCF iteration");
 return 1;
}
static int sector_parallel_rows(int p){int row=(int)sqrt((double)p);while(row>1&&p%row)row--;return row;}
static void sector_parallel_ranks(int*p){
 if(!sector_parallel.decision_ready)die("SECTOR_PARALLEL ranks requested before collective agreement");
 memcpy(p,sector_parallel.counts,sizeof(sector_parallel.counts));if(p[0]+p[1]+p[2]!=sector_parallel.size)die("SECTOR_RANKS must sum to the actual MPI world size");
}
static void sector_parallel_prepare(double*S,double*H,double*W){
 if(sector_parallel.loaded){
  if(sector_parallel.n!=B.N||sector_parallel.S!=S||sector_parallel.H!=H||sector_parallel.W!=W||sector_parallel.U!=B.U||memcmp(sector_parallel.desc,descH,sizeof(sector_parallel.desc))||omp_get_max_threads()!=sector_parallel.threads||mkl_get_max_threads()!=sector_parallel.mkl_threads)
   die("SECTOR_PARALLEL buffer/distribution changed after initialization");
  return;
 }
 int rank,size;MPI_Comm_rank(mpi_comm_level1,&rank);MPI_Comm_size(mpi_comm_level1,&size);
 if(!((size==256&&sector_parallel.threads==1)||(size==64&&sector_parallel.threads==4))||sector_parallel.threads!=sector_parallel.mkl_threads)die("SECTOR_PARALLEL scope requires4nodes with256MPIx1 or64MPIx4");
 if(omp_get_max_threads()!=sector_parallel.threads||mkl_get_max_threads()!=sector_parallel.mkl_threads)die("SECTOR_PARALLEL actual thread counts changed after collective agreement");
 MPI_Comm localcomm;MPI_Comm_split_type(mpi_comm_level1,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&localcomm);
 int local_rank,local_size;MPI_Comm_rank(localcomm,&local_rank);MPI_Comm_size(localcomm,&local_size);
 int expected_local=64/sector_parallel.threads,leader=local_rank==0,node_count=0,minimum_rank,maximum_rank,all_ok=0;
 MPI_Allreduce(&rank,&minimum_rank,1,MPI_INT,MPI_MIN,localcomm);MPI_Allreduce(&rank,&maximum_rank,1,MPI_INT,MPI_MAX,localcomm);
 int okay=local_size==expected_local&&minimum_rank%expected_local==0&&maximum_rank-minimum_rank+1==expected_local;
 MPI_Allreduce(&leader,&node_count,1,MPI_INT,MPI_SUM,mpi_comm_level1);MPI_Allreduce(&okay,&all_ok,1,MPI_INT,MPI_MIN,mpi_comm_level1);
 if(!all_ok||!sector_scope_valid(size,node_count,local_size,sector_parallel.threads,sector_parallel.mkl_threads))die("SECTOR_PARALLEL requires4 contiguous shared-node rank domains and64cores per node");
 CpuNSectorThreadProof thread_proof;sector_thread_probe(sector_parallel.threads,&thread_proof);
 okay=sector_node_layout(localcomm,expected_local,sector_parallel.threads,64,&thread_proof);MPI_Allreduce(&okay,&all_ok,1,MPI_INT,MPI_MIN,mpi_comm_level1);
 if(!all_ok)die("SECTOR_PARALLEL actual OMP/MKL thread affinity does not cover64 distinct CPU cores per node");
 char hostname[MPI_MAX_PROCESSOR_NAME];int hostname_length;MPI_Get_processor_name(hostname,&hostname_length);
 printf("CPU_NSOC_SECTOR_THREAD_LAYOUT rank=%d hostname=%.*s local_ranks=%d threads=%d actual_team=%d mkl_threads=%d node_cores=64 singleton_threads=1 disjoint_rank_domains=1 cpu_ids=",rank,hostname_length,hostname,local_size,sector_parallel.threads,thread_proof.actual,thread_proof.mkl_threads);
 for(int thread=0;thread<sector_parallel.threads;thread++)printf("%s%d",thread?",":"",thread_proof.cpu[thread]);printf("\n");fflush(stdout);
 MPI_Comm_free(&localcomm);
 if(!S||!H||!W||!B.U||S==H||S==W||S==B.U||H==W||H==B.U||W==B.U)die("SECTOR_PARALLEL requires distinct live S/H/W/U buffers");
 int counts[3],first=0;sector_parallel_ranks(counts);sector_parallel.assigned=-1;
 const int sectors[3]={0,2,3};
 for(int j=0;j<3;j++){
  CpuNSectorGroup*g=sector_parallel.group+j;memset(g,0,sizeof(*g));g->sector=sectors[j];g->p=counts[j];g->pr=sector_parallel_rows(g->p);g->pc=g->p/g->pr;g->root=first;
  g->member=rank>=first&&rank<first+g->p;g->comm=g->rows=g->cols=MPI_COMM_NULL;g->row=g->col=-1;
  int map[256];for(int row=0;row<g->pr;row++)for(int col=0;col<g->pc;col++)map[col*g->pr+row]=first+row*g->pc+col;
  g->ctxt=Csys2blacs_handle(mpi_comm_level1);Cblacs_gridmap(&g->ctxt,map,g->pr,g->pr,g->pc);
  int qr,qc;Cblacs_gridinfo(g->ctxt,&qr,&qc,&g->row,&g->col);
  if(g->member&&(qr!=g->pr||qc!=g->pc||g->row<0||g->col<0))die("SECTOR_PARALLEL actual BLACS membership mismatch");
  if(!g->member&&g->ctxt!=-1)die("SECTOR_PARALLEL inactive BLACS context must be−1");
  int key=g->member?g->col*g->pr+g->row:0;MPI_Comm_split(mpi_comm_level1,g->member?1:MPI_UNDEFINED,key,&g->comm);
  Sector*t=B.sec+g->sector;g->nb=t->nb<64?t->nb:64;
  while(g->nb>1&&(g->nb*g->pr>t->n||g->nb*g->pc>t->n))g->nb/=2;
  g->desc[0]=1;g->desc[1]=-1;g->desc[2]=g->desc[3]=t->n;g->desc[4]=g->desc[5]=g->nb;g->ld=g->desc[8]=1;
  if(g->member){
   if(t->n<g->pr||t->n<g->pc)die("SECTOR_PARALLEL group too large for nonempty ELPA owners");
   int zero=0,info;g->nr=numroc_(&t->n,&g->nb,&g->row,&zero,&g->pr);g->nc=numroc_(&t->n,&g->nb,&g->col,&zero,&g->pc);g->ld=g->nr;
   if(g->nr<1||g->nc<1)die("SECTOR_PARALLEL installed ELPA requires nonempty local rows and columns");
   descinit_(g->desc,&t->n,&t->n,&g->nb,&g->nb,&zero,&zero,&g->ctxt,&g->ld,&info);if(info)die("SECTOR_PARALLEL descriptor initialization");
   if((size_t)g->nr>SIZE_MAX/(size_t)g->nc)die("SECTOR_PARALLEL local size overflow");g->count=(size_t)g->nr*g->nc;
   if(g->count>B.cap||g->count>B.active)die("SECTOR_PARALLEL local owner buffers exceed existing W/U/S capacities before migration");
   MPI_Comm_split(g->comm,g->col,g->row,&g->rows);MPI_Comm_split(g->comm,g->row,g->col,&g->cols);sector_parallel.assigned=j;
   printf("CPU_NSOC_SECTOR_PARALLEL_LAYOUT rank=%d sector=%d group_ranks=%d grid=%d,%d coord=%d,%d block=%d local=%d,%d count=%zu WU_cap=%zu S_cap=%zu extra_dense_matrices=0\n",rank,g->sector,g->p,g->pr,g->pc,g->row,g->col,g->nb,g->nr,g->nc,g->count,B.cap,B.active);fflush(stdout);
  }
  first+=g->p;
 }
 okay=sector_parallel.assigned>=0;MPI_Allreduce(&okay,&all_ok,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!all_ok)die("SECTOR_PARALLEL incomplete membership");
 sector_parallel.n=B.N;sector_parallel.rank=rank;sector_parallel.size=size;sector_parallel.S=S;sector_parallel.H=H;sector_parallel.W=W;sector_parallel.U=B.U;memcpy(sector_parallel.desc,descH,sizeof(sector_parallel.desc));sector_parallel.loaded=1;
}
static void sector_parallel_migrate_factors(double*S,double*W){
 if(!sector_parallel.loaded||sector_parallel.migrated)die("SECTOR_PARALLEL invalid fixed-S migration state");
 int one=1;
 /* Every Ss world-factor source remains intact through all three transfers.
  * A process belongs to one target only, so W holds exactly its own factor. */
 for(int j=0;j<3;j++){
  CpuNSectorGroup*g=sector_parallel.group+j;Sector*t=B.sec+g->sector;
  pdgemr2d_(&t->n,&t->n,S+t->soff,&one,&one,t->desc,W,&one,&one,g->desc,&ictxt1);
 }
 MPI_Barrier(mpi_comm_level1);
 CpuNSectorGroup*g=sector_parallel.group+sector_parallel.assigned;memcpy(S,W,g->count*sizeof(double));sector_parallel.migrated=1;
 if(rank0()){printf("CPU_NSOC_SECTOR_PARALLEL_FACTORS migrated=1 extra_dense_matrices=0 world_sources_consumed_before_overwrite=1\n");fflush(stdout);}
}
static void sector_parallel_input(double*H,double*W,int j,double*destination,int ofs){
 CpuNSectorGroup*g=sector_parallel.group+j;Sector*t=B.sec+g->sector;int one=1,at=ofs+1;
 pdgemr2d_(&t->n,&t->n,H,&at,&at,descH,destination,&one,&one,g->desc,&ictxt1);
 (void)W;
}
static void sector_parallel_solve(int step,double*H,double*S,double*W,double*eval,int phase){
 if(!sector_parallel.loaded||!sector_parallel.migrated)die("SECTOR_PARALLEL missing cached fixed-S factors");
 double begin=MPI_Wtime();
 sector_parallel_input(H,W,0,W,B.sec[0].ofs);sector_parallel_input(H,W,0,B.U,B.sec[1].ofs);
 CpuNSectorGroup*e=sector_parallel.group;double local[2]={0,0},all[2];
 if(e->member)for(size_t i=0;i<e->count;i++){double d=W[i]-B.U[i];local[0]+=d*d;local[1]+=W[i]*W[i]+B.U[i]*B.U[i];}
 MPI_Allreduce(local,all,2,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 if(B.version==3)paired_verify(sqrt(all[0]/fmax(all[1],1e-300)),"D3 Hamiltonian equivalent blocks disagree in sector subgroup");
 if(e->member)for(size_t i=0;i<e->count;i++)W[i]=.5*(W[i]+B.U[i]);
 sector_parallel_input(H,W,1,W,B.sec[2].ofs);sector_parallel_input(H,W,2,W,B.sec[3].ofs);
 MPI_Barrier(mpi_comm_level1);double input_end=MPI_Wtime();
 /* No world collectives occur between this point and the solve completion.
  * W=H, B.U=Q, Ss=fixed factor, each with this rank's private group descriptor. */
 CpuNSectorGroup*g=sector_parallel.group+sector_parallel.assigned;Sector*t=B.sec+g->sector;
 int n=t->n,one=1,nev=n,rows=MPI_Comm_c2f(g->rows),cols=MPI_Comm_c2f(g->cols),comm=MPI_Comm_c2f(g->comm);double a=1,b=0;
 pdtrsm_("L","L","N","N",&n,&n,&a,S,&one,&one,g->desc,W,&one,&one,g->desc);
 pdtrsm_("R","L","T","N",&n,&n,&a,S,&one,&one,g->desc,W,&one,&one,g->desc);
 pdtran_(&n,&n,&a,W,&one,&one,g->desc,&b,B.U,&one,&one,g->desc);
 for(size_t i=0;i<g->count;i++)W[i]=.5*(W[i]+B.U[i]);double transform_end=MPI_Wtime();
 /* ELPA_SCOPE_PATCH_BEGIN */
 CpuNSocElpaOmpScope elpa_omp_scope;
 if(sector_parallel.elpa_omp_threads&&!CpuNSocElpaOmpBegin(&elpa_omp_scope))die("ELPA OMP scope requires actual outsideOMP4/MKL4 and insideOMP1/MKL4");
/* ELPA_SCOPE_PATCH_END */
 elpa_solve_evp_real_2stage_double_impl_(&n,&nev,W,&g->ld,eval+t->ofs,B.U,&g->ld,&g->nb,&g->nc,&rows,&cols,&comm);
/* ELPA_SCOPE_PATCH_BEGIN */
 if(sector_parallel.elpa_omp_threads){
  if(!CpuNSocElpaOmpEnd(&elpa_omp_scope))die("ELPA OMP scope failed restoration to outsideOMP4/MKL4");
  sector_parallel.elpa_omp_calls++;
  if(sector_parallel.rank==g->root){printf("CPU_NSOC_ELPA_OMP_SCOPE step=%d sector=%d rank=%d call=%llu outside_omp=%d inside_omp=%d restored_omp=%d outside_mkl=%d inside_mkl=%d restored_mkl=%d\n",step,g->sector,sector_parallel.rank,sector_parallel.elpa_omp_calls,elpa_omp_scope.before_omp,elpa_omp_scope.inside_omp,elpa_omp_scope.restored_omp,elpa_omp_scope.before_mkl,elpa_omp_scope.inside_mkl,elpa_omp_scope.restored_mkl);fflush(stdout);}
 }
/* ELPA_SCOPE_PATCH_END */

 if(nev!=n)die("SECTOR_PARALLEL ELPA changed full returned state count");double eigen_end=MPI_Wtime();
 pdtrsm_("L","L","T","N",&n,&n,&a,S,&one,&one,g->desc,B.U,&one,&one,g->desc);double recover_end=MPI_Wtime();
 MPI_Barrier(mpi_comm_level1);double solved_end=MPI_Wtime();
 for(int j=0;j<3;j++){
  CpuNSectorGroup*ret=sector_parallel.group+j;Sector*dst=B.sec+ret->sector;int at=dst->ofs+1;
  pdgemr2d_(&dst->n,&dst->n,B.U,&one,&one,ret->desc,H,&at,&at,descH,&ictxt1);
  if(j==0){at=B.sec[1].ofs+1;pdgemr2d_(&dst->n,&dst->n,B.U,&one,&one,ret->desc,H,&at,&at,descH,&ictxt1);}
  MPI_Bcast(eval+dst->ofs,dst->n,MPI_DOUBLE,ret->root,mpi_comm_level1);
  if(j==0)memcpy(eval+B.sec[1].ofs,eval+dst->ofs,dst->n*sizeof(double));
 }
 double end=MPI_Wtime();
 if(phase||step==1){double local_times[7]={input_end-begin,transform_end-input_end,eigen_end-transform_end,recover_end-eigen_end,solved_end-input_end,end-solved_end,end-begin},times[7];MPI_Reduce(local_times,times,7,MPI_DOUBLE,MPI_MAX,0,mpi_comm_level1);if(rank0()){printf("CPU_NSOC_SECTOR_PARALLEL_PHASE step=%d input=%.9f transform_max=%.9f eigen_max=%.9f recover_max=%.9f solve_critical=%.9f return=%.9f total=%.9f full_nev=1 paired_reuse=1\n",step,times[0],times[1],times[2],times[3],times[4],times[5],times[6]);fflush(stdout);}}
}
static void sector_parallel_release(void){
/* ELPA_SCOPE_PATCH_BEGIN */
 if(sector_parallel.elpa_omp_threads){printf("CPU_NSOC_ELPA_OMP_CALLS rank=%d calls=%llu restored_omp=%d restored_mkl=%d\n",sector_parallel.rank,sector_parallel.elpa_omp_calls,omp_get_max_threads(),mkl_get_max_threads());fflush(stdout);}
/* ELPA_SCOPE_PATCH_END */

 if(sector_parallel.loaded)for(int j=0;j<3;j++){CpuNSectorGroup*g=sector_parallel.group+j;if(g->member){MPI_Comm_free(&g->rows);MPI_Comm_free(&g->cols);MPI_Comm_free(&g->comm);Cblacs_gridexit(g->ctxt);}}
 memset(&sector_parallel,0,sizeof(sector_parallel));
}
