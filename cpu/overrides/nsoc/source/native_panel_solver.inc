/* Included after the unchanged sparse congruence/paired solver definitions.
 * The only large output scratch is the existing borrowed B.U allocation. */
#include "rep_panel.h"
#include "panel_budget.h"
#include <errno.h>
#include "native_panel_metadata.h"
typedef struct {int n,nr,nc,desc[9];size_t offset,size;} DirectSector;
static RpaPlan *direct_plan=NULL;static DirectSector direct_sector[3];
static RpaContext direct_context;static double*direct_values=NULL;
static int direct_registered=0;
#include "native_workspace_proof.inc"
static uint64_t direct_hash(uint64_t h,const void*data,size_t bytes){const unsigned char*x=data;for(size_t k=0;k<bytes;k++){h^=x[k];h*=UINT64_C(1099511628211);}return h;}
void CpuNSym_DirectConfigure(void){
 if(direct_mode>=0)return;const char*text=getenv("OPENMX_CPU_NSOC_DIRECT_SECTORS");int mode=0;if(text){char*end;long v=strtol(text,&end,10);if(end==text||*end||v<0||v>2)die("DIRECT_SECTORS must be0,1,2");mode=(int)v;}
 int lo,hi;MPI_Allreduce(&mode,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(&mode,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);if(lo!=hi)die("DIRECT_SECTORS differs across ranks");direct_mode=mode;CpuRSym_PanelEnable(mode!=0);
}
int CpuNSym_DirectRegister(int step,int n,double*****h,double****s){
 CpuNSym_DirectConfigure();direct_active=0;direct_registered=0;if(!direct_mode)return 0;
 if(!CpuNSym_FullProjected()||!CpuNSym_Prepare(step)||SpinP_switch!=0||Cnt_switch||option("OPENMX_CPU_NSOC_PARTIAL_EVP",0)!=0)die("DIRECT_SECTORS requires full projected scalar uncontracted Gamma path");
 if(!CpuRSym_PanelFields(h[0],s))die("DIRECT_SECTORS caller H/S differs from the sealed projected field");
 basis_init(n);if(!d3_basis()||B.version!=3||B.nsec!=4)die("DIRECT_SECTORS currently requires paired real D3 basis version3");
 int nat,ns,ng;const int*types,*norb;const int*const*perms;const double*const*rot;
 if(!CpuRSym_PanelAOMetadata(&nat,&ns,&ng,&types,&norb,&perms,&rot))die("DIRECT_SECTORS missing fresh native Reynolds/AO metadata seal");
 if(!direct_plan){
  uint64_t bh=UINT64_C(1469598103934665603);bh=direct_hash(bh,B.cp,(size_t)(n+1)*8);bh=direct_hash(bh,B.ri,(size_t)B.nnz*4);bh=direct_hash(bh,B.cv,(size_t)B.nnz*8);
  int np;MPI_Comm_size(mpi_comm_level1,&np);uint64_t dh=UINT64_C(1469598103934665603);dh=direct_hash(dh,B.rankmap,(size_t)np*4);int layout[3]={np_rows,np_cols,32};dh=direct_hash(dh,layout,sizeof(layout));direct_context=CpuRSym_PanelContext(bh?bh:1,dh?dh:1);
  RpaBasis basis={0};basis.n=n;basis.natoms=nat;basis.species=ns;basis.ng=ng;for(int k=0;k<4;k++)basis.sectors[k]=B.sz[k];basis.cp=B.cp;basis.ri=B.ri;basis.cv=B.cv;basis.types=types;basis.norb=norb;basis.permutations=perms;basis.rotations=rot;basis.context=direct_context;
  int *coords=mem(2*np,sizeof(int));for(int r=0;r<np_rows;r++)for(int c=0;c<np_cols;c++){int rank=B.rankmap[r*np_cols+c];coords[2*rank]=r;coords[2*rank+1]=c;}
  MPI_Comm local_comm;MPI_Comm_split_type(mpi_comm_level1,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,&local_comm);
  int local_ranks;MPI_Comm_size(local_comm,&local_ranks);size_t headroom=OpenMX_PanelHeadroom();unsigned long long budget=OpenMX_PanelBudget(headroom,local_ranks),node_budget;
  MPI_Allreduce(&budget,&node_budget,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,local_comm);MPI_Comm_free(&local_comm);budget=node_budget;
  const char*override=getenv("OPENMX_CPU_NSOC_PANEL_BUDGET_BYTES");
  if(override){char*end;errno=0;unsigned long long requested=strtoull(override,&end,10);if(errno||end==override||*end||*override=='-'||requested>SIZE_MAX)die("invalid PANEL_BUDGET_BYTES");if(requested<budget)budget=requested;}
  if(rank0()){printf("CPU_NSOC_PANEL_BUDGET bytes_per_rank=%llu local_ranks=%d measured_headroom_bytes=%zu reserve_factor=8 source=node_and_cgroup_available\n",budget,local_ranks,headroom);fflush(stdout);}
  int count;RpaBlock*blocks=CpuRSym_PanelBlocks(&count);char error[256];int failure=RPA_CREATE_OK;
  direct_plan=Rpa_CreateChecked(mpi_comm_level1,&basis,count,blocks,32,np_rows,np_cols,coords,(size_t)budget,error,sizeof(error),&failure);free(blocks);free(coords);
  if(!direct_plan){
   if(failure!=RPA_CREATE_RESOURCE&&failure!=RPA_CREATE_RANGE)die(error);
   if(rank0()){printf("CPU_NSOC_DIRECT_FALLBACK step=%d reason=%s failure_kind=%d path=projected_D3_full_assembly paired_D3_solver=retained symmetry_disabled=0\n",step,error,failure);fflush(stdout);}
   /* Current fields are already Reynolds-projected. Returning0 makes the
    * existing caller assemble them through the certified full-AO D3 path.
    * Keep the seal provider enabled: changing it after its first revision is
    * forbidden, and its future hashes remain valid diagnostics. */
   direct_mode=0;direct_active=0;direct_registered=0;return 0;
  }
  size_t used=0;int index[3]={0,2,3};for(int k=0;k<3;k++){DirectSector*t=direct_sector+k;t->n=B.sec[index[k]].n;int block=32,zero=0,info;t->nr=numroc_(&t->n,&block,&my_prow,&zero,&np_rows);t->nc=numroc_(&t->n,&block,&my_pcol,&zero,&np_cols);int ld=t->nr?t->nr:1;descinit_(t->desc,&t->n,&t->n,&block,&block,&zero,&zero,&ictxt1,&ld,&info);if(info)die("DIRECT_SECTORS scratch descriptor");t->offset=used;t->size=Rpa_LocalSize(direct_plan,k);if(t->size!=(size_t)t->nr*t->nc||used>B.cap||t->size>B.cap-used)die("DIRECT_SECTORS outputs exceed existing eigen scratch");used+=t->size;}
  direct_values=mem(Rpa_InputSize(direct_plan),8);RpaStats stats;Rpa_GetStats(direct_plan,&stats);
  int rank;MPI_Comm_rank(mpi_comm_level1,&rank);printf("CPU_NSOC_DIRECT_PLAN rank=%d local_reps=%zu endpoints=%zu panels=%zu routes=%zu retained_bytes=%zu peak_bytes=%zu send_H_values=%zu recv_H_values=%zu output_existing_scratch=%zu buffer_send_values=%zu buffer_recv_values=%zu value_batches=%d extra_dense_matrices=0 strong_D3_intertwiner=%.12e\n",rank,stats.local_representatives,stats.endpoint_count,stats.panels,stats.routes,stats.retained_bytes,stats.peak_bytes,stats.send_values,stats.recv_values,used,stats.send_capacity_values,stats.recv_capacity_values,stats.value_collectives,stats.intertwiner_residual);fflush(stdout);
 }
 RpaStamp stamp;if(step==1&&!CpuRSym_PanelPack(2,0,direct_context,direct_values,Rpa_InputSize(direct_plan),&stamp))die("DIRECT_SECTORS overlap seal stale before assembly bypass");if(!CpuRSym_PanelPack(1,step,direct_context,direct_values,Rpa_InputSize(direct_plan),&stamp))die("DIRECT_SECTORS H seal stale before assembly bypass");
 direct_registered=1;direct_active=1;direct_step=step;(void)h;(void)s;
 if(rank0()){printf("CPU_NSOC_DIRECT_REGISTER step=%d mode=%d raw_assembly_bypassed=%d certificate=actual_D3_AO_intertwiner_distinct_A_fresh_Reynolds_hash\n",step,direct_mode,direct_mode==1);fflush(stdout);}return direct_mode==1;
}
static void direct_assemble(int step,int kind,double*full,int*desc,double*work){
 if(!direct_active||!direct_registered||step!=direct_step)die("DIRECT_SECTORS missing matching registration");RpaStamp stamp;int field_step=kind==1?step:0;if(!CpuRSym_PanelPack(kind,field_step,direct_context,direct_values,Rpa_InputSize(direct_plan),&stamp))die("DIRECT_SECTORS projected field seal stale at solve");
 double*out[3];size_t capacity[3];for(int k=0;k<3;k++){out[k]=B.U+direct_sector[k].offset;capacity[k]=direct_sector[k].size;}
 char error[256];double start=MPI_Wtime();if(Rpa_Apply(direct_plan,&stamp,direct_values,out,capacity,error,sizeof(error)))die(error);double applied=MPI_Wtime();
 if(direct_mode==2){
  /* Native reference was transformed before the direct result overwrites it.
   * Work is existing Cs scratch; no new complete AO matrix is retained. */
  int one=1,index[3]={0,2,3};double local_error=0,maximum;
  for(int k=0;k<3;k++){DirectSector*t=direct_sector+k;int a=B.sec[index[k]].ofs+1;pdgemr2d_(&t->n,&t->n,full,&a,&a,desc,work,&one,&one,t->desc,&ictxt1);for(size_t i=0;i<t->size;i++){double d=fabs(work[i]-out[k][i]);if(!isfinite(d))die("DIRECT_SECTORS nonfinite native oracle");if(d>local_error)local_error=d;}if(k==0){a=B.sec[1].ofs+1;pdgemr2d_(&t->n,&t->n,full,&a,&a,desc,work,&one,&one,t->desc,&ictxt1);for(size_t i=0;i<t->size;i++){double d=fabs(work[i]-out[k][i]);if(!isfinite(d))die("DIRECT_SECTORS nonfinite E-partner oracle");if(d>local_error)local_error=d;}}}
  MPI_Allreduce(&local_error,&maximum,1,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);if(rank0()){printf("CPU_NSOC_DIRECT_NATIVE_AUDIT step=%d field=%c max_abs=%.12e tolerance=1e-10 all_E_partners=1 source=original_full_AO_assembly_and_sparse_congruence\n",step,kind==1?'H':'S',maximum);fflush(stdout);}if(!isfinite(maximum)||maximum>1e-10)die("DIRECT_SECTORS native full-matrix reference differs");
 }
 int packed=kind==1&&CpuNSym_NativeReturnEnabled();
 int packed_s=kind==2&&direct_mode==1&&CpuNSym_NativeReturnEnabled();
 if(!packed&&!packed_s)memset(full,0,B.active*8);
 int one=1,index[3]={0,2,3};for(int k=0;k<3;k++){
  DirectSector*t=direct_sector+k;Sector*target=B.sec+index[k];int a=target->ofs+1;
  if(packed_s)pdgemr2d_(&t->n,&t->n,out[k],&one,&one,t->desc,full+target->soff,&one,&one,target->desc,&ictxt1);
  else if(packed)pdgemr2d_(&t->n,&t->n,out[k],&one,&one,t->desc,full+B.hoff[index[k]],&one,&one,target->desc,&ictxt1);
  else pdgemr2d_(&t->n,&t->n,out[k],&one,&one,t->desc,full,&a,&a,desc,&ictxt1);
  if(k==0&&!packed_s){if(packed)memcpy(full+B.hoff[1],full+B.hoff[0],(size_t)target->nr*target->nc*sizeof(double));else{a=B.sec[1].ofs+1;pdgemr2d_(&t->n,&t->n,out[k],&one,&one,t->desc,full,&a,&a,desc,&ictxt1);}}
 }
 if(rank0()){printf("CPU_NSOC_DIRECT_ASSEMBLY step=%d field=%c apply_rank0=%.9f total_rank0=%.9f payload=H_or_S_only symmetry_cross_blocks=construction_certificate extra_dense_matrices=0\n",step,kind==1?'H':'S',applied-start,MPI_Wtime()-start);fflush(stdout);}
}
static void direct_release(void){Rpa_Destroy(direct_plan);direct_plan=NULL;free(direct_values);direct_values=NULL;direct_active=0;direct_registered=0;memset(direct_workspace,0,sizeof(direct_workspace));memset(direct_workspace_count,0,sizeof(direct_workspace_count));}
