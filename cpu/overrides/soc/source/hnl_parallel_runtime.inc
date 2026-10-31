typedef struct {int enabled,count;size_t extra_bytes;CpuHNLWorkspace serial,*w;} CpuHNLPool;
static int CpuHNLParallelMode(void){
 static int cached=-1;if(cached>=0)return cached;const char*s=getenv("OPENMX_CPU_HNL_REMOTE_PARALLEL");int mode=!s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:-1,in[2]={mode,-mode},out[2];
 MPI_Allreduce(in,out,2,MPI_INT,MPI_MIN,mpi_comm_level1);if(out[0]<0||out[0]!=-out[1])CpuF3TileFail("HNL_REMOTE_PARALLEL must agree across ranks and be exactly0or1");
 if(mode){int okay=!omp_in_parallel()&&MD_switch==0&&Solver==2&&SpinP_switch==3&&SO_switch==1&&Cnt_switch==0&&DFTSCF_loop==1&&Hub_U_switch==0&&Constraint_NCS_switch==0&&Zeeman_NCS_switch==0&&Zeeman_NCO_switch==0,all;
  MPI_Allreduce(&okay,&all,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!all)CpuF3TileFail("HNL_REMOTE_PARALLEL requires serial-entry fixed one-step uncontracted nonmagnetic Cluster SOC");
 }
 cached=mode;return mode;
}
static void CpuHNLPoolInit(CpuHNLPool*p,dcomplex***x0,dcomplex***y0,dcomplex***z0,dcomplex***x1,dcomplex***y1,dcomplex***z1){
 memset(p,0,sizeof(*p));p->enabled=CpuHNLParallelMode();p->count=1;p->w=&p->serial;
 dcomplex***view[6]={x0,y0,z0,x1,y1,z1};for(int k=0;k<6;k++)p->serial.m[k]=view[k];
 if(p->enabled&&Matomnum){
  int nt=omp_get_max_threads(),width=List_YOUSO[7];if(nt<1||nt>64||width<1||width>128)CpuF3TileFail("HNL private workspace width/thread count outside bounded scope");
  size_t per=(size_t)18*width*width*sizeof(dcomplex)+(size_t)18*width*sizeof(dcomplex*)+18*sizeof(dcomplex**);
  size_t table=(size_t)nt*sizeof(*p->w),limit=(size_t)64*1024*1024;
  if(table>limit||(size_t)(nt-1)>(limit-table)/per)CpuF3TileFail("HNL private workspaces exceed64MiB rank bound");
  for(int m=1;m<=Matomnum;m++){int g=M2G[m];if(g<1||g>atomnum)CpuF3TileFail("HNL local atom map outside global range");for(int q=1;q<m;q++)if(M2G[q]==g)CpuF3TileFail("HNL local atom map not injective");}
  p->w=calloc((size_t)nt,sizeof(*p->w));if(!p->w)CpuF3TileFail("HNL workspace table allocation");p->count=nt;p->extra_bytes=table;p->w[0]=p->serial;
  for(int t=1;t<nt;t++){if(CpuHNLWorkspaceAllocate(p->w+t,width))CpuF3TileFail("HNL private workspace allocation");p->extra_bytes+=p->w[t].bytes;}
 }
}
#include "hnl_remote_atom.inc"
static void CpuHNLRemoteRun(CpuHNLPool*p,int IDR,int Original_Mc_AN,double*****CDM0,double*****iDM0){
 if(p->enabled&&Matomnum){
#pragma omp parallel for schedule(dynamic,1)
  for(int m=1;m<=Matomnum;m++){int thread=omp_get_thread_num();if(thread<0||thread>=p->count)abort();CpuHNLRemoteAtom(m,IDR,Original_Mc_AN,CDM0,iDM0,p->w+thread);}
  /* omp parallel-for join completes every read before ring bookkeeping and
   * the next receive overwrite DS_NL[...][Matomnum+1]. */
 }else for(int m=1;m<=Matomnum;m++)CpuHNLRemoteAtom(m,IDR,Original_Mc_AN,CDM0,iDM0,p->w);
}
static void CpuHNLPoolFree(CpuHNLPool*p){
 if(p->enabled||CpuFP.active){
  int rank;MPI_Comm_rank(mpi_comm_level1,&rank);char line[1024];int used=snprintf(line,sizeof(line),"CPU_HNL_REMOTE_THREADS rank=%d enabled=%d threads=%d extra_workspace_bytes=%zu",rank,p->enabled,p->count,p->extra_bytes);
  if(used<0||(size_t)used>=sizeof(line))abort();
  if(p->count<=4){for(int t=0;t<p->count;t++){int n=snprintf(line+used,sizeof(line)-(size_t)used," t%d_atoms=%llu t%d_pairs=%llu t%d_dhnl_calls=%llu",t,p->w[t].atoms,t,p->w[t].pairs,t,p->w[t].calls);if(n<0||(size_t)n>=sizeof(line)-(size_t)used)abort();used+=n;}}
  if((size_t)used+1>=sizeof(line))abort();line[used++]='\n';fflush(stdout);if(fwrite(line,1,(size_t)used,stdout)!=(size_t)used)abort();fflush(stdout);
  if(p->count>4)for(int t=0;t<p->count;t++){int n=snprintf(line,sizeof(line),"CPU_HNL_REMOTE_THREAD rank=%d thread=%d atoms=%llu pairs=%llu dhnl_calls=%llu\n",rank,t,p->w[t].atoms,p->w[t].pairs,p->w[t].calls);if(n<0||(size_t)n>=sizeof(line))abort();if(fwrite(line,1,(size_t)n,stdout)!=(size_t)n)abort();fflush(stdout);}
 }
 for(int t=1;t<p->count;t++)CpuHNLWorkspaceFree(p->w+t);if(p->w!=&p->serial)free(p->w);memset(p,0,sizeof(*p));
}
