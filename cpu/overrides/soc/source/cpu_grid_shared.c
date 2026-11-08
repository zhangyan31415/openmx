#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "openmx_common.h"
#include "cpu_grid_shared.h"
#include <errno.h>
#include <sched.h>
#include <sys/syscall.h>
#include <sys/mman.h>
#include <unistd.h>

static MPI_Comm node_comm=MPI_COMM_NULL;
static MPI_Comm grid_world=MPI_COMM_NULL; /* borrowed: native calculation owns it */
static MPI_Win grid_win=MPI_WIN_NULL;
static unsigned char *grid_base=NULL;
static size_t *atom_offsets=NULL;
static int *atom_grids=NULL, *atom_norbs=NULL;
static int atom_count=0;
static size_t grid_bytes=0;
static MPI_Comm numa_comm=MPI_COMM_NULL;
static MPI_Win replica_win=MPI_WIN_NULL;
static unsigned char *replica_base=NULL;
static int numa_enabled=0,local_numa=-1,numa_rank=0,numa_size=0,numa_domains=0;
static int node_leader=-1,world_nodes=0;
static int *atom_owner_numa=NULL;
static size_t *replica_offsets=NULL;
static size_t legacy_remote_copy_bytes=0;
static int replica_atoms=0;
static size_t replica_mapped=0,replica_copy_bytes=0;
static unsigned long long publish_epoch=0;
static double last_copy_rankmax=0.0,last_publish_rankmax=0.0;
static int grid_ready=0,grid_suspended=0,grid_rebuilding=0,last_suspend_removed=0;
static unsigned long long suspend_epoch=0;

static void fail(MPI_Comm comm, const char *message)
{
  int rank=0;
  MPI_Comm_rank(comm,&rank);
  if(rank==0){ fprintf(stderr,"[CPU_GRID_SHARED] %s\n",message); fflush(stderr); }
  MPI_Abort(comm,1);
  abort();
}

static void collective_require(MPI_Comm comm, int condition, const char *message)
{
  int all_ok=0;
  MPI_Allreduce(&condition,&all_ok,1,MPI_INT,MPI_MIN,comm);
  if(!all_ok) fail(comm,message);
}

static int strict_flag(const char *name)
{
  const char *value=getenv(name);int initialized=0;
  if(!value||!strcmp(value,"0"))return 0;
  if(!strcmp(value,"1"))return 1;
  MPI_Initialized(&initialized);
  if(initialized)fail(MPI_COMM_WORLD,"grid shared/NUMA flags must be0or1");
  fprintf(stderr,"%s must be0or1\n",name);exit(EXIT_FAILURE);return 0;
}
int CpuGridShared_Requested(void)
{
  int shared=strict_flag("OPENMX_CPU_GRID_SHARED");
  int numa=strict_flag("OPENMX_CPU_GRID_NUMA"),initialized=0;
  /* No collective here: Requested may be called differing numbers of times. */
  if(numa&&!shared){
    MPI_Initialized(&initialized);
    if(initialized)fail(MPI_COMM_WORLD,"OPENMX_CPU_GRID_NUMA=1 requires OPENMX_CPU_GRID_SHARED=1");
    fprintf(stderr,"OPENMX_CPU_GRID_NUMA=1 requires OPENMX_CPU_GRID_SHARED=1\n");exit(EXIT_FAILURE);
  }
  return shared;
}

int CpuGridShared_PlanBytes(size_t ngrid,size_t norb,size_t previous,size_t *next)
{
  size_t values, bytes;
  if(next==NULL || (norb && ngrid>SIZE_MAX/norb)) return 0;
  values=ngrid*norb;
  if(values>SIZE_MAX/sizeof(Type_Orbs_Grid)) return 0;
  bytes=values*sizeof(Type_Orbs_Grid);
  if(previous>SIZE_MAX-bytes) return 0;
  *next=previous+bytes;
  /* MPI_Aint and pointer differences must both represent the allocation. */
  if(*next>(size_t)PTRDIFF_MAX || (MPI_Aint)*next<0 || (size_t)(MPI_Aint)*next!=*next) return 0;
  return 1;
}

static void init_numa(MPI_Comm comm,int natoms,const int *owners,int node_size);

void CpuGridShared_Init(MPI_Comm comm,int natoms,const int *owners,
                        const int *grid_counts,const int *species,
                        const int *species_norb)
{
  int rank,nrank,node_rank,node_size,*members=NULL,*on_node=NULL;
  int atom,valid=1,rc,allocation_ok;
  MPI_Aint query_bytes=0;
  int disp_unit=0;
  void *segment=NULL;
  size_t bytes=0, length;
  if(node_comm!=MPI_COMM_NULL) CpuGridShared_Destroy();
  grid_world=comm;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&nrank);
  const char *domain=getenv("OPENMX_CPU_GRID_DOMAIN");
  collective_require(comm,!domain||!strcmp(domain,"node"),"replicated NUMA variant requires canonical GRID_DOMAIN=node");
  collective_require(comm,natoms>=0 && sizeof(Type_Orbs_Grid)==4 && owners && grid_counts && species && species_norb,
                     "invalid atom metadata or unsupported orbital grid scalar type");
  length=(size_t)natoms+1;
  collective_require(comm,length<=SIZE_MAX/sizeof(size_t),"atom metadata size overflow");
  rc=MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node_comm);
  collective_require(comm,rc==MPI_SUCCESS && node_comm!=MPI_COMM_NULL,"shared-memory communicator creation failed");
  MPI_Comm_set_errhandler(node_comm,MPI_ERRORS_RETURN);
  MPI_Comm_rank(node_comm,&node_rank);
  MPI_Comm_size(node_comm,&node_size);
  members=(int*)malloc((size_t)node_size*sizeof(int));
  on_node=(int*)calloc((size_t)nrank,sizeof(int));
  atom_offsets=(size_t*)malloc(length*sizeof(size_t));
  atom_grids=(int*)malloc(length*sizeof(int));
  atom_norbs=(int*)malloc(length*sizeof(int));
  collective_require(comm,members && on_node && atom_offsets && atom_grids && atom_norbs,
                     "allocation of atom metadata failed");
  MPI_Allgather(&rank,1,MPI_INT,members,1,MPI_INT,node_comm);
  for(atom=0;atom<node_size;atom++) on_node[members[atom]]=1;
  for(atom=0;atom<=natoms;atom++){
    atom_offsets[atom]=SIZE_MAX;
    atom_grids[atom]=atom_norbs[atom]=0;
  }
  for(atom=1;atom<=natoms;atom++){
    int norb;
    if(owners[atom]<0 || owners[atom]>=nrank || grid_counts[atom]<0 || species[atom]<0){ valid=0; continue; }
    norb=species_norb[species[atom]];
    if(norb<0){ valid=0; continue; }
    if(!on_node[owners[atom]]) continue;
    atom_offsets[atom]=bytes;
    atom_grids[atom]=grid_counts[atom];
    atom_norbs[atom]=norb;
    if(!CpuGridShared_PlanBytes((size_t)grid_counts[atom],(size_t)norb,bytes,&bytes)) valid=0;
  }
  free(members); free(on_node);
  collective_require(comm,valid,"invalid metadata or shared orbital cache byte-size overflow");
  /* Exactly one physical cache per shared-memory domain, with owner ranks
     writing disjoint atom ranges; other ranks allocate no payload. */
  rc=MPI_Win_allocate_shared(node_rank==0?(MPI_Aint)bytes:0,1,MPI_INFO_NULL,node_comm,&segment,&grid_win);
  allocation_ok=(rc==MPI_SUCCESS && grid_win!=MPI_WIN_NULL);
  collective_require(comm,allocation_ok,"MPI shared orbital cache allocation failed");
  MPI_Win_set_errhandler(grid_win,MPI_ERRORS_RETURN);
  rc=MPI_Win_shared_query(grid_win,0,&query_bytes,&disp_unit,&segment);
  collective_require(comm,rc==MPI_SUCCESS && query_bytes==(MPI_Aint)bytes && (!bytes || segment!=NULL),
                     "MPI shared orbital cache query failed");
  grid_base=(unsigned char*)segment;
  grid_bytes=bytes;
  atom_count=natoms;
  rc=MPI_Win_lock_all(0,grid_win);
  collective_require(comm,rc==MPI_SUCCESS,"shared orbital cache access epoch failed");
  init_numa(comm,natoms,owners,node_size);
  if(node_rank==0){
    printf("[CPU_GRID_SHARED] node_leader_rank=%d node_ranks=%d atoms=%d unique_payload_bytes=%zu metadata_bytes_per_rank=%zu physical_payload_copies=1 scope=canonical_only numa_replicas_reported_separately=1\n",
           rank,node_size,natoms,grid_bytes,length*(sizeof(size_t)+2*sizeof(int)));
    fflush(stdout);
  }
}

int CpuGridShared_IsAtom(int global_atom)
{
  return node_comm!=MPI_COMM_NULL && global_atom>0 && global_atom<=atom_count && atom_offsets[global_atom]!=SIZE_MAX;
}

Type_Orbs_Grid *CpuGridShared_Row(int global_atom,int grid_index)
{
  if(!CpuGridShared_IsAtom(global_atom) || grid_index<0 || grid_index>=atom_grids[global_atom] || atom_norbs[global_atom]==0) return NULL;
  return (Type_Orbs_Grid*)(grid_base+atom_offsets[global_atom]+(size_t)grid_index*(size_t)atom_norbs[global_atom]*sizeof(Type_Orbs_Grid));
}

void CpuGridShared_Sync(void)
{
  if(grid_win!=MPI_WIN_NULL){
    int rc=MPI_Win_sync(grid_win);
    if(rc!=MPI_SUCCESS) fail(node_comm,"shared orbital cache writer synchronization failed");
    MPI_Barrier(node_comm);
    rc=MPI_Win_sync(grid_win);
    if(rc!=MPI_SUCCESS) fail(node_comm,"shared orbital cache reader synchronization failed");
  }
}

size_t CpuGridShared_Bytes(void) { return grid_bytes; }

static size_t atom_bytes(int atom)
{
  return (size_t)atom_grids[atom]*(size_t)atom_norbs[atom]*sizeof(Type_Orbs_Grid);
}
static void init_numa(MPI_Comm comm,int natoms,const int *owners,int node_size)
{
  int flag=strict_flag("OPENMX_CPU_GRID_NUMA"),low,high,rank,nrank,i,j,valid=1,rc=MPI_SUCCESS;
  int *rank_domain=NULL,identity[2];unsigned cpu=0,node=0;cpu_set_t affinity;
  MPI_Aint queried=0;int disp=0;void *segment=NULL;
  MPI_Allreduce(&flag,&low,1,MPI_INT,MPI_MIN,comm);
  MPI_Allreduce(&flag,&high,1,MPI_INT,MPI_MAX,comm);
  collective_require(comm,low==high,"OPENMX_CPU_GRID_NUMA must agree on every rank");
  if(!flag)return;
  collective_require(comm,strict_flag("OPENMX_CPU_GRID_SHARED")==1,"NUMA grid requires shared grid mode at Init");
  MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&nrank);
#if defined(__linux__) && defined(SYS_getcpu)
  CPU_ZERO(&affinity);
  if(sched_getaffinity(0,sizeof(affinity),&affinity)!=0||CPU_COUNT(&affinity)!=1)valid=0;
  if(syscall(SYS_getcpu,&cpu,&node,NULL)!=0||cpu>=CPU_SETSIZE||node>INT_MAX||!CPU_ISSET(cpu,&affinity))valid=0;
#else
  valid=0;
#endif
  collective_require(comm,valid,"NUMA grid requires reliable Linux getcpu and one-CPU rank affinity");
  local_numa=(int)node;
  node_leader=rank;
  MPI_Allreduce(MPI_IN_PLACE,&node_leader,1,MPI_INT,MPI_MIN,node_comm);
  rank_domain=malloc((size_t)nrank*2*sizeof(*rank_domain));
  atom_owner_numa=malloc(((size_t)natoms+1)*sizeof(*atom_owner_numa));
  replica_offsets=malloc(((size_t)natoms+1)*sizeof(*replica_offsets));
  collective_require(comm,rank_domain&&atom_owner_numa&&replica_offsets,"NUMA grid owner metadata allocation failed");
  /* Linux NUMA IDs are host-local: the physical-node identity must participate
     before comparing them. No address or row for another host is exposed. */
  identity[0]=node_leader;identity[1]=local_numa;
  MPI_Allgather(identity,2,MPI_INT,rank_domain,2,MPI_INT,comm);
  numa_domains=0;world_nodes=0;
  for(i=0;i<nrank;i++){
    if(rank_domain[2*i]==i)world_nodes++;
    if(rank_domain[2*i]!=node_leader)continue;
    for(j=0;j<i;j++)if(rank_domain[2*j]==node_leader&&rank_domain[2*j+1]==rank_domain[2*i+1])break;
    if(j==i)numa_domains++;
  }
  atom_owner_numa[0]=-1;replica_copy_bytes=legacy_remote_copy_bytes=0;replica_atoms=0;valid=1;j=0;
  unsigned char *needed=calloc((size_t)natoms+1,1);
  collective_require(comm,needed&&FNAN&&natn,"NUMA grid demand graph allocation or neighbor metadata missing");
  for(i=0;i<=natoms;i++)replica_offsets[i]=SIZE_MAX;
  for(i=0;i<nrank;i++)if(rank_domain[2*i]==node_leader)j++;
  if(j!=node_size)valid=0;
  /* ReadRow consumers bind first-neighbor orbital rows for their owned
     centers. Build the union for this physical node + NUMA identity before
     freeing the rank-domain map. No off-node pointer is exposed. */
  for(i=1;i<=natoms;i++){
    atom_owner_numa[i]=-1;
    if(!CpuGridShared_IsAtom(i))continue;
    if(rank_domain[2*owners[i]]!=node_leader){valid=0;continue;}
    atom_owner_numa[i]=rank_domain[2*owners[i]+1];
    if(atom_owner_numa[i]!=local_numa)continue;
    if(FNAN[i]<0||!natn[i]){valid=0;continue;}
    for(int h=0;h<=FNAN[i];h++){
      int atom=natn[i][h];if(atom<1||atom>natoms){valid=0;continue;}
      if(CpuGridShared_IsAtom(atom))needed[atom]=1;
    }
  }
  for(i=1;i<=natoms;i++){
    if(!CpuGridShared_IsAtom(i)||atom_owner_numa[i]==local_numa)continue;
    legacy_remote_copy_bytes+=atom_bytes(i);
    if(!needed[i])continue;
    size_t next;
    if(!CpuGridShared_PlanBytes((size_t)atom_grids[i],(size_t)atom_norbs[i],replica_copy_bytes,&next)){valid=0;continue;}
    replica_offsets[i]=replica_copy_bytes;replica_copy_bytes=next;replica_atoms++;
  }
  free(needed);free(rank_domain);
  collective_require(comm,valid,"invalid node/NUMA neighbor demand or replica byte-size overflow");
  rc=MPI_Comm_split(node_comm,local_numa,rank,&numa_comm);
  collective_require(comm,rc==MPI_SUCCESS&&numa_comm!=MPI_COMM_NULL,"NUMA communicator creation failed");
  MPI_Comm_set_errhandler(numa_comm,MPI_ERRORS_RETURN);
  MPI_Comm_rank(numa_comm,&numa_rank);MPI_Comm_size(numa_comm,&numa_size);
  replica_mapped=replica_copy_bytes; /* packed demand, without canonical/unused holes */
  valid=1;
  if(replica_mapped){
    rc=MPI_Win_allocate_shared(numa_rank==0?(MPI_Aint)replica_mapped:0,1,MPI_INFO_NULL,numa_comm,&segment,&replica_win);
    valid=rc==MPI_SUCCESS&&replica_win!=MPI_WIN_NULL;
  }
  collective_require(comm,valid,"NUMA replica window allocation failed");
  if(replica_win!=MPI_WIN_NULL){
    MPI_Win_set_errhandler(replica_win,MPI_ERRORS_RETURN);
    rc=MPI_Win_shared_query(replica_win,0,&queried,&disp,&segment);
    valid=rc==MPI_SUCCESS&&queried==(MPI_Aint)replica_mapped&&segment!=NULL;
  }
  collective_require(comm,valid,"NUMA replica window query failed");
  replica_base=(unsigned char*)segment;
  if(replica_win!=MPI_WIN_NULL)rc=MPI_Win_lock_all(0,replica_win);
  collective_require(comm,rc==MPI_SUCCESS,"NUMA replica access epoch failed");
  numa_enabled=1;publish_epoch=0;
  if(numa_rank==0){
    printf("[CPU_GRID_NUMA] active=1 numa_id=%d domains=%d numa_ranks=%d canonical_bytes=%zu replica_mapped_capacity_bytes=%zu replica_copy_payload_bytes=%zu locality_verified=0 node_leader_rank=%d physical_nodes=%d scope=node_owned remote_ghosts=native_evaluation replica_cached_atoms=%d all_remote_copy_bytes=%zu policy=first_neighbor_demand packing=compact\n",
      local_numa,numa_domains,numa_size,grid_bytes,replica_mapped,replica_copy_bytes,node_leader,world_nodes,replica_atoms,legacy_remote_copy_bytes);fflush(stdout);
  }
}
Type_Orbs_Grid *CpuGridShared_ReadRow(int atom,int grid)
{
  if(grid_suspended||grid_rebuilding)return NULL;
  Type_Orbs_Grid *canonical=CpuGridShared_Row(atom,grid);
  if(!canonical||!numa_enabled||atom_owner_numa[atom]==local_numa)return canonical;
  /* Unanticipated consumers retain a correct canonical row. Replication is
     a locality optimization, not a condition for numerical correctness. */
  if(!replica_offsets||replica_offsets[atom]==SIZE_MAX)return canonical;
  if(!replica_base||publish_epoch==0)return NULL;
  return (Type_Orbs_Grid*)(replica_base+replica_offsets[atom]+(size_t)grid*(size_t)atom_norbs[atom]*sizeof(Type_Orbs_Grid));
}
static void replica_sync(void)
{
  if(replica_win!=MPI_WIN_NULL){
    int rc=MPI_Win_sync(replica_win);
    if(rc!=MPI_SUCCESS)fail(node_comm,"NUMA replica writer synchronization failed");
    MPI_Barrier(numa_comm);
    rc=MPI_Win_sync(replica_win);
    if(rc!=MPI_SUCCESS)fail(node_comm,"NUMA replica reader synchronization failed");
  }
}
/* The kernel reports whole VMAs, not arbitrary subranges. Only call this a
   replica locality measurement when the VMA matches this window to page slack. */
static void report_numa_pages(void)
{
  FILE *f;char line[8192],*p;unsigned long lo=0,hi=0,a,b,pages=0,local=0,value;
  unsigned long wanted=(unsigned long)(uintptr_t)replica_base;int id,matched=0,dedicated=0;long page=sysconf(_SC_PAGESIZE);
  if(numa_rank!=0||!replica_base||!replica_mapped)return;
  f=fopen("/proc/self/maps","r");
  if(f){while(fgets(line,sizeof(line),f))if(sscanf(line,"%lx-%lx",&a,&b)==2&&a<=wanted&&wanted<b){lo=a;hi=b;matched=1;break;}fclose(f);}
  if(!matched||page<=0||wanted>ULONG_MAX-replica_mapped||hi<wanted+replica_mapped){
    printf("[CPU_GRID_NUMA] locality numa_id=%d verified=0 reason=unmatched_replica_vma\n",local_numa);return;
  }
  dedicated=(wanted-lo<=(unsigned long)page*2&&hi-(wanted+replica_mapped)<=(unsigned long)page*2);
  f=fopen("/proc/self/numa_maps","r");matched=0;
  if(f){while(fgets(line,sizeof(line),f)){
    if(sscanf(line,"%lx",&a)!=1||a!=lo)continue;
    matched=1;p=line;
    while(*p){if(*p=='N'&&(p==line||p[-1]==' '||p[-1]=='\t')&&sscanf(p,"N%d=%lu",&id,&value)==2){pages+=value;if(id==local_numa)local+=value;}p++;}
    break;
  }fclose(f);}
  printf("[CPU_GRID_NUMA] locality numa_id=%d verified=%d scope=whole_vma vma_begin=0x%lx vma_end=0x%lx reported_pages=%lu local_pages=%lu vma_local_page_fraction=%.9f replica_mapped_capacity_bytes=%zu reason=%s\n",
    local_numa,matched&&dedicated&&pages>0,lo,hi,pages,local,pages?(double)local/pages:-1.,replica_mapped,
    !matched?"numa_maps_unavailable":!dedicated?"vma_contains_unattributed_storage":!pages?"no_resident_pages":"dedicated_replica_vma");fflush(stdout);
}
void CpuGridShared_Publish(void)
{
  double started,copy_started,copy_seconds,publish_seconds;size_t copied=0,ordinal=0;int atom,node_rank;
  unsigned long long local_bytes,all_bytes,mapped,all_mapped;
  if(node_comm==MPI_COMM_NULL)return;
  collective_require(grid_world,!grid_suspended,"Publish requires BeginRebuild after Suspend");
  if(!numa_enabled){CpuGridShared_Sync();grid_ready=1;grid_rebuilding=0;return;}
  started=MPI_Wtime();CpuGridShared_Sync();copy_started=MPI_Wtime();
  for(atom=1;atom<=atom_count;atom++){
    if(!CpuGridShared_IsAtom(atom))continue;
    if(!replica_offsets||replica_offsets[atom]==SIZE_MAX)continue;
    if(ordinal++%(size_t)numa_size!=(size_t)numa_rank)continue;
    size_t bytes=atom_bytes(atom);
    if(bytes){memcpy(replica_base+replica_offsets[atom],grid_base+atom_offsets[atom],bytes);copied+=bytes;}
  }
  copy_seconds=MPI_Wtime()-copy_started;replica_sync();publish_epoch++;
  publish_seconds=MPI_Wtime()-started;
  MPI_Allreduce(&copy_seconds,&last_copy_rankmax,1,MPI_DOUBLE,MPI_MAX,node_comm);
  MPI_Allreduce(&publish_seconds,&last_publish_rankmax,1,MPI_DOUBLE,MPI_MAX,node_comm);
  local_bytes=(unsigned long long)copied;mapped=numa_rank==0?(unsigned long long)replica_mapped:0;
  MPI_Reduce(&local_bytes,&all_bytes,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,node_comm);
  MPI_Reduce(&mapped,&all_mapped,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,node_comm);
  MPI_Comm_rank(node_comm,&node_rank);
  if(node_rank==0){
    printf("[CPU_GRID_NUMA] publish epoch=%llu copied_payload_bytes=%llu replica_mapped_capacity_sum_bytes=%llu copy_wall_rankmax_seconds=%.9f publish_wall_rankmax_seconds=%.9f includes_stats_collectives=0\n",
      publish_epoch,all_bytes,all_mapped,last_copy_rankmax,last_publish_rankmax);fflush(stdout);
  }
  report_numa_pages();
  grid_ready=1;grid_rebuilding=0;
}

/* MADV_DONTNEED would leave the shared backing store resident. MADV_REMOVE
 * punches the backing-store hole and makes all aliases observe zero pages.
 * Keep allocator/MPI boundaries intact by discarding only whole interior
 * pages. Failure leaves that mapping intact; it never changes the algorithm. */
static size_t remove_shared_payload(unsigned char *base,size_t bytes,
                                   size_t *attempted,int *error)
{
  *attempted=0;*error=0;
  if(!base||!bytes)return 0;
  long page=sysconf(_SC_PAGESIZE);
  if(page<=0){*error=EINVAL;return 0;}
  uintptr_t begin=(uintptr_t)base,end;
  if(bytes>UINTPTR_MAX-begin || begin>UINTPTR_MAX-(uintptr_t)(page-1)){
    *error=EOVERFLOW;return 0;
  }
  end=begin+bytes;
  uintptr_t lo=((begin+(uintptr_t)page-1)/(uintptr_t)page)*(uintptr_t)page;
  uintptr_t hi=(end/(uintptr_t)page)*(uintptr_t)page;
  if(hi<=lo)return 0;
  *attempted=(size_t)(hi-lo);
#if defined(CPU_GRID_TEST_REMOVE_UNSUPPORTED)
  *error=ENOSYS;
#elif defined(__linux__) && defined(MADV_REMOVE)
  if(madvise((void*)lo,(size_t)(hi-lo),MADV_REMOVE)==0)return (size_t)(hi-lo);
  *error=errno;
#else
  *error=ENOSYS;
#endif
  return 0;
}
int CpuGridShared_IsSuspended(void){return grid_suspended;}
int CpuGridShared_Suspend(void)
{
  if(node_comm==MPI_COMM_NULL)return 0;
  int state[2]={grid_suspended,grid_ready&&!grid_rebuilding},lo[2],hi[2];
  MPI_Allreduce(state,lo,2,MPI_INT,MPI_MIN,grid_world);
  MPI_Allreduce(state,hi,2,MPI_INT,MPI_MAX,grid_world);
  collective_require(grid_world,lo[0]==hi[0]&&lo[1]==hi[1],"Suspend state differs across ranks");
  if(grid_suspended)return last_suspend_removed;
  collective_require(grid_world,grid_ready&&!grid_rebuilding,"Suspend requires published stable grids");
  CpuGridShared_Sync();replica_sync();MPI_Barrier(grid_world);
  int rank,node_rank,error=0,error_max;
  MPI_Comm_rank(grid_world,&rank);MPI_Comm_rank(node_comm,&node_rank);
  unsigned long long local[6]={0,0,0,0,0,0},sum[6];
  if(node_rank==0){
    size_t tried=0;int why=0;size_t removed=remove_shared_payload(grid_base,grid_bytes,&tried,&why);
    local[0]=tried;local[1]=removed;local[4]+=why!=0;local[5]+=tried!=0;
    if(why>error)error=why;
  }
  if(numa_enabled&&numa_rank==0){
    size_t tried=0;int why=0;size_t removed=remove_shared_payload(replica_base,replica_mapped,&tried,&why);
    local[2]=tried;local[3]=removed;local[4]+=why!=0;local[5]+=tried!=0;
    if(why>error)error=why;
  }
  CpuGridShared_Sync();replica_sync();MPI_Barrier(grid_world);
  MPI_Allreduce(local,sum,6,MPI_UNSIGNED_LONG_LONG,MPI_SUM,grid_world);
  MPI_Allreduce(&error,&error_max,1,MPI_INT,MPI_MAX,grid_world);
  grid_ready=0;grid_suspended=1;grid_rebuilding=0;suspend_epoch++;
  last_suspend_removed=(sum[1]||sum[3])?1:0;
  if(rank==0){
    printf("[CPU_GRID_SUSPEND] epoch=%llu canonical_attempted_bytes=%llu canonical_removed_bytes=%llu replica_attempted_bytes=%llu replica_removed_bytes=%llu windows=%llu failed_windows=%llu errno_max=%d result=%s method=MADV_REMOVE virtual_capacity_retained=1 native_rebuild_required=1\n",
      suspend_epoch,sum[0],sum[1],sum[2],sum[3],sum[5],sum[4],error_max,
      !last_suspend_removed?"safe_noop":sum[4]?"partial":"removed");fflush(stdout);
  }
  return last_suspend_removed;
}
void CpuGridShared_BeginRebuild(void)
{
  if(node_comm==MPI_COMM_NULL)return;
  collective_require(grid_world,grid_suspended,"BeginRebuild requires suspended grids");
  /* No reader can become ready until all canonical owner writes and replica
   * copies have passed the ordinary native publication barriers. */
  MPI_Barrier(grid_world);grid_suspended=0;grid_rebuilding=1;grid_ready=0;
}
void CpuGridShared_NumaStats(CpuGridSharedNumaStats *out)
{
  if(!out)return;out->enabled=numa_enabled;out->local_numa=local_numa;out->domains=numa_domains;
  out->mapped_bytes=replica_mapped;out->copy_bytes=replica_copy_bytes;out->epochs=publish_epoch;
  out->copy_rankmax_seconds=last_copy_rankmax;out->publish_rankmax_seconds=last_publish_rankmax;
}
void CpuGridShared_Destroy(void)
{
  /* All ranks finish every reader before either canonical or replica unmaps. */
  if(node_comm!=MPI_COMM_NULL)CpuGridShared_Sync();
  if(replica_win!=MPI_WIN_NULL){replica_sync();MPI_Win_unlock_all(replica_win);MPI_Win_free(&replica_win);}
  if(numa_comm!=MPI_COMM_NULL)MPI_Comm_free(&numa_comm);
  if(grid_win!=MPI_WIN_NULL){MPI_Win_unlock_all(grid_win);MPI_Win_free(&grid_win);}
  if(node_comm!=MPI_COMM_NULL)MPI_Comm_free(&node_comm);
  free(atom_offsets);free(atom_grids);free(atom_norbs);free(atom_owner_numa);free(replica_offsets);
  atom_offsets=NULL;atom_grids=NULL;atom_norbs=NULL;atom_owner_numa=NULL;replica_offsets=NULL;legacy_remote_copy_bytes=0;replica_atoms=0;
  grid_base=NULL;grid_bytes=0;atom_count=0;replica_base=NULL;replica_mapped=replica_copy_bytes=0;
  numa_enabled=0;local_numa=-1;numa_rank=numa_size=numa_domains=0;publish_epoch=0;
  node_leader=-1;world_nodes=0;
  last_copy_rankmax=last_publish_rankmax=0.;
  grid_world=MPI_COMM_NULL;grid_ready=grid_suspended=grid_rebuilding=last_suspend_removed=0;suspend_epoch=0;
}
