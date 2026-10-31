#include "cpu_vna_shared.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
/* Raw projector values retain the native FP32 representation. Heavy payload is
 * stored once on each node needing the atom, never once per world rank. */
#define PANEL_FLOATS ((size_t)262144)
static struct {
 int active,rank,np,atoms,projectors,scratch_fan,scratch_norb,epoch,writing;
 int local_rank,nodes,node_id,raw0;
 int *owner,*fan,*norb,*neighbor,*image,*rank_node;
 unsigned char *need;
 size_t *block,*offset,*elements;
 float *base;
 size_t owned_bytes,cache_bytes,unique_bytes;
 MPI_Comm comm,node,leaders;MPI_Win window;
} V;
static void fail(const char*s){int initialized=0,rank=-1;MPI_Initialized(&initialized);if(initialized)MPI_Comm_rank(MPI_COMM_WORLD,&rank);fprintf(stderr,"[CPU_VNA_SHARED] ERROR rank=%d %s\n",rank,s);fflush(stderr);if(initialized)MPI_Abort(V.active?V.comm:MPI_COMM_WORLD,96);abort();}
static void*alloc(size_t n,size_t item){if(item&&n>SIZE_MAX/item)fail("metadata allocation overflow");void*p=calloc(n?n:1,item);if(!p)fail("allocation failed");return p;}
static size_t plus(size_t a,size_t b){if(b>SIZE_MAX-a)fail("size addition overflow");return a+b;}
static size_t times(size_t a,size_t b){if(a&&b>SIZE_MAX/a)fail("size product overflow");return a*b;}
static unsigned char*needed(int node){return V.need+times((size_t)node,(size_t)V.atoms+1);}
static int atom_node(int a){return V.rank_node[V.owner[a]];}
static void sync_node(void){MPI_Win_sync(V.window);MPI_Barrier(V.node);MPI_Win_sync(V.window);}
int CpuVnaShared_Requested(void){const char*s=getenv("OPENMX_CPU_VNA_SHARED");if(!s||!strcmp(s,"0"))return 0;if(strcmp(s,"1"))fail("OPENMX_CPU_VNA_SHARED must be exactly0 or1");return 1;}
static int raw0_requested(void){const char*s=getenv("OPENMX_CPU_VNA_RAW0");if(!s||!strcmp(s,"0"))return 0;if(strcmp(s,"1"))fail("OPENMX_CPU_VNA_RAW0 must be exactly0 or1");return 1;}
int CpuVnaShared_Agree(MPI_Comm comm){int wanted[2]={CpuVnaShared_Requested(),raw0_requested()},total[2]={0,0},np;MPI_Comm_size(comm,&np);MPI_Allreduce(wanted,total,2,MPI_INT,MPI_SUM,comm);for(int k=0;k<2;k++)if(total[k]!=0&&total[k]!=np)fail(k?"OPENMX_CPU_VNA_RAW0 differs across ranks":"OPENMX_CPU_VNA_SHARED differs across ranks");if(total[1]&&total[0]!=np)fail("OPENMX_CPU_VNA_RAW0 requires OPENMX_CPU_VNA_SHARED=1");return total[0]==np;}
int CpuVnaShared_Active(void){return V.active;}
size_t CpuVnaShared_OwnedBytes(void){return V.active?V.owned_bytes:0;}
void CpuVnaShared_Init(MPI_Comm comm,int atoms,const int*owners,const int*fan,int**natn,int**ncn,int**image4,const int*norb,int projectors,int scratch_fan,int scratch_norb){
 if(V.active)fail("reinitialization before destroy");
 V.comm=comm;MPI_Comm_rank(comm,&V.rank);MPI_Comm_size(comm,&V.np);V.raw0=raw0_requested();
 if(atoms<0||projectors<1||scratch_fan<0||scratch_norb<1||!owners||!fan||!natn||!ncn||!image4||!norb)fail("invalid initialization input");
 MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,V.rank,MPI_INFO_NULL,&V.node);MPI_Comm_rank(V.node,&V.local_rank);
 int leader=V.rank;MPI_Bcast(&leader,1,MPI_INT,0,V.node);int*world_leader=alloc(V.np,sizeof(int));MPI_Allgather(&leader,1,MPI_INT,world_leader,1,MPI_INT,comm);
 V.rank_node=alloc(V.np,sizeof(int));int*roots=alloc(V.np,sizeof(int));for(int p=0;p<V.np;p++)if(world_leader[p]==p)roots[V.nodes++]=p;
 for(int p=0;p<V.np;p++){int j=0;while(j<V.nodes&&roots[j]!=world_leader[p])j++;if(j==V.nodes)fail("invalid shared node leader map");V.rank_node[p]=j;}V.node_id=V.rank_node[V.rank];free(roots);free(world_leader);
 MPI_Comm_split(comm,V.local_rank==0?0:MPI_UNDEFINED,V.rank,&V.leaders);
 V.atoms=atoms;V.projectors=projectors;V.scratch_fan=scratch_fan;V.scratch_norb=scratch_norb;V.epoch=-1;
 V.owner=alloc((size_t)atoms+1,sizeof(int));V.fan=alloc((size_t)atoms+1,sizeof(int));V.norb=alloc((size_t)atoms+1,sizeof(int));
 V.block=alloc((size_t)atoms+2,sizeof(size_t));V.offset=alloc((size_t)atoms+1,sizeof(size_t));V.elements=alloc((size_t)atoms+1,sizeof(size_t));
 for(int a=1;a<=atoms;a++){
  if(owners[a]<0||owners[a]>=V.np||fan[a]<0||fan[a]>scratch_fan||norb[a]<1||norb[a]>scratch_norb||!natn[a]||!ncn[a])fail("invalid atom or scratch capacity");
  V.owner[a]=owners[a];V.fan[a]=fan[a];V.norb[a]=norb[a];V.block[a+1]=plus(V.block[a],(size_t)fan[a]+1);
  V.elements[a]=times(times(times(4,(size_t)fan[a]+1),norb[a]),projectors);V.offset[a]=SIZE_MAX;
 }
 size_t blocks=V.block[atoms+1];V.neighbor=alloc(blocks,sizeof(int));V.image=alloc(times(blocks,3),sizeof(int));
 for(int a=1;a<=atoms;a++)for(int h=0;h<=fan[a];h++){
  size_t at=V.block[a]+h;int b=natn[a][h],cell=ncn[a][h];if(b<1||b>atoms||cell<0||!image4[cell])fail("invalid neighbor or periodic image");V.neighbor[at]=b;
  for(int k=0;k<3;k++)V.image[3*at+k]=image4[cell][k+1];
 }
 /* Replicated metadata must agree before deriving matching network streams. */
 unsigned long long hash=1469598103934665603ULL,low,high;
 #define HASH_VALUE(x) do{hash^=(unsigned long long)(x);hash*=1099511628211ULL;}while(0)
 HASH_VALUE(atoms);HASH_VALUE(projectors);HASH_VALUE(V.raw0);for(int a=1;a<=atoms;a++){HASH_VALUE(V.owner[a]);HASH_VALUE(V.fan[a]);HASH_VALUE(V.norb[a]);for(int h=0;h<=V.fan[a];h++){size_t b=V.block[a]+h;HASH_VALUE(V.neighbor[b]);for(int k=0;k<3;k++)HASH_VALUE(V.image[3*b+k]);}}
 #undef HASH_VALUE
 MPI_Allreduce(&hash,&low,1,MPI_UNSIGNED_LONG_LONG,MPI_MIN,comm);MPI_Allreduce(&hash,&high,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);if(low!=high)fail("inconsistent replicated atom/neighbor metadata");
 V.need=alloc(times((size_t)V.nodes,(size_t)atoms+1),sizeof(unsigned char));
 for(int a=1;a<=atoms;a++){unsigned char*need=needed(atom_node(a));need[a]=1;for(int h=0;h<=V.fan[a];h++)need[V.neighbor[V.block[a]+h]]=1;}
 size_t count=0;for(int a=1;a<=atoms;a++){if(needed(V.node_id)[a]){V.offset[a]=count;count=plus(count,V.elements[a]);}if(atom_node(a)==V.node_id)V.unique_bytes=plus(V.unique_bytes,times(V.elements[a],sizeof(float)));}
 V.cache_bytes=times(count,sizeof(float));V.owned_bytes=V.local_rank==0?V.cache_bytes:0;
 if(V.owned_bytes>(size_t)PTRDIFF_MAX)fail("MPI_Aint allocation overflow");
 float*local=NULL;int rc=MPI_Win_allocate_shared((MPI_Aint)V.owned_bytes,sizeof(float),MPI_INFO_NULL,V.node,&local,&V.window);if(rc!=MPI_SUCCESS)fail("MPI shared allocation failed");
 MPI_Aint size;int disp;void*base=NULL;MPI_Win_shared_query(V.window,0,&size,&disp,&base);if(size<0||(size_t)size!=V.cache_bytes||disp!=sizeof(float)||(size&&!base))fail("shared segment query mismatch");V.base=base;
 MPI_Win_lock_all(0,V.window);V.active=1;
 unsigned long long local_bytes=V.owned_bytes,total=0;MPI_Reduce(&local_bytes,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,comm);
 if(!V.local_rank){char host[MPI_MAX_PROCESSOR_NAME];int n;MPI_Get_processor_name(host,&n);printf("[CPU_VNA_SHARED_NODE] node=%d host=%s cached_payload_bytes=%zu unique_owner_payload_bytes=%zu needed_atoms=",V.node_id,host,V.cache_bytes,V.unique_bytes);int nneed=0;for(int a=1;a<=atoms;a++)nneed+=needed(V.node_id)[a];printf("%d total_atoms=%d\n",nneed,atoms);fflush(stdout);}
 if(V.raw0&&!V.rank){printf("[CPU_VNA_RAW0] enabled=1 raw_components=1 weighted_derivative_components=3 cache_components=4 capacity_unchanged=1\n");fflush(stdout);}
 if(!V.rank){printf("[CPU_VNA_SHARED] active=1 ranks=%d nodes=%d projectors=%d owned_payload_bytes=%llu scratch_payload_bytes=0 FP32_components=4 needed_node_cache=1\n",V.np,V.nodes,projectors,total);fflush(stdout);}
}
float*CpuVnaShared_Row(int k,int a,int h,int i){
 if(!V.active||k<0||k>=4||a<1||a>V.atoms||h<0||h>V.fan[a]||i<0||i>=V.norb[a])return NULL;
 if(V.offset[a]==SIZE_MAX)fail("requested atom outside local node needed set");
 size_t off=V.offset[a]+(((size_t)k*(V.fan[a]+1)+h)*V.norb[a]+i)*V.projectors;
 return V.base+off;
}
/* Iterate matching global-atom streams without a full packed peer payload. */
struct stream {int source,dest,atom,epoch;size_t off;};
static int included(const struct stream*s,int a){return atom_node(a)==s->source&&needed(s->dest)[a];}
/* Raw contraction reads only comp0. Owner derivatives stay resident until
 * their unchanged epoch1 publication. Capacity/layout remain four components. */
static size_t raw_length(int a){return V.raw0?V.elements[a]/4:V.elements[a];}
static size_t stream_total(int source,int dest,int epoch){size_t n=0;for(int a=1;a<=V.atoms;a++)if(atom_node(a)==source&&needed(dest)[a])n=plus(n,epoch?times(V.elements[a]/4,3):raw_length(a));return n;}
static void stream_copy(struct stream*s,float*panel,size_t n,int receive){
 size_t pos=0;while(pos<n){while(s->atom<=V.atoms&&!included(s,s->atom))s->atom++;if(s->atom>V.atoms)fail("network stream exceeded planned atom metadata");int a=s->atom;size_t component=V.elements[a]/4;size_t length=s->epoch?times(component,3):raw_length(a);size_t remain=length-s->off,take=n-pos;if(take>remain)take=remain;float*p=V.base+V.offset[a]+(s->epoch?component:0)+s->off;if(receive)memcpy(p,panel+pos,times(take,sizeof(float)));else memcpy(panel+pos,p,times(take,sizeof(float)));pos+=take;s->off+=take;if(s->off==length){s->atom++;s->off=0;}}
}
static unsigned long long exchange(int epoch){
 if(V.local_rank||V.nodes<2)return 0;
 float*send=alloc(PANEL_FLOATS,sizeof(float)),*recv=alloc(PANEL_FLOATS,sizeof(float));unsigned long long sent=0;
 for(int step=1;step<V.nodes;step++){
  int dest=(V.node_id+step)%V.nodes,src=(V.node_id-step+V.nodes)%V.nodes;
  size_t nt=stream_total(V.node_id,dest,epoch),nr=stream_total(src,V.node_id,epoch),largest=nt>nr?nt:nr;
  unsigned long long rounds=(unsigned long long)(largest/PANEL_FLOATS+(largest%PANEL_FLOATS!=0)),maxrounds=0;MPI_Allreduce(&rounds,&maxrounds,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,V.leaders);
  struct stream sw={V.node_id,dest,1,epoch,0},rw={src,V.node_id,1,epoch,0};
  for(unsigned long long round=0;round<maxrounds;round++){
   int ns=(int)(nt>PANEL_FLOATS?PANEL_FLOATS:nt),nc=(int)(nr>PANEL_FLOATS?PANEL_FLOATS:nr);
   stream_copy(&sw,send,(size_t)ns,0);
   int rc=MPI_Sendrecv(send,ns,MPI_FLOAT,dest,96,recv,nc,MPI_FLOAT,src,96,V.leaders,MPI_STATUS_IGNORE);if(rc!=MPI_SUCCESS)fail("needed-node projector transfer failed");
   stream_copy(&rw,recv,(size_t)nc,1);nt-=ns;nr-=nc;sent+=(unsigned long long)ns*sizeof(float);
  }
 }
 free(send);free(recv);return sent;
}
void CpuVnaShared_Publish(int epoch){
 if(!V.active||epoch<0||epoch>1||epoch!=V.epoch+1||(epoch==1&&!V.writing))fail("invalid raw/weighted publication order");
 double start=MPI_Wtime();sync_node();MPI_Barrier(V.comm);
 unsigned long long sent=exchange(epoch),total=0;sync_node();MPI_Barrier(V.comm);V.epoch=epoch;V.writing=0;
 double local=MPI_Wtime()-start,critical=0;MPI_Reduce(&local,&critical,1,MPI_DOUBLE,MPI_MAX,0,V.comm);MPI_Reduce(&sent,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,V.comm);
 if(!V.rank){printf("[CPU_VNA_SHARED] publish=%s epoch=%d network_payload_bytes=%llu critical_seconds=%.9f\n",epoch?"weighted_derivatives":"raw",epoch,total,critical);fflush(stdout);}
}
void CpuVnaShared_BeginWeighted(void){if(!V.active||V.epoch!=0||V.writing)fail("invalid weighted writer epoch");/* Cross-node raw readers must finish before any owner mutates its source. */MPI_Barrier(V.comm);sync_node();V.writing=1;}
void CpuVnaShared_BindAtom(float*****ds,int scratch,int a){
 if(!V.active||V.epoch<0||V.writing||!ds||scratch<0||a<1||a>V.atoms)fail("invalid complete-atom view");
 for(int k=0;k<4;k++)for(int h=0;h<=V.fan[a];h++)for(int i=0;i<V.norb[a];i++)ds[k][scratch][h][i]=(V.raw0&&V.epoch==0&&k)?NULL:CpuVnaShared_Row(k,a,h,i);
}
void CpuVnaShared_BindCenter(float*****ds,int scratch,int a,int**reverse){
 if(!V.active||V.epoch!=1||V.writing||!ds||scratch<0||a<1||a>V.atoms||!reverse)fail("invalid weighted center view");
 for(int h=0;h<=V.fan[a];h++){
  size_t at=V.block[a]+h;int b=V.neighbor[at];if(!reverse[h])fail("null reverse neighbor row");int rh=reverse[h][0];if(rh<0||rh>V.fan[b])fail("invalid reverse neighbor slot");size_t rev=V.block[b]+rh;
  if(V.neighbor[rev]!=a)fail("reverse neighbor has wrong atom");
  for(int xyz=0;xyz<3;xyz++)if((int64_t)V.image[3*at+xyz]+V.image[3*rev+xyz]!=0)fail("reverse neighbor has wrong periodic image");
  for(int k=0;k<4;k++)for(int i=0;i<V.norb[b];i++)ds[k][scratch][h][i]=CpuVnaShared_Row(k,b,rh,i);
 }
}
void CpuVnaShared_Destroy(void){
 if(!V.active)return;
 MPI_Barrier(V.comm);sync_node();MPI_Win_unlock_all(V.window);MPI_Win_free(&V.window);if(V.leaders!=MPI_COMM_NULL)MPI_Comm_free(&V.leaders);MPI_Comm_free(&V.node);
 free(V.owner);free(V.fan);free(V.norb);free(V.block);free(V.offset);free(V.elements);free(V.rank_node);free(V.need);free(V.neighbor);free(V.image);memset(&V,0,sizeof(V));
}
