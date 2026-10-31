#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "cpu_grid_fields_shared.unified.h"

int cpu_grid_fields_active=0;
int cpu_grid_fields_width=0;
static int requested=0,rank_global=0,world_size=0,node_rank=0,node_size=0;
static int node_id=0,node_count=0,*rank_nodes=NULL,*counts=NULL,*displacements=NULL;
static MPI_Comm world_comm=MPI_COMM_NULL,node_comm=MPI_COMM_NULL,leaders_comm=MPI_COMM_NULL;
static MPI_Win field_win=MPI_WIN_NULL;
static double *fields=NULL,*send_panel=NULL,*receive_panel=NULL;
static size_t points=0,columns=0,third=0,owned_first=0,owned_count=0,field_bytes=0;
static unsigned long generation[CPU_GF_MAX_FIELDS];
enum { PANEL_POINTS=131072 };
static int components_for(int field){return field==CPU_GF_DENSITY||field==CPU_GF_VXC||field==CPU_GF_VPOT?CPU_GF_VXC:field==CPU_GF_REF||field==CPU_GF_HARTREE?1:0;}
static const char *name_for(int field){return field==CPU_GF_DENSITY?"Density":field==CPU_GF_VXC?"Vxc":field==CPU_GF_VPOT?"Vpot":field==CPU_GF_REF?"RefVxc":"dVHart";}
static int collective_bad(int bad,MPI_Comm comm){int any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,comm);return any;}
/* Avoid multiplying a potentially large global grid by the MPI rank. */
static size_t split_columns(int rank){return columns/(size_t)world_size*(size_t)rank+((columns%(size_t)world_size)*(size_t)rank+(size_t)world_size-1)/(size_t)world_size;}
static void rank_range(int rank,size_t *first,size_t *last){*first=split_columns(rank)*third;*last=split_columns(rank+1)*third;}
static void sync_fields(void){MPI_Win_sync(field_win);MPI_Barrier(node_comm);MPI_Win_sync(field_win);}
static void free_metadata(void)
{
 free(rank_nodes);free(counts);free(displacements);free(send_panel);free(receive_panel);
 rank_nodes=counts=displacements=NULL;send_panel=receive_panel=NULL;
 if(leaders_comm!=MPI_COMM_NULL)MPI_Comm_free(&leaders_comm);
 if(node_comm!=MPI_COMM_NULL)MPI_Comm_free(&node_comm);
 if(world_comm!=MPI_COMM_NULL)MPI_Comm_free(&world_comm);
}
int CpuGridFields_Configure(MPI_Comm world,int components)
{
 int local[2]={components,cpu_grid_fields_active},lo[2],hi[2];
 MPI_Allreduce(local,lo,2,MPI_INT,MPI_MIN,world);
 MPI_Allreduce(local,hi,2,MPI_INT,MPI_MAX,world);
 if(lo[0]!=hi[0]||(lo[0]!=2&&lo[0]!=4)||hi[1])return -1;
 cpu_grid_fields_width=components;
 return 0;
}
int CpuGridFields_Agree(MPI_Comm world,int valid)
{
 const char *s=getenv("OPENMX_CPU_GRID_FIELDS_SHARED");
 int mode=!s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:-1;
 int local[2]={mode,!!valid},lo[2],hi[2];
 MPI_Allreduce(local,lo,2,MPI_INT,MPI_MIN,world);MPI_Allreduce(local,hi,2,MPI_INT,MPI_MAX,world);
 if(lo[0]<0||lo[0]!=hi[0]||(hi[0]&&lo[1]!=1)){requested=0;return -1;}
 requested=mode;return mode;
}
int CpuGridFields_Init(MPI_Comm world,int n1,int n2,int n3,size_t local_B,size_t local_C)
{
 int disp,local[3]={n1,n2,n3},lo[3],hi[3],bad=0;
 MPI_Aint query_bytes;void *base=NULL;unsigned long long old_local,old_total;
 MPI_Comm_rank(world,&rank_global);MPI_Comm_size(world,&world_size);
 MPI_Allreduce(local,lo,3,MPI_INT,MPI_MIN,world);MPI_Allreduce(local,hi,3,MPI_INT,MPI_MAX,world);
 if(!requested||cpu_grid_fields_active||n1<=0||n2<=0||n3<=0||
    (cpu_grid_fields_width!=2&&cpu_grid_fields_width!=4))bad=1;
 for(int i=0;i<3;i++)if(lo[i]!=hi[i])bad=1;
 if(!bad&&((size_t)n1>SIZE_MAX/(size_t)n2||(size_t)n1*n2>SIZE_MAX/(size_t)n3))bad=1;
 if(collective_bad(bad,world))return -1;
 columns=(size_t)n1*n2;third=(size_t)n3;points=columns*third;
 if(points>SIZE_MAX/CPU_GF_FIELDS/sizeof(double)||points>PTRDIFF_MAX/CPU_GF_FIELDS/sizeof(double))bad=1;
 size_t owned_last;rank_range(rank_global,&owned_first,&owned_last);owned_count=owned_last-owned_first;
 if(local_B!=owned_count||local_C>ULLONG_MAX/(unsigned)world_size/CPU_GF_FIELDS/sizeof(double))bad=1;
 if(collective_bad(bad,world))return -1;
 field_bytes=CPU_GF_FIELDS*points*sizeof(double);
 if((MPI_Aint)field_bytes<0||(size_t)(MPI_Aint)field_bytes!=field_bytes)bad=1;
 if(collective_bad(bad,world))return -1;
 MPI_Comm_dup(world,&world_comm);
 MPI_Comm_split_type(world_comm,MPI_COMM_TYPE_SHARED,rank_global,MPI_INFO_NULL,&node_comm);
 MPI_Comm_rank(node_comm,&node_rank);MPI_Comm_size(node_comm,&node_size);
 MPI_Comm_split(world_comm,node_rank==0?0:MPI_UNDEFINED,rank_global,&leaders_comm);
 if(node_rank==0){MPI_Comm_rank(leaders_comm,&node_id);MPI_Comm_size(leaders_comm,&node_count);}
 MPI_Bcast(&node_id,1,MPI_INT,0,node_comm);MPI_Bcast(&node_count,1,MPI_INT,0,node_comm);
 if(node_count>1){
  rank_nodes=(int*)malloc((size_t)world_size*sizeof(int));
  if(node_rank==0){counts=(int*)malloc((size_t)node_count*sizeof(int));displacements=(int*)malloc((size_t)node_count*sizeof(int));send_panel=(double*)malloc((size_t)PANEL_POINTS*sizeof(double));receive_panel=(double*)malloc((size_t)PANEL_POINTS*sizeof(double));}
  bad=!rank_nodes||(node_rank==0&&(!counts||!displacements||!send_panel||!receive_panel));
  if(collective_bad(bad,world_comm)){free_metadata();return -2;}
  MPI_Allgather(&node_id,1,MPI_INT,rank_nodes,1,MPI_INT,world_comm);
 }
 MPI_Win_allocate_shared(node_rank==0?(MPI_Aint)field_bytes:0,sizeof(double),MPI_INFO_NULL,node_comm,&base,&field_win);
 MPI_Win_shared_query(field_win,0,&query_bytes,&disp,&base);
 bad=query_bytes!=(MPI_Aint)field_bytes||disp!=(int)sizeof(double)||!base;
 if(collective_bad(bad,world_comm)){MPI_Win_free(&field_win);field_win=MPI_WIN_NULL;free_metadata();return -3;}
 fields=(double*)base;MPI_Win_lock_all(0,field_win);
 /* On every physical node initialize each GN once, including remote-B ranges. */
 size_t first=points/(size_t)node_size*(size_t)node_rank+(points%(size_t)node_size)*(size_t)node_rank/(size_t)node_size;
 size_t last=points/(size_t)node_size*(size_t)(node_rank+1)+(points%(size_t)node_size)*(size_t)(node_rank+1)/(size_t)node_size;
 for(int c=0;c<components_for(CPU_GF_DENSITY);c++)memset(fields+(size_t)c*points+first,0,(last-first)*sizeof(double));
 memset(fields+(size_t)CPU_GF_HARTREE*points+first,0,(last-first)*sizeof(double));
 memset(generation,0,sizeof(generation));sync_fields();cpu_grid_fields_active=1;
 old_local=(unsigned long long)local_C*CPU_GF_FIELDS*sizeof(double);
 MPI_Reduce(&old_local,&old_total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,world_comm);
 if(rank_global==0){printf("[CPU_GRID_FIELDS_SHARED] active=1 fields=%d points=%zu shared_bytes=%zu old_C_bytes=%llu single_node=%d nodes=%d leader_scratch_bytes=%zu\n",CPU_GF_FIELDS,points,field_bytes,old_total,node_count==1,node_count,node_count==1?(size_t)0:2*(size_t)PANEL_POINTS*sizeof(double));fflush(stdout);}
 return 0;
}
int CpuGridFields_Active(void){return cpu_grid_fields_active;}
int CpuGridFields_Requested(void){return requested;}
size_t CpuGridFields_Points(void){return points;}
size_t CpuGridFields_Bytes(void){return field_bytes;}
double *CpuGridFields_Row(int field,int component)
{
 if(!cpu_grid_fields_active||component<0||component>=components_for(field))return NULL;
 return fields+(size_t)(field+component)*points;
}
int CpuGridFields_Owns(const double *ptr)
{
 if(!cpu_grid_fields_active||!ptr)return 0;
 uintptr_t address=(uintptr_t)ptr,start=(uintptr_t)fields;
 return address>=start&&address-start<field_bytes;
}
static int find_field(const double *ptr)
{
 if(!cpu_grid_fields_active)return -1;
 for(int field=0;field<CPU_GF_FIELDS;field++)if(ptr==fields+(size_t)field*points)return field;
 return -1;
}
static void exchange_panels(int field,int components)
{
 if(node_count==1||node_rank!=0)return;
 for(int c=0;c<components;c++)for(size_t panel=0;panel<points;){
  size_t width=points-panel<PANEL_POINTS?points-panel:PANEL_POINTS,last=panel+width,sent=0;
  memset(counts,0,(size_t)node_count*sizeof(int));
  for(int r=0;r<world_size;r++){
   size_t first,end;rank_range(r,&first,&end);if(first<panel)first=panel;if(end>last)end=last;
   if(end>first){size_t n=end-first;counts[rank_nodes[r]]+=(int)n;if(rank_nodes[r]==node_id){memcpy(send_panel+sent,fields+(size_t)(field+c)*points+first,n*sizeof(double));sent+=n;}}
  }
  int total=0;for(int n=0;n<node_count;n++){displacements[n]=total;total+=counts[n];}
  /* All counts/displacements are bounded by width <= PANEL_POINTS < INT_MAX. */
  if(total!=(int)width||sent!=(size_t)counts[node_id])MPI_Abort(world_comm,97);
  MPI_Allgatherv(send_panel,counts[node_id],MPI_DOUBLE,receive_panel,counts,displacements,MPI_DOUBLE,leaders_comm);
  memcpy(counts,displacements,(size_t)node_count*sizeof(int));
  for(int r=0;r<world_size;r++){
   size_t first,end;rank_range(r,&first,&end);if(first<panel)first=panel;if(end>last)end=last;
   if(end>first){size_t n=end-first;int node=rank_nodes[r];memcpy(fields+(size_t)(field+c)*points+first,receive_panel+counts[node],n*sizeof(double));counts[node]+=(int)n;}
  }
  panel=last;
 }
}
static int publish(double **source,double **target,int components,int vector_call)
{
 if(!cpu_grid_fields_active)return 0;
 int field=target?find_field(target[0]):-1,bad=0;
 unsigned long meta[3]={(unsigned long)(field+1),(unsigned long)components,field>=0?generation[field]:0},lo[3],hi[3];
 MPI_Allreduce(meta,lo,3,MPI_UNSIGNED_LONG,MPI_MIN,world_comm);MPI_Allreduce(meta,hi,3,MPI_UNSIGNED_LONG,MPI_MAX,world_comm);
 for(int k=0;k<3;k++)if(lo[k]!=hi[k])return -1;
 if(field<0)return 0;
 int expected=components_for(field),duplicate=vector_call&&CPU_GF_FIELDS==8&&expected==2&&components==1;
 if(!source||(!duplicate&&expected!=components)||!expected)bad=1;
 if(!bad)for(int c=0;c<expected;c++)if((owned_count&&!source[c])||target[c]!=fields+(size_t)(field+c)*points)bad=1;
 if(collective_bad(bad,world_comm))return -1;
 if(hi[2]==ULONG_MAX)return -1;
 MPI_Barrier(node_comm); /* prior readers have left this epoch */
 if(owned_count){
  if(duplicate)memcpy(source[1],source[0],owned_count*sizeof(double));
  for(int c=0;c<expected;c++)memcpy(fields+(size_t)(field+c)*points+owned_first,source[c],owned_count*sizeof(double));
 }
 sync_fields(); /* publish B owners to the node leader before its packing */
 exchange_panels(field,expected);
 if(node_count>1)sync_fields(); /* leader writes become visible to all consumers */
 generation[field]++;
 if(rank_global==0){printf("[CPU_GRID_FIELDS_SHARED_WRITE] field=%s components=%d generation=%lu points=%zu source=B snapshot=1 nodes=%d\n",name_for(field),expected,generation[field],points,node_count);fflush(stdout);}
 return 1;
}
int CpuGridFields_Copy2(double **source,double **target,int components){return publish(source,target,components,1);}
int CpuGridFields_Copy1(double *source,double *target){double *s[1]={source},*t[1]={target};return publish(s,t,1,0);}
void CpuGridFields_Destroy(void)
{
 if(!cpu_grid_fields_active)return;
 sync_fields();MPI_Barrier(node_comm);MPI_Win_unlock_all(field_win);MPI_Win_free(&field_win);
 field_win=MPI_WIN_NULL;free_metadata();fields=NULL;points=columns=third=owned_first=owned_count=field_bytes=0;
 node_count=node_id=node_rank=node_size=world_size=0;cpu_grid_fields_active=0;
}
