#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <malloc.h>
#include <sys/mman.h>
#include <time.h>
#include "openmx_common.h"
#include "cpu_memory_phase.h"
#include "cpu_memory_inventory.h"
extern size_t CpuC3_CompactBytes(void) __attribute__((weak));
extern void CpuAllocationSnapshot(const char*) __attribute__((weak));
extern dcomplex *EVec1_NonCol,*Ss2_Cx,*Hs2_Cx,*Cs2_Cx;
extern double *CDM1,*Work1;
extern double *Ss_Re,*Cs_Re,*rHs11_Re,*rHs12_Re,*rHs22_Re,*iHs11_Re,*iHs12_Re,*iHs22_Re;
static unsigned long long resident(void *p,size_t bytes)
{
  long page=sysconf(_SC_PAGESIZE);uintptr_t lo;size_t span,pages,i;unsigned char *vec;unsigned long long total=0;
  if(!p||!bytes)return 0;
  lo=(uintptr_t)p & ~((uintptr_t)page-1);
  if(bytes>SIZE_MAX-((uintptr_t)p-lo)-(size_t)page)return 0;
  span=bytes+((uintptr_t)p-lo);pages=(span+(size_t)page-1)/(size_t)page;
  vec=mmap(NULL,pages,PROT_READ|PROT_WRITE,MAP_PRIVATE|MAP_ANONYMOUS,-1,0);
  if(vec==MAP_FAILED)return 0;
  if(mincore((void*)lo,pages*(size_t)page,vec)==0)for(i=0;i<pages;i++)if(vec[i]&1)total+=(unsigned long long)page;
  munmap(vec,pages);return total;
}
static void pss(unsigned long long x[6])
{
  FILE *f=fopen("/proc/self/smaps","r");char line[8192],name[4096]="",mode[8];unsigned long a,b,v;int cat=3;
  memset(x,0,6*sizeof(*x));if(!f)return;
  while(fgets(line,sizeof(line),f)){
    int n=0;
    if(sscanf(line,"%lx-%lx %7s %*s %*s %*s %n",&a,&b,mode,&n)==3 && n>0){
      char *q=line+n;while(*q==' '||*q=='\t')q++;
      strncpy(name,q,sizeof(name)-1);name[sizeof(name)-1]=0;
      cat=strstr(name,"[heap]")?0:mode[3]=='s'?2:(*q=='\n'||*q=='\0'||*q=='[')?1:3;
    }else if(sscanf(line,"Pss: %lu kB",&v)==1){x[0]+=(unsigned long long)v*1024;x[cat+1]+=(unsigned long long)v*1024;}
    else if(sscanf(line,"Rss: %lu kB",&v)==1)x[5]+=(unsigned long long)v*1024;
  }
  fclose(f);
}
void CpuMemoryPhase(const char *phase)
{
  const char *on=getenv("OPENMX_MEMORY_PROBE");int rank;size_t i;double started;struct timespec now;
  void *ptrs[14];const char *labels[14]={"EVec1_NonCol","Ss2_Cx","Hs2_Cx","Cs2_Cx","Ss_Re","Cs_Re","rHs11","rHs12","rHs22","iHs11","iHs12","iHs22","CDM1","Work1"};
  unsigned long long local[47]={0},sum[47],max[47];FILE *f;char path[4096];
  if(!on||strcmp(on,"1"))return;
  started=MPI_Wtime();MPI_Barrier(mpi_comm_level1);MPI_Comm_rank(mpi_comm_level1,&rank);
  ptrs[0]=EVec1_NonCol;ptrs[1]=Ss2_Cx;ptrs[2]=Hs2_Cx;ptrs[3]=Cs2_Cx;ptrs[4]=Ss_Re;ptrs[5]=Cs_Re;
  ptrs[6]=rHs11_Re;ptrs[7]=rHs12_Re;ptrs[8]=rHs22_Re;ptrs[9]=iHs11_Re;ptrs[10]=iHs12_Re;ptrs[11]=iHs22_Re;ptrs[12]=CDM1;ptrs[13]=Work1;
  pss(local);
  for(i=0;i<14;i++){size_t bytes=ptrs[i]?malloc_usable_size(ptrs[i]):0;local[6+2*i]=bytes;local[7+2*i]=resident(ptrs[i],bytes);}
  int inventory=!strcmp(phase,"after_cluster_allocation")||!strcmp(phase,"after_force");
  if(inventory){
    CpuMemoryHistoryInventory(local+34,local+40);
    snprintf(path,sizeof(path),"memory_allocator_%s_rank%04d.xml",phase,rank);
    f=fopen(path,"w");if(f){malloc_info(0,f);fclose(f);}
  }
  if(CpuC3_CompactBytes)local[46]=CpuC3_CompactBytes();
  clock_gettime(CLOCK_REALTIME,&now);
  snprintf(path,sizeof(path),"memory_phase_rank%04d.jsonl",rank);f=fopen(path,"a");
  if(f){fprintf(f,"{\"stage\":\"%s\",\"unix_seconds\":%ld.%09ld,\"rank\":%d,\"PSS_bytes\":%llu,\"heap_PSS_bytes\":%llu,\"anon_PSS_bytes\":%llu,\"shared_PSS_bytes\":%llu,\"file_PSS_bytes\":%llu,\"RSS_bytes\":%llu,\"arrays\":{",phase,(long)now.tv_sec,now.tv_nsec,rank,local[0],local[1],local[2],local[3],local[4],local[5]);
    for(i=0;i<14;i++)fprintf(f,"%s\"%s\":{\"allocated_bytes\":%llu,\"resident_page_span_upper_bound_bytes\":%llu}",i?",":"",labels[i],local[6+2*i],local[7+2*i]);
    fprintf(f,"},\"history_inventory_measured\":%d,\"history_usable_bytes\":[",inventory);
    for(i=0;i<6;i++)fprintf(f,"%s%llu",i?",":"",local[34+i]);
    fprintf(f,"],\"history_allocation_counts\":[");for(i=0;i<6;i++)fprintf(f,"%s%llu",i?",":"",local[40+i]);
    fprintf(f,"],\"compact_vector_payload_bytes\":%llu}\n",local[46]);fclose(f);}
  MPI_Reduce(local,sum,47,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,mpi_comm_level1);
  MPI_Reduce(local,max,47,MPI_UNSIGNED_LONG_LONG,MPI_MAX,0,mpi_comm_level1);
  if(rank==0){
    f=fopen("memory_phase_global.jsonl","a");if(f){
      fprintf(f,"{\"stage\":\"%s\",\"unix_seconds\":%ld.%09ld,\"PSS_bytes\":%llu,\"heap_PSS_bytes\":%llu,\"anon_PSS_bytes\":%llu,\"shared_PSS_bytes\":%llu,\"file_PSS_bytes\":%llu,\"RSS_sum_bytes\":%llu,\"max_rank_PSS_bytes\":%llu,\"arrays\":{",phase,(long)now.tv_sec,now.tv_nsec,sum[0],sum[1],sum[2],sum[3],sum[4],sum[5],max[0]);
      for(i=0;i<14;i++)fprintf(f,"%s\"%s\":{\"allocated_bytes\":%llu,\"resident_page_span_upper_bound_bytes\":%llu}",i?",":"",labels[i],sum[6+2*i],sum[7+2*i]);
      fprintf(f,"},\"history_inventory_measured\":%d,\"history_usable_bytes\":[",inventory);
      for(i=0;i<6;i++)fprintf(f,"%s%llu",i?",":"",sum[34+i]);
      fprintf(f,"],\"history_allocation_counts\":[");for(i=0;i<6;i++)fprintf(f,"%s%llu",i?",":"",sum[40+i]);
      fprintf(f,"],\"compact_vector_payload_bytes\":%llu,\"probe_wall_rank0_seconds\":%.9f}\n",sum[46],MPI_Wtime()-started);fclose(f);}
    printf("CPU_MEMORY_PHASE stage=%s PSS_GiB=%.6f heap_PSS_GiB=%.6f resident_EVec_GiB=%.6f resident_Q_GiB=%.6f probe_rank0_seconds=%.6f\n",phase,sum[0]/1073741824.,sum[1]/1073741824.,sum[7]/1073741824.,sum[13]/1073741824.,MPI_Wtime()-started);fflush(stdout);
  }
  if(CpuAllocationSnapshot && (inventory||!strcmp(phase,"before_force")))CpuAllocationSnapshot(phase);
}
