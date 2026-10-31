#ifndef CPU_NSOC_SECTOR_LAYOUT_H
#define CPU_NSOC_SECTOR_LAYOUT_H
/* Fixed-size CPU/MPI proof metadata only; no dense electronic matrices. */
#include <mpi.h>
#include <omp.h>
#include <mkl.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
extern long syscall(long,...);
#define NSLP_MAX_THREADS 8
#define NSLP_MASK_WORDS 16
#define NSLP_MAX_LOCAL_RANKS 64
typedef struct {
 int requested,actual,mkl_threads,valid,cpu[NSLP_MAX_THREADS];
} CpuNSectorThreadProof;
static int sector_scope_valid(int world_size,int nodes,int local_ranks,int threads,int mkl_threads){
 return nodes==4&&threads==mkl_threads&&
  ((world_size==256&&local_ranks==64&&threads==1)||
   (world_size==64&&local_ranks==16&&threads==4));
}
static int sector_thread_probe(int threads,CpuNSectorThreadProof*proof){
 memset(proof,0,sizeof(*proof));proof->requested=threads;proof->mkl_threads=mkl_get_max_threads();
 for(int i=0;i<NSLP_MAX_THREADS;i++)proof->cpu[i]=-1;
 if(threads<1||threads>NSLP_MAX_THREADS||omp_get_max_threads()!=threads||proof->mkl_threads!=threads||omp_in_parallel())return 0;
 unsigned long masks[NSLP_MAX_THREADS][NSLP_MASK_WORDS];int okay[NSLP_MAX_THREADS];memset(masks,0,sizeof(masks));memset(okay,0,sizeof(okay));
 #pragma omp parallel num_threads(threads) shared(masks,okay,proof)
 {
  int tid=omp_get_thread_num();
  #pragma omp single
  proof->actual=omp_get_num_threads();
  if(tid<NSLP_MAX_THREADS)okay[tid]=syscall(SYS_sched_getaffinity,0,sizeof(masks[tid]),masks[tid])>0;
 }
 if(proof->actual!=threads)return 0;
 for(int t=0;t<threads;t++){
  if(!okay[t])return 0;int count=0;
  for(int word=0;word<NSLP_MASK_WORDS;word++)for(unsigned bit=0;bit<8*sizeof(unsigned long);bit++)if(masks[t][word]&(1UL<<bit)){count++;proof->cpu[t]=word*8*sizeof(unsigned long)+bit;}
  if(count!=1)return 0;
  for(int q=0;q<t;q++)if(proof->cpu[q]==proof->cpu[t])return 0;
 }
 proof->valid=1;return 1;
}
static int sector_node_layout(MPI_Comm node,int expected_ranks,int threads,int expected_cores,const CpuNSectorThreadProof*proof){
 int peers;MPI_Comm_size(node,&peers);
 if(peers<1||peers>NSLP_MAX_LOCAL_RANKS||expected_ranks!=peers||threads<1||threads>NSLP_MAX_THREADS||expected_cores!=peers*threads||expected_cores>64)return 0;
 int packet[2+NSLP_MAX_THREADS],all[NSLP_MAX_LOCAL_RANKS][2+NSLP_MAX_THREADS];
 packet[0]=proof->valid;packet[1]=proof->actual;for(int t=0;t<NSLP_MAX_THREADS;t++)packet[2+t]=proof->cpu[t];
 MPI_Allgather(packet,2+NSLP_MAX_THREADS,MPI_INT,all,2+NSLP_MAX_THREADS,MPI_INT,node);
 int used[64]={0};
 for(int r=0;r<peers;r++){
  if(all[r][0]!=1||all[r][1]!=threads)return 0;
  for(int t=0;t<threads;t++){int cpu=all[r][2+t];if(cpu<0||cpu>=expected_cores||used[cpu])return 0;used[cpu]=1;}
 }
 for(int cpu=0;cpu<expected_cores;cpu++)if(used[cpu]!=1)return 0;
 return 1;
}
#endif
