#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sched.h>
#include <unistd.h>
#include <sys/syscall.h>
#include "openmx_common.h"
#include "cpu_grid_shared.h"
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,np,lr,leader;MPI_Comm local;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&local);MPI_Comm_rank(local,&lr);leader=rank;MPI_Allreduce(MPI_IN_PLACE,&leader,1,MPI_INT,MPI_MIN,local);
 if(argc<2)MPI_Abort(MPI_COMM_WORLD,2);char*list=strdup(argv[1]),*word=strtok(list,",");for(int i=0;i<lr&&word;i++)word=strtok(NULL,",");if(!word)MPI_Abort(MPI_COMM_WORLD,3);int cpu=atoi(word);cpu_set_t mask;CPU_ZERO(&mask);CPU_SET(cpu,&mask);if(sched_setaffinity(0,sizeof(mask),&mask))MPI_Abort(MPI_COMM_WORLD,4);free(list);unsigned actual,node;syscall(SYS_getcpu,&actual,&node,NULL);
 int id[2]={leader,(int)node},*domains=malloc((size_t)np*2*sizeof(int));MPI_Allgather(id,2,MPI_INT,domains,2,MPI_INT,MPI_COMM_WORLD);int numa_count=0;for(int r=0;r<np;r++){int seen=0;for(int q=0;q<r;q++)if(domains[2*q]==domains[2*r]&&domains[2*q+1]==domains[2*r+1])seen=1;if(!seen&&domains[2*r]==leader)numa_count++;}if(numa_count<2){if(!rank)fprintf(stderr,"test needs two physical NUMA domains per node\n");MPI_Abort(MPI_COMM_WORLD,5);}
 int n=4*np;int*owners=calloc(n+1,4),*grids=calloc(n+1,4),*species=calloc(n+1,4),norbs[2]={2,3};FNAN=calloc(n+1,4);natn=calloc(n+1,sizeof(int*));
 for(int a=1;a<=n;a++){owners[a]=(a-1)/4;grids[a]=3+a*2;species[a]=a%2;natn[a]=calloc(n+1,4);}
 setenv("OPENMX_CPU_GRID_SHARED","1",1);setenv("OPENMX_CPU_GRID_NUMA","1",1);int errors=0;size_t total_expected=0,total_old=0;
 for(int kind=0;kind<3;kind++){
  for(int a=1;a<=n;a++){natn[a][0]=a;FNAN[a]=kind==0?1:kind==1?0:n-1;if(kind==0)natn[a][1]=a%n+1;else if(kind==2){int h=1;for(int b=1;b<=n;b++)if(b!=a)natn[a][h++]=b;}}
  unsigned char*needed=calloc(n+1,1);for(int a=1;a<=n;a++)if(domains[2*owners[a]]==leader&&domains[2*owners[a]+1]==(int)node)for(int h=0;h<=FNAN[a];h++)needed[natn[a][h]]=1;
  size_t expected=0,old=0;for(int a=1;a<=n;a++)if(domains[2*owners[a]]==leader&&domains[2*owners[a]+1]!=(int)node){size_t z=(size_t)grids[a]*norbs[species[a]]*sizeof(Type_Orbs_Grid);old+=z;if(needed[a])expected+=z;}
  CpuGridShared_Init(MPI_COMM_WORLD,n,owners,grids,species,norbs);CpuGridSharedNumaStats stats;CpuGridShared_NumaStats(&stats);
  if(stats.copy_bytes!=expected||stats.mapped_bytes!=expected){if(!rank)fprintf(stderr,"FAIL demand kind=%d expected=%zu actual_copy=%zu actual_mapped=%zu\n",kind,expected,stats.copy_bytes,stats.mapped_bytes);errors++;}
  for(int epoch=1;epoch<=2;epoch++){
   for(int a=1;a<=n;a++)if(owners[a]==rank)for(int g=0;g<grids[a];g++){Type_Orbs_Grid*row=CpuGridShared_Row(a,g);if(!row)MPI_Abort(MPI_COMM_WORLD,6);for(int o=0;o<norbs[species[a]];o++)row[o]=10000*epoch+100*a+10*g+o;}
   CpuGridShared_Publish();
   for(int a=1;a<=n;a++)for(int g=0;g<grids[a];g++){
    Type_Orbs_Grid*canonical=CpuGridShared_Row(a,g),*read=CpuGridShared_ReadRow(a,g);
    if(domains[2*owners[a]]!=leader){if(canonical||read)errors++;continue;}
    int replicated=needed[a]&&domains[2*owners[a]+1]!=(int)node;
    if(!canonical||!read){errors++;continue;}if((read!=canonical)!=replicated)errors++;
    for(int o=0;o<norbs[species[a]];o++)if(read[o]!=10000*epoch+100*a+10*g+o)errors++;
   }
   /* Every reader must finish before owners overwrite the canonical rows. */
   MPI_Barrier(MPI_COMM_WORLD);
  }
  CpuGridShared_Destroy();free(needed);total_expected+=expected;total_old+=old;
 }
 int all;MPI_Allreduce(&errors,&all,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);if(!rank)printf("%s NUMA demand values/pointers/reinit rank0_expected_bytes=%zu legacy_remote_bytes=%zu ranks=%d domains_per_node=%d errors=%d\n",all?"FAIL":"PASS",total_expected,total_old,np,numa_count,all);
 for(int a=1;a<=n;a++)free(natn[a]);free(natn);free(FNAN);free(owners);free(grids);free(species);free(domains);MPI_Comm_free(&local);MPI_Finalize();return all?1:0;
}
