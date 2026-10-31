#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "cpu_grid_fields_shared.unified.h"
static size_t split(size_t n,int r,int p){return n/p*r+((n%p)*(size_t)r+p-1)/p;}
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,np,errors=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
 setenv("OPENMX_CPU_GRID_FIELDS_SHARED","1",1);
 for(int width=2;width<=4;width+=2){
  if(CpuGridFields_Configure(MPI_COMM_WORLD,width)!=0)errors++;
  if(CpuGridFields_Agree(MPI_COMM_WORLD,1)!=1)errors++;
  size_t cols=7,points=21,first=split(cols,rank,np)*3,count=(split(cols,rank+1,np)*3)-first;
  if(CpuGridFields_Init(MPI_COMM_WORLD,1,7,3,count,count))MPI_Abort(MPI_COMM_WORLD,2);
  if(CpuGridFields_Bytes()!=(size_t)(3*width+2)*points*sizeof(double))errors++;
  double*src[4],*dst[4];for(int c=0;c<4;c++)src[c]=calloc(count?count:1,sizeof(double));
  int fields[5]={CPU_GF_DENSITY,CPU_GF_VXC,CPU_GF_VPOT,CPU_GF_REF,CPU_GF_HARTREE};
  for(int epoch=0;epoch<2;epoch++)for(int f=0;f<5;f++){
   int nc=f<3?width:1;
   for(int c=0;c<nc;c++){dst[c]=CpuGridFields_Row(fields[f],c);if(!dst[c])MPI_Abort(MPI_COMM_WORLD,3);for(size_t i=0;i<count;i++)src[c][i]=10000*epoch+1000*f+100*c+first+i;}
   int rc=nc==1?CpuGridFields_Copy1(src[0],dst[0]):CpuGridFields_Copy2(src,dst,nc);
   if(rc!=1)errors++;
   for(int c=0;c<nc;c++)for(size_t i=rank;i<points;i+=np)if(dst[c][i]!=10000*epoch+1000*f+100*c+i)errors++;
  }
  if(width==2){
   for(int c=0;c<2;c++)dst[c]=CpuGridFields_Row(CPU_GF_DENSITY,c);
   for(size_t i=0;i<count;i++){src[0][i]=first+i+0.5;src[1][i]=-1;}
   if(CpuGridFields_Copy2(src,dst,1)!=1)errors++;
   for(size_t i=rank;i<points;i+=np)if(dst[0][i]!=i+0.5||dst[1][i]!=i+0.5)errors++;
  }
  if(CpuGridFields_Configure(MPI_COMM_WORLD,width)!=-1)errors++;
  CpuGridFields_Destroy();for(int c=0;c<4;c++)free(src[c]);
 }
 if(CpuGridFields_Configure(MPI_COMM_WORLD,rank==0?2:4)!=-1)errors++;
 if(CpuGridFields_Configure(MPI_COMM_WORLD,3)!=-1)errors++;
 int total;MPI_Allreduce(&errors,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
 if(rank==0)printf("UNIFIED_GRID_FIELDS widths=2,4 publication=all_fields scalar_duplicate=checked layout_bytes=8,14 rank_disagreement=rejected active_reconfigure=rejected errors=%d\n",total);
 MPI_Finalize();return total?1:0;
}
