#ifndef OPENMX_CPU_SECTOR_NATIVE_PACK_H
#define OPENMX_CPU_SECTOR_NATIVE_PACK_H
#include <mpi.h>
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
typedef struct {int n,ofs,nr,nc,nb;size_t offset;} CpuNativeSector;
typedef struct {int row,col;double value;} CpuNativePackEntry;
static int CpuNativePackOwner(int g,int nb,int p){return(g/nb)%p;}
static int CpuNativePackLocal(int g,int nb,int p){return(g/(nb*p))*nb+g%nb;}
static int CpuNativePackGlobal(int i,int nb,int coordinate,int p){return(i/nb)*nb*p+coordinate*nb+i%nb;}
/* Convert reset-origin sector descriptors into compact arrays with the
 * original full-matrix row/column owners. This is a value-preserving bounded
 * permutation, needed once after solves; source and destination are distinct. */
static int CpuNativePack(MPI_Comm world,int n,int nb,int nr,int nc,int pr,int pc,
 int prow,int pcol,const int*rankmap,const CpuNativeSector*sectors,int ns,
 const double*input,double*output,size_t capacity,size_t*offsets,size_t*used,
 int tile,size_t*temporary_bytes)
{
 int np;MPI_Comm_size(world,&np);int valid=n>0&&nb>0&&nr>=0&&nc>=0&&pr>0&&pc>0&&pr*pc==np&&
  prow>=0&&prow<pr&&pcol>=0&&pcol<pc&&rankmap&&sectors&&ns>0&&ns<=48&&
  input&&output&&input!=output&&offsets&&used&&tile>0&&tile<=1024;
 int rs[48],cs[48],rn[48],cn[48],racc=0,cacc=0,previous=0;size_t need=0;
 if(valid)for(int s=0;s<ns;s++){
  const CpuNativeSector*t=sectors+s;if(t->ofs!=previous||t->n<=0||t->ofs+t->n>n||t->nb<=0||t->nr<0||t->nc<0){valid=0;break;}
  rs[s]=racc;cs[s]=cacc;rn[s]=cn[s]=0;
  for(int g=t->ofs;g<t->ofs+t->n;g++){rn[s]+=CpuNativePackOwner(g,nb,pr)==prow;cn[s]+=CpuNativePackOwner(g,nb,pc)==pcol;}
  offsets[s]=need;need+=(size_t)rn[s]*cn[s];racc+=rn[s];cacc+=cn[s];previous+=t->n;
 }
 if(previous!=n||racc!=nr||cacc!=nc||need>capacity)valid=0;
 int all;MPI_Allreduce(&valid,&all,1,MPI_INT,MPI_MIN,world);if(!all)return -1;
 int*counts=calloc((size_t)np*5+2,sizeof(int));if(!counts){MPI_Abort(world,96);abort();}
 int*sc=counts,*rc=sc+np,*sd=rc+np,*rd=sd+np+1,*cursor=rd+np+1;
 CpuNativePackEntry*send=NULL,*recv=NULL;size_t sendcap=0,recvcap=0,peak=0,total_received=0;
 MPI_Datatype entry;MPI_Type_contiguous((int)sizeof(CpuNativePackEntry),MPI_BYTE,&entry);MPI_Type_commit(&entry);
 int bad=0;
 for(int begin=0;begin<n;begin+=tile){
  int end=n-begin<tile?n:begin+tile;for(int r=0;r<np;r++)sc[r]=0;
  for(int s=0;s<ns;s++){
   const CpuNativeSector*t=sectors+s;
   for(int j=0;j<t->nc;j++){
    int gj=t->ofs+CpuNativePackGlobal(j,t->nb,pcol,pc);if(gj<t->ofs||gj>=t->ofs+t->n||gj<begin||gj>=end)continue;
    for(int i=0;i<t->nr;i++){
     int gi=t->ofs+CpuNativePackGlobal(i,t->nb,prow,pr);if(gi<t->ofs||gi>=t->ofs+t->n)continue;
     int dest=rankmap[CpuNativePackOwner(gi,nb,pr)*pc+CpuNativePackOwner(gj,nb,pc)];
     if(sc[dest]==INT_MAX){MPI_Abort(world,96);abort();}sc[dest]++;
    }
   }
  }
  MPI_Alltoall(sc,1,MPI_INT,rc,1,MPI_INT,world);sd[0]=rd[0]=0;
  for(int r=0;r<np;r++){
   if((int64_t)sd[r]+sc[r]>INT_MAX||(int64_t)rd[r]+rc[r]>INT_MAX){MPI_Abort(world,96);abort();}
   sd[r+1]=sd[r]+sc[r];rd[r+1]=rd[r]+rc[r];cursor[r]=sd[r];
  }
  if((size_t)sd[np]>sendcap){free(send);sendcap=sd[np];send=malloc((sendcap?sendcap:1)*sizeof(*send));}
  if((size_t)rd[np]>recvcap){free(recv);recvcap=rd[np];recv=malloc((recvcap?recvcap:1)*sizeof(*recv));}
  if(!send){send=malloc(sizeof(*send));}if(!recv){recv=malloc(sizeof(*recv));}
  if(!send||!recv){MPI_Abort(world,96);abort();}
  size_t bytes=(sendcap+recvcap)*sizeof(*send)+(size_t)(np*5+2)*sizeof(int);if(bytes>peak)peak=bytes;
  for(int s=0;s<ns;s++){
   const CpuNativeSector*t=sectors+s;
   for(int j=0;j<t->nc;j++){
    int gj=t->ofs+CpuNativePackGlobal(j,t->nb,pcol,pc);if(gj<t->ofs||gj>=t->ofs+t->n||gj<begin||gj>=end)continue;
    for(int i=0;i<t->nr;i++){
     int gi=t->ofs+CpuNativePackGlobal(i,t->nb,prow,pr);if(gi<t->ofs||gi>=t->ofs+t->n)continue;
     int dest=rankmap[CpuNativePackOwner(gi,nb,pr)*pc+CpuNativePackOwner(gj,nb,pc)];
     send[cursor[dest]++]=(CpuNativePackEntry){gi,gj,input[t->offset+(size_t)j*t->nr+i]};
    }
   }
  }
  MPI_Alltoallv(send,sc,sd,entry,recv,rc,rd,entry,world);
  total_received+=(size_t)rd[np];
  for(int k=0;k<rd[np];k++){
   CpuNativePackEntry*v=recv+k;int s=0;while(s<ns&&v->row>=sectors[s].ofs+sectors[s].n)s++;
   if(s>=ns||v->row<sectors[s].ofs||v->col<sectors[s].ofs||v->col>=sectors[s].ofs+sectors[s].n||
      CpuNativePackOwner(v->row,nb,pr)!=prow||CpuNativePackOwner(v->col,nb,pc)!=pcol){bad=1;continue;}
   int i=CpuNativePackLocal(v->row,nb,pr)-rs[s],j=CpuNativePackLocal(v->col,nb,pc)-cs[s];
   if(i<0||i>=rn[s]||j<0||j>=cn[s]){bad=1;continue;}
   output[offsets[s]+(size_t)j*rn[s]+i]=v->value;
  }
 }
 if(total_received!=need)bad=1;
 MPI_Type_free(&entry);free(send);free(recv);free(counts);
 MPI_Allreduce(&bad,&all,1,MPI_INT,MPI_MAX,world);if(all)return -2;
 *used=need;if(temporary_bytes)*temporary_bytes=peak;return 0;
}
#endif
