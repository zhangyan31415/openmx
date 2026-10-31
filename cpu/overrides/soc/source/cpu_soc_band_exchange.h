#ifndef OPENMX_CPU_SOC_BAND_EXCHANGE_H
#define OPENMX_CPU_SOC_BAND_EXCHANGE_H
#include <mpi.h>
#include <stdint.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
/* Transport complete complex band vectors in bounded contiguous packets.
 * local[dp[dest]+j] identifies the source band; receive rows remain grouped by
 * source rank. No byte displacement is passed to MPI. ev must stay unchanged
 * until every packet has arrived in the existing receive workspace. */
static int CpuSocExchangeBands(int ns,const void *ev,const int *local,
 const int *ct,const int *dp,const int *rc,const int *rdp,void *recv,
 size_t limit_bytes,MPI_Comm comm)
{
 int rank,p;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&p);
 if(ns<1||(size_t)ns>SIZE_MAX/(2*sizeof(double)))return -1;
 size_t bytes=(size_t)ns*2*sizeof(double),capacity=limit_bytes/bytes;
 size_t count_bound=(size_t)INT_MAX/((size_t)ns*2);
 if(capacity>count_bound)capacity=count_bound;
 if(!capacity||capacity>INT_MAX)return -1;
 unsigned char *send=(unsigned char*)malloc(capacity*bytes);
 int failed=send==NULL,any_failed;
 MPI_Allreduce(&failed,&any_failed,1,MPI_INT,MPI_MAX,comm);
 if(any_failed){free(send);return -2;}
 for(int shift=0;shift<p;shift++){
  int dest=(rank+shift)%p,src=(rank-shift+p)%p;
  if(ct[dest]<0||rc[src]<0||dp[dest]<0||rdp[src]<0){free(send);return -1;}
  int nsends=(int)(((size_t)ct[dest]+capacity-1)/capacity);
  int nrecvs=(int)(((size_t)rc[src]+capacity-1)/capacity);
  int rounds=nsends>nrecvs?nsends:nrecvs,global_rounds;
  MPI_Allreduce(&rounds,&global_rounds,1,MPI_INT,MPI_MAX,comm);
  for(int round=0;round<global_rounds;round++){
   size_t base=(size_t)round*capacity;
   int sn=base<(size_t)ct[dest]?(int)((size_t)ct[dest]-base):0;
   int rn=base<(size_t)rc[src]?(int)((size_t)rc[src]-base):0;
   if((size_t)sn>capacity)sn=(int)capacity;if((size_t)rn>capacity)rn=(int)capacity;
   for(int band=0;band<sn;band++){
    int at=local[(size_t)dp[dest]+base+band];if(at<0){free(send);return -1;}
    memcpy(send+(size_t)band*bytes,(const unsigned char*)ev+(size_t)at*bytes,bytes);
   }
   /* For a zero-count receive, pass the workspace base rather than forming a
    * pointer beyond its object during another rank's additional rounds. */
   void *target=rn?(unsigned char*)recv+((size_t)rdp[src]+base)*bytes:recv;
   int result=MPI_Sendrecv(send,sn*ns*2,MPI_DOUBLE,dest,17931,
                         target,rn*ns*2,MPI_DOUBLE,src,17931,comm,MPI_STATUS_IGNORE);
   if(result){free(send);return -3;}
  }
 }
 free(send);return 0;
}
#endif
