#ifndef OPENMX_CFS_SELECTED_LEFT_H
#define OPENMX_CFS_SELECTED_LEFT_H
#include <mpi.h>
#include <stdint.h>
#include <stddef.h>
typedef struct CfsSelectedPlan CfsSelectedPlan;
typedef struct {double pack,collective,multiply;unsigned long long logical_values;int columns,batches;size_t plan_bytes;} CfsSelectedTimes;
CfsSelectedPlan*CfsSelectedCreate(int n,int nb,int nr,int nc,int pr,int pc,int row,int col,const int64_t*cp,const int*ri,const double*cv,MPI_Comm comm,char*error,size_t ez);
int CfsSelectedLeft(CfsSelectedPlan*p,const double*A,int lda,double*W,int ldw,const int*columns,int count,int profile,CfsSelectedTimes*t,char*error,size_t ez);
void CfsSelectedDestroy(CfsSelectedPlan*p);
#endif
