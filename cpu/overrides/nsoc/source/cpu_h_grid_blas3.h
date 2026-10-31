#ifndef CPU_H_GRID_BLAS3_H
#define CPU_H_GRID_BLAS3_H
#include <stddef.h>
typedef struct {
 int m,n,count,ldh,b_per_overlap;
 const float *const *a,*const *b;
 const int *nc,*nh,*legacy_grid,*global_grid;
 const double *potential;
 double gridvol;
 size_t (*index)(size_t,unsigned long long);
} CpuHGridView;
typedef struct {int ao,tile;size_t bytes;double *left,*right;} CpuHGridWorkspace;
typedef struct {double pack_seconds,gemm_seconds;} CpuHGridTiming;
/* Strict optional0/1; absent meansoff, malformed means-1. */
int CpuHGridFlag(const char*);
int CpuHGridInit(CpuHGridWorkspace*,int,int);
void CpuHGridDestroy(CpuHGridWorkspace*);
int CpuHGridSingleThreadBegin(void);
void CpuHGridSingleThreadEnd(int);
int CpuHGridAccumulate(const CpuHGridView*,CpuHGridWorkspace*,double*,CpuHGridTiming*);
/* Remaining-length arithmetic is safe even when total==INT_MAX. */
static inline int CpuHGridTileLength(int total,int base,int tile){
 if(total<0||base<0||base>total||tile<1)return -1;
 int remain=total-base;return remain<tile?remain:tile;
}
#endif
