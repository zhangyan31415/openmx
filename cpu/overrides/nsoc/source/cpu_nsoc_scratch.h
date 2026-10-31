#ifndef OPENMX_CPU_NSOC_SCRATCH_H
#define OPENMX_CPU_NSOC_SCRATCH_H
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
/* The caller owns the allocation. This descriptor owns no storage. Its alias
 * may be used only while previous native eigenvectors are dead, ending before
 * native redistribution writes fresh AO-major coefficients. */
typedef struct { double *ptr; size_t capacity; } CpuNSocEigenScratch;
static int CpuNSocScratch_Reserve(CpuNSocEigenScratch *s,double **p,size_t current,size_t required){
 if(!s||!p||!*p||!current||!required||current>SIZE_MAX/sizeof(double)||required>SIZE_MAX/sizeof(double))return -1;
 if(required>current){double *q=realloc(*p,required*sizeof(double));if(!q)return -1;*p=q;current=required;}
 s->ptr=*p;s->capacity=current;return 0;
}
static double *CpuNSocScratch_Borrow(CpuNSocEigenScratch *s,size_t required){
 return s&&s->ptr&&required&&required<=s->capacity?s->ptr:NULL;
}
static void CpuNSocScratch_Release(CpuNSocEigenScratch *s){if(s){s->ptr=NULL;s->capacity=0;}}
#endif
