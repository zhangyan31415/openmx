#ifndef CPU_NSOC_ELPA_OMP_SCOPE_H
#define CPU_NSOC_ELPA_OMP_SCOPE_H
#include <omp.h>
#include <mkl.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
typedef struct {int before_omp,before_mkl,inside_omp,inside_mkl,restored_omp,restored_mkl;} CpuNSocElpaOmpScope;
static int CpuNSocElpaOmpParse(const char*text,int*setting){
 *setting=0;if(!text||!*text)return 1;char*end;double value=strtod(text,&end);
 if(end==text||*end||!isfinite(value)||(value!=0&&value!=1))return 0;*setting=(int)value;return 1;
}
static int CpuNSocElpaOmpBegin(CpuNSocElpaOmpScope*scope){
 memset(scope,0,sizeof(*scope));scope->before_omp=omp_get_max_threads();scope->before_mkl=mkl_get_max_threads();
 if(scope->before_omp!=4||scope->before_mkl!=4||omp_in_parallel())return 0;
 omp_set_num_threads(1);scope->inside_omp=omp_get_max_threads();scope->inside_mkl=mkl_get_max_threads();
 if(scope->inside_omp!=1||scope->inside_mkl!=4){omp_set_num_threads(scope->before_omp);return 0;}
 return 1;
}
static int CpuNSocElpaOmpEnd(CpuNSocElpaOmpScope*scope){
 int inside_omp=omp_get_max_threads(),inside_mkl=mkl_get_max_threads();
 omp_set_num_threads(scope->before_omp);scope->restored_omp=omp_get_max_threads();scope->restored_mkl=mkl_get_max_threads();
 return inside_omp==1&&inside_mkl==4&&scope->restored_omp==scope->before_omp&&scope->restored_mkl==scope->before_mkl;
}
#endif
