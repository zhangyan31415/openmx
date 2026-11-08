#define _POSIX_C_SOURCE 200809L
#include "cpu_h_grid_blas3.h"
#include <mkl.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
static double stamp(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
int CpuHGridFlag(const char*s){if(!s||!strcmp(s,"0"))return 0;return !strcmp(s,"1")?1:-1;}
int CpuHGridInit(CpuHGridWorkspace*w,int ao,int tile){
 if(!w)return -1;memset(w,0,sizeof(*w));if(ao<1||tile<1)return -1;
 if((size_t)ao>(1024*1024)/(2*sizeof(double))/(size_t)tile)return -2;
 size_t bytes=2*(size_t)ao*(size_t)tile*sizeof(double);void*p=0;
 if(posix_memalign(&p,64,bytes))return -3;
 w->ao=ao;w->tile=tile;w->bytes=bytes;w->left=p;w->right=w->left+(size_t)ao*tile;return 0;
}
void CpuHGridDestroy(CpuHGridWorkspace*w){if(w){free(w->left);memset(w,0,sizeof(*w));}}
int CpuHGridSingleThreadBegin(void){return mkl_set_num_threads_local(1);}
void CpuHGridSingleThreadEnd(int old){mkl_set_num_threads_local(old);}
int CpuHGridAccumulate(const CpuHGridView*v,CpuHGridWorkspace*w,double*h,CpuHGridTiming*tm){
 if(!v||!w||!w->left||!w->right||!h||v->m<1||v->n<1||v->m>w->ao||v->n>w->ao||v->ldh<v->n||v->count<0||w->tile<1||(v->b_per_overlap!=0&&v->b_per_overlap!=1))return -1;
 if(!v->count)return 0;
 if(!v->a||(!v->b&&!v->overlap_row)||((!v->nc||!v->nh)&&!v->overlap_index)||!v->legacy_grid||!v->global_grid||!v->potential||!v->index)return -1;
 for(int base=0;base<v->count;){
  int k=CpuHGridTileLength(v->count,base,w->tile);double t=tm?stamp():0;
  for(int q=0;q<k;q++){
   int l=base+q,c=v->overlap_index?v->overlap_index(0,v->center,v->neighbor,l):v->nc[l];size_t p=v->index((size_t)v->legacy_grid[c],(unsigned long long)v->global_grid[c]);
   int nh=v->overlap_index?v->overlap_index(1,v->center,v->neighbor,l):v->nh[l];
   const float*a=v->a[c],*b=(v->b_per_overlap&&v->overlap_row)?v->overlap_row(v->center,v->neighbor,l):v->b[v->b_per_overlap?l:nh];double gv=v->gridvol*v->potential[p];
   for(int i=0;i<v->m;i++)w->left[(size_t)q*v->m+i]=gv*(double)a[i];
   for(int j=0;j<v->n;j++)w->right[(size_t)q*v->n+j]=(double)b[j];
  }
  if(tm){tm->pack_seconds+=stamp()-t;t=stamp();}
  cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,v->m,v->n,k,1,w->left,v->m,w->right,v->n,1,h,v->ldh);
  if(tm)tm->gemm_seconds+=stamp()-t;
  base+=k;/* k<=count-base, so this cannot overflow an int. */
 }
 return 0;
}
