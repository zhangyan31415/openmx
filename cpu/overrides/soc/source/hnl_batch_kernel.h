#ifndef OPENMX_HNL_BATCH_KERNEL_H
#define OPENMX_HNL_BATCH_KERNEL_H
#include <mkl.h>
#include "hnl_cache.h"
/* A stacks three derivative directions and concatenates both distinct j
 * branches. B stacks four native spin components. One FP64 GEMM replaces
 * AO-pair/radial-channel scalar calls, with original angular weights. */
static void hnl_batch_compute(const HnlWeights*w,int nl,int nr,int K,int radial,
 const int*ell,const double*ep,const double*em,double esign,
 double**left[2][3],double**right[2],double*A,double*B,double*C){
 const int rows=3*nl,cols=4*nr,inner=2*K;
 for(int d=0;d<3;d++)for(int m=0;m<nl;m++)for(int b=0;b<2;b++)
  memcpy(A+(size_t)(d*nl+m)*inner+b*K,left[b][d][m],(size_t)K*sizeof(double));
 int offset=0;
 for(int r=0;r<radial;r++){
  const int L=ell[r],nang=2*L+1;double energy[2]={esign*ep[r],esign*em[r]};double pf[2]={(L+1.)/(2*L+1.),L/(2*L+1.)};
  for(int branch=0;branch<2;branch++)for(int i=0;i<nang;i++){
   double *row=B+(size_t)(branch*K+offset+i)*cols;
   for(int n=0;n<nr;n++){
    const double*x=right[branch][n]+offset;
    row[n]=pf[branch]*energy[branch]*x[i];
    for(int c=0;c<3;c++){
     double q=0.;for(int j=0;j<nang;j++)q+=w->w[L][branch][c][i][j]*x[j];
     row[(c+1)*nr+n]=energy[branch]*q;
    }
   }
  }
  offset+=nang;
 }
 /* The caller validates offset==K and dimensions before this kernel. */
 int previous=mkl_set_num_threads_local(1);
 cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,rows,cols,inner,1.,A,inner,B,cols,0.,C,cols);
 mkl_set_num_threads_local(previous);
}
static void hnl_batch_value(const double*C,int nl,int nr,int m,int n,double fugou,double out[18]){
 for(int d=0;d<3;d++){
  const double*row=C+(size_t)(d*nl+m)*4*nr;
  out[d]=out[3+d]=row[n];out[9+d]=fugou*row[nr+n];out[12+d]=-out[9+d];
  out[6+d]=fugou*row[2*nr+n];out[15+d]=fugou*row[3*nr+n];
 }
}
#endif
