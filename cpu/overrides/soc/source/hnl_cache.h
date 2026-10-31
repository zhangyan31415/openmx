#ifndef OPENMX_EXACT_HNL_CACHE_H
#define OPENMX_EXACT_HNL_CACHE_H
/* Exact small angular operator applied to the right projector once per epoch.
 * Tables are generated with unchanged native one-hot probes. No SOC signs,
 * j-branch equivalence, positive energies, or orthogonal AO assumptions. */
#include <stdint.h>
#include <string.h>
#include <stddef.h>
#define HNL_CACHE_SLOTS 512
#define HNL_ANGULAR_MAX 7

typedef struct {double w[4][2][3][7][7];} HnlWeights;
typedef struct {
 uint64_t epoch;
 const double *right[2];
 double pf[2],energy[2];
 int ell;
 double q[2][4][7];
} HnlSlot;
typedef struct {
 HnlSlot slot[HNL_CACHE_SLOTS];
 uint64_t epoch,hits,misses,calls;
} HnlCache;
static void hnl_weights_init(HnlWeights*w,void(*oracle)(int,int,int,int,double[3])){
 memset(w,0,sizeof(*w));
 for(int ell=0;ell<4;ell++)for(int branch=0;branch<2;branch++)for(int a=0;a<=2*ell;a++)for(int b=0;b<=2*ell;b++){
  double out[3];oracle(ell,branch,a,b,out);for(int c=0;c<3;c++)w->w[ell][branch][c][a][b]=out[c];
 }
}
static void hnl_begin(HnlCache*c){
 if(++c->epoch==0){memset(c->slot,0,sizeof(c->slot));c->epoch=1;}
}
static uint64_t hnl_mix(uint64_t x){x^=x>>30;x*=UINT64_C(0xbf58476d1ce4e5b9);x^=x>>27;x*=UINT64_C(0x94d049bb133111eb);return x^(x>>31);}
static HnlSlot*hnl_right(HnlCache*c,const HnlWeights*w,int ell,double pf0,double pf1,double e0,double e1,const double*const b[2]){
 uint64_t hash=hnl_mix((uintptr_t)b[0])^hnl_mix((uintptr_t)b[1]+(unsigned)ell);
 HnlSlot*s=&c->slot[hash&(HNL_CACHE_SLOTS-1)];c->calls++;
 if(s->epoch==c->epoch&&s->right[0]==b[0]&&s->right[1]==b[1]&&s->ell==ell&&s->pf[0]==pf0&&s->pf[1]==pf1&&s->energy[0]==e0&&s->energy[1]==e1){c->hits++;return s;}
 c->misses++;s->epoch=c->epoch;s->right[0]=b[0];s->right[1]=b[1];s->ell=ell;s->pf[0]=pf0;s->pf[1]=pf1;s->energy[0]=e0;s->energy[1]=e1;
 int n=2*ell+1;
 for(int branch=0;branch<2;branch++)for(int i=0;i<n;i++){
  s->q[branch][0][i]=s->pf[branch]*s->energy[branch]*b[branch][i];
  for(int comp=0;comp<3;comp++){
   double x=0.;for(int j=0;j<n;j++)x+=w->w[ell][branch][comp][i][j]*b[branch][j];
   s->q[branch][comp+1][i]=s->energy[branch]*x;
  }
 }
 return s;
}
static void hnl_apply(HnlCache*c,const HnlWeights*w,int ell,double pf0,double pf1,double e0,double e1,double fugou,const double*a[2][3],const double*b[2],double sums[18]){
 HnlSlot*s=hnl_right(c,w,ell,pf0,pf1,e0,e1,b);int n=2*ell+1;
 for(int d=0;d<3;d++){
  double rr=0.,ii=0.,ur=0.,ui=0.;
  for(int branch=0;branch<2;branch++){
   const double *x=a[branch][d],*q0=s->q[branch][0],*q1=s->q[branch][1],*q2=s->q[branch][2],*q3=s->q[branch][3];
   for(int i=0;i<n;i++){rr+=x[i]*q0[i];ii+=x[i]*q1[i];ur+=x[i]*q2[i];ui+=x[i]*q3[i];}
  }
  sums[d]+=rr;sums[3+d]+=rr;
  sums[9+d]+=fugou*ii;sums[12+d]-=fugou*ii;
  sums[6+d]+=fugou*ur;sums[15+d]+=fugou*ui;
 }
}
#endif
