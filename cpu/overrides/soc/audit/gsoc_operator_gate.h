#ifndef GSOC_OPERATOR_GATE_H
#define GSOC_OPERATOR_GATE_H
/* Independent C99, allocation-bounded audit for a distributed native SOC solver.
 * Callbacks are collective. Complex layout is two consecutive doubles r,i.
 * This is a sampled operator guard, NOT an entrywise covariance proof.
 * For raw covariance/Hermiticity diagnostics, apply uses RAW directed-neighbor
 * H and matching S. Eigenpair audits may instead use the solved projection.
 * forward/inverse apply the complete unitary CSC U/U^dagger in native spin order.
 * Runtime memory O(nspin) plus O(nspin * sampled_bands) for final eigen checks.
 */
#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
typedef struct { double r,i; } GsocAuditComplex;
typedef int (*GsocAuditApply)(void *,const GsocAuditComplex *,GsocAuditComplex *,GsocAuditComplex *);
typedef int (*GsocAuditBasis)(void *,const GsocAuditComplex *,GsocAuditComplex *);
typedef int (*GsocAuditEigenvector)(void *,int,GsocAuditComplex *,double *);
typedef struct {
 int n,nsector; const int *offset,*dtype;
 void *ctx; GsocAuditApply apply; GsocAuditBasis forward,inverse;
} GsocAuditOperator;
typedef struct {
 double comm_abs[2],comm_rel[2],tr_abs[2],tr_rel[2],real_abs[2],real_rel[2];
 int probes,nonfinite; uint64_t seed;
} GsocAuditGateResult;
/* Spatial sector/real-gauge defects measure the explicitly requested projection.
 * Nonmagnetic time reversal is still required; spatial grid discretization is
 * not a reason to discard a magnetic term. Nonfinite measurements always fail. */
static int gsoc_audit_gate_status(const GsocAuditGateResult*r,int projected,double tolerance){
 if(!r||r->nonfinite||!isfinite(tolerance)||tolerance<=0)return -1;
 for(int field=0;field<2;field++){
  if(!isfinite(r->comm_rel[field])||!isfinite(r->tr_rel[field])||!isfinite(r->real_rel[field]))return -1;
  if(r->tr_rel[field]>tolerance || (!projected &&
     (r->comm_rel[field]>tolerance||r->real_rel[field]>tolerance)))return 1;
 }
 return 0;
}
typedef struct {double raw_abs,raw_relative,Sorth;int count,nonfinite;} GsocAuditEigenResult;
static GsocAuditComplex gsoc_audit_c(double r,double i){GsocAuditComplex z;z.r=r;z.i=i;return z;}
static GsocAuditComplex gsoc_audit_mul(GsocAuditComplex a,GsocAuditComplex b){return gsoc_audit_c(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);}
static GsocAuditComplex gsoc_audit_conj(GsocAuditComplex a){a.i=-a.i;return a;}
static double gsoc_audit_norm(const GsocAuditComplex*x,int n){double s=0.;int i;for(i=0;i<n;i++)s+=x[i].r*x[i].r+x[i].i*x[i].i;return sqrt(s);}
static uint64_t gsoc_audit_rng(uint64_t*s){uint64_t x=*s;x^=x>>12;x^=x<<25;x^=x>>27;*s=x;return x*UINT64_C(2685821657736338717);}
/* Physical time reversal J K, J=i sigma_y, native order up[0:N],down[N:2N]. */
static void gsoc_audit_time_reverse(int n,const GsocAuditComplex*x,GsocAuditComplex*y){int i,N=n/2;for(i=0;i<N;i++){y[i]=gsoc_audit_conj(x[N+i]);y[N+i]=gsoc_audit_c(-x[i].r,x[i].i);}}
static void gsoc_audit_metric(const GsocAuditComplex*a,const GsocAuditComplex*b,int n,double*ab,double*rel,int*nonfinite){int i;double ss=0,aa=0,bb=0;for(i=0;i<n;i++){double r=a[i].r-b[i].r,t=a[i].i-b[i].i;ss+=r*r+t*t;aa+=a[i].r*a[i].r+a[i].i*a[i].i;bb+=b[i].r*b[i].r+b[i].i*b[i].i;if(!isfinite(ss)||!isfinite(aa)||!isfinite(bb))*nonfinite=1;}double d=sqrt(ss),scale=fmax(sqrt(aa),sqrt(bb));if(d>*ab)*ab=d;if(d/fmax(scale,1e-300)>*rel)*rel=d/fmax(scale,1e-300);}
static int gsoc_audit_project(const GsocAuditOperator*o,int sector,const GsocAuditComplex*x,GsocAuditComplex*y,GsocAuditComplex*t){int i;if(o->inverse(o->ctx,x,t))return -2;for(i=0;i<o->n;i++)if(i<o->offset[sector]||i>=o->offset[sector+1])t[i]=gsoc_audit_c(0,0);return o->forward(o->ctx,t,y);}
/* Returns 0 for successful measurement, negative for callback/schema/allocation
 * failure. Caller separately applies predeclared tolerances; there is no hidden
 * projection, relaxed default, or automatic fallback to an unsymmetrized solve.
 */
static int gsoc_audit_operator_gate(const GsocAuditOperator*o,int probes,uint64_t seed,GsocAuditGateResult*r){
 int p,s,k,i,n=o->n,rc=0;GsocAuditComplex *v,*x,*h,*ov,*px,*hp,*sp,*expect,*t,*q;
 if(!r||n<=0||n%2||o->nsector<1||!o->offset||!o->dtype||!o->apply||!o->forward||!o->inverse||probes<1||!seed)return -1;
 if(o->offset[0]!=0||o->offset[o->nsector]!=n)return -1;
 for(s=0;s<o->nsector;s++)if(o->offset[s]>=o->offset[s+1]||(o->dtype[s]!=1&&o->dtype[s]!=2))return -1;
 memset(r,0,sizeof(*r));r->probes=probes;r->seed=seed;
 v=(GsocAuditComplex*)calloc((size_t)n*9,sizeof(*v));if(!v)return -3;
 x=v;h=x+n;ov=h+n;px=ov+n;hp=px+n;sp=hp+n;expect=sp+n;t=expect+n;q=t+n;
 for(p=0;p<probes;p++){
  double scale=1./sqrt(2.*n);for(i=0;i<n;i++)x[i]=gsoc_audit_c((gsoc_audit_rng(&seed)&1)?scale:-scale,(gsoc_audit_rng(&seed)&1)?scale:-scale);
  if((rc=o->apply(o->ctx,x,h,ov)))goto done;
  for(s=0;s<o->nsector;s++){
   if((rc=gsoc_audit_project(o,s,x,px,t))||(rc=o->apply(o->ctx,px,hp,sp)))goto done;
   for(k=0;k<2;k++){GsocAuditComplex*a=k?ov:h,*b=k?sp:hp;if((rc=gsoc_audit_project(o,s,a,expect,t)))goto done;gsoc_audit_metric(b,expect,n,&r->comm_abs[k],&r->comm_rel[k],&r->nonfinite);}
  }
  gsoc_audit_time_reverse(n,x,px);if((rc=o->apply(o->ctx,px,hp,sp)))goto done;
  for(k=0;k<2;k++){gsoc_audit_time_reverse(n,k?ov:h,expect);gsoc_audit_metric(k?sp:hp,expect,n,&r->tr_abs[k],&r->tr_rel[k],&r->nonfinite);}
  for(s=0;s<o->nsector;s++)if(o->dtype[s]==1){
   memset(t,0,(size_t)n*sizeof(*t));scale=1./sqrt((double)(o->offset[s+1]-o->offset[s]));for(i=o->offset[s];i<o->offset[s+1];i++)t[i].r=(gsoc_audit_rng(&seed)&1)?scale:-scale;
   if((rc=o->forward(o->ctx,t,px))||(rc=o->apply(o->ctx,px,hp,sp)))goto done;
   for(k=0;k<2;k++){double imag=0,all=0;if((rc=o->inverse(o->ctx,k?sp:hp,q)))goto done;for(i=o->offset[s];i<o->offset[s+1];i++){imag+=q[i].i*q[i].i;all+=q[i].r*q[i].r+q[i].i*q[i].i;}if(!isfinite(imag)||!isfinite(all))r->nonfinite=1;imag=sqrt(imag);if(imag>r->real_abs[k])r->real_abs[k]=imag;if(imag/fmax(sqrt(all),1e-300)>r->real_rel[k])r->real_rel[k]=imag/fmax(sqrt(all),1e-300);}
  }
 }
 done:free(v);return rc? -2:0;
}
/* Residual checks use the operator supplied by the caller: raw for diagnostics,
 * or exactly the solved projected H/S for projected-solver acceptance. The raw_*
 * result member names are retained for API compatibility. S orthogonality
 * includes distinct sampled Kramers partners; getvec returns physical spinors.
 */
static int gsoc_audit_eigenpairs(const GsocAuditOperator*o,GsocAuditEigenvector getvec,const int*bands,int count,GsocAuditEigenResult*r){
 int a,b,i,n=o->n,rc=0;double e;GsocAuditComplex *c,*sc,*hc,*tmp;
 if(!r||!getvec||!bands||count<1||n<=0)return -1;
 memset(r,0,sizeof(*r));r->count=count;
 c=(GsocAuditComplex*)calloc((size_t)n*(2*count+2),sizeof(*c));if(!c)return -3;sc=c+(size_t)n*count;hc=sc+(size_t)n*count;tmp=hc+n;
 for(a=0;a<count;a++){
  GsocAuditComplex*x=c+(size_t)a*n,*sx=sc+(size_t)a*n;
  if((rc=getvec(o->ctx,bands[a],x,&e))||(rc=o->apply(o->ctx,x,hc,sx)))goto done;
  double hh=gsoc_audit_norm(hc,n),ss=gsoc_audit_norm(sx,n),rs=0;
  for(i=0;i<n;i++){tmp[i]=gsoc_audit_c(hc[i].r-e*sx[i].r,hc[i].i-e*sx[i].i);double z=hypot(tmp[i].r,tmp[i].i);if(z>r->raw_abs)r->raw_abs=z;rs+=z*z;}
  double rel=sqrt(rs)/fmax(hh+fabs(e)*ss,1e-300);if(rel>r->raw_relative)r->raw_relative=rel;
  if(!isfinite(e)||!isfinite(rs)||!isfinite(rel)||!isfinite(hh)||!isfinite(ss))r->nonfinite=1;
 }
 for(a=0;a<count;a++)for(b=0;b<count;b++){double re=0,im=0;for(i=0;i<n;i++){GsocAuditComplex x=c[(size_t)a*n+i],y=sc[(size_t)b*n+i];re+=x.r*y.r+x.i*y.i;im+=x.r*y.i-x.i*y.r;}double d=hypot(re-(a==b),im);if(d>r->Sorth)r->Sorth=d;if(!isfinite(d))r->nonfinite=1;}
 done:free(c);return rc?-2:0;
}
/* Independent Hermiticity bilinear probe. Sector/TR/real-gauge checks alone
 * cannot detect a real nonsymmetric matrix inside an invariant sector. */
typedef struct {double absolute[2],relative[2];int nonfinite,probes;} GsocAuditHermitianResult;
static int gsoc_audit_hermiticity(const GsocAuditOperator*o,int probes,uint64_t seed,GsocAuditHermitianResult*r){
 if(!o||!r||o->n<1||!o->apply||probes<1||!seed)return -1;
 int n=o->n,p,i,k,rc=0;memset(r,0,sizeof(*r));r->probes=probes;
 GsocAuditComplex*v=(GsocAuditComplex*)calloc((size_t)n*6,sizeof(*v));if(!v)return -3;
 GsocAuditComplex*x=v,*y=x+n,*hx=y+n,*sx=hx+n,*hy=sx+n,*sy=hy+n;
 double scale=1./sqrt(2.*n);
 for(p=0;p<probes;p++){
  for(i=0;i<n;i++){x[i]=gsoc_audit_c((gsoc_audit_rng(&seed)&1)?scale:-scale,(gsoc_audit_rng(&seed)&1)?scale:-scale);y[i]=gsoc_audit_c((gsoc_audit_rng(&seed)&1)?scale:-scale,(gsoc_audit_rng(&seed)&1)?scale:-scale);}
  if((rc=o->apply(o->ctx,x,hx,sx))||(rc=o->apply(o->ctx,y,hy,sy)))break;
  for(k=0;k<2;k++){GsocAuditComplex *ax=k?sx:hx,*ay=k?sy:hy;double re=0,im=0;
   for(i=0;i<n;i++){re+=x[i].r*ay[i].r+x[i].i*ay[i].i-ax[i].r*y[i].r-ax[i].i*y[i].i;im+=x[i].r*ay[i].i-x[i].i*ay[i].r-ax[i].r*y[i].i+ax[i].i*y[i].r;}
   double a=hypot(re,im),den=gsoc_audit_norm(x,n)*gsoc_audit_norm(ay,n)+gsoc_audit_norm(ax,n)*gsoc_audit_norm(y,n),b=a/fmax(den,1e-300);
   if(!isfinite(a)||!isfinite(den)||!isfinite(b))r->nonfinite=1;
   if(a>r->absolute[k])r->absolute[k]=a;
   if(b>r->relative[k])r->relative[k]=b;
  }
 }
 free(v);return rc?-2:0;
}
#endif
