#ifndef GSOC_EXPORT_AUDIT_H
#define GSOC_EXPORT_AUDIT_H
#include "../audit/gsoc_operator_gate.h"
/* Apply exactly the solver's sector restriction, conjugate-pair averaging and
 * real-sector projection to a Hermitian raw operator. This is independent of
 * real-space group averaging; comparing the two catches incomplete commutants.
 * O(nspin) scratch; callbacks remain collective. */
static int gsoc_export_projected_apply(const GsocAuditOperator*o,const int*partner,
 const GsocAuditComplex*x,GsocAuditComplex*h,GsocAuditComplex*s){
 int n=o->n,rc=0;size_t count=(size_t)n*9;
 if(n<1||count/9!=(size_t)n||!partner)return -1;
 GsocAuditComplex*v=calloc(count,sizeof(*v));if(!v)return -2;
 GsocAuditComplex*xb=v,*z=xb+n,*physical=z+n,*hp=physical+n,*sp=hp+n,*hb=sp+n,*sb=hb+n,*ho=sb+n,*so=ho+n;
 if((rc=o->inverse(o->ctx,x,xb)))goto done;
 for(int k=0;k<o->nsector;k++){
  int peer=partner[k];if(peer<0)for(int j=0;j<o->nsector;j++)if(partner[j]==k){peer=j;break;}
  int size=o->offset[k+1]-o->offset[k],choices=peer>=0?2:1,variants=o->dtype[k]==1?2:1;
  if(peer>=0&&o->offset[peer+1]-o->offset[peer]!=size){rc=-3;goto done;}
  for(int choice=0;choice<choices;choice++)for(int variant=0;variant<variants;variant++){
   int src=choice?peer:k,conjugate=choice^variant;double weight=1./(choices*variants);
   memset(z,0,(size_t)n*sizeof(*z));
   for(int j=0;j<size;j++){z[o->offset[src]+j]=xb[o->offset[k]+j];if(conjugate)z[o->offset[src]+j].i*=-1.;}
   if((rc=o->forward(o->ctx,z,physical))||(rc=o->apply(o->ctx,physical,hp,sp))||
      (rc=o->inverse(o->ctx,hp,hb))||(rc=o->inverse(o->ctx,sp,sb)))goto done;
   for(int j=0;j<size;j++){int from=o->offset[src]+j,to=o->offset[k]+j;double sign=conjugate?-1.:1.;
    ho[to].r+=weight*hb[from].r;ho[to].i+=weight*sign*hb[from].i;
    so[to].r+=weight*sb[from].r;so[to].i+=weight*sign*sb[from].i;}
  }
 }
 if((rc=o->forward(o->ctx,ho,h)))goto done;rc=o->forward(o->ctx,so,s);
 done:free(v);return rc;
}
typedef struct {const GsocAuditOperator*raw;const int*partner;} GsocProjectedAudit;
static int gsoc_export_projected_callback(void*ctx,const GsocAuditComplex*x,
 GsocAuditComplex*h,GsocAuditComplex*s){
 GsocProjectedAudit*p=(GsocProjectedAudit*)ctx;
 return gsoc_export_projected_apply(p->raw,p->partner,x,h,s);
}
#endif
