/* Fixed scalar/spinor VNA force: one native projector-gradient pair at a time.
 * Keep the original FP32 projector products and FP64 density/force sums. */
#include <stdio.h>
#include "openmx_common.h"
#include "cpu_vna_shared.h"
#include "cpu_vna_stream.h"
#include "cpu_overlap_index.h"
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static void stream_fail(const char *message){int rank;MPI_Comm_rank(mpi_comm_level1,&rank);fprintf(stderr,"CPU_VNA_FORCE_STREAM_ERROR rank=%d %s\n",rank,message);MPI_Abort(mpi_comm_level1,96);abort();}
int CpuVnaForceStream_Mode(void){
 const char *s=getenv("OPENMX_CPU_VNA_FORCE_STREAM");int mode=0;
 if(s&&*s){if(!strcmp(s,"1"))mode=1;else if(!strcmp(s,"2"))mode=2;else if(strcmp(s,"0"))stream_fail("mode must be 0(off),1(keep raw),2(drop then restore raw)");}
 int lo,hi;MPI_Allreduce(&mode,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(&mode,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);if(lo!=hi)stream_fail("mode differs across ranks");
 int spin_eligible=(SpinP_switch==0&&SO_switch==0)||CpuSOCMemoryEligible();
 int eligible=CpuVnaShared_Deferred()&&MD_switch==0&&Solver==2&&spin_eligible&&Cnt_switch==0&&ProExpn_VNA==1&&!scf_stress_flag&&omp_get_max_threads()==1;
 int all;MPI_Allreduce(&eligible,&all,1,MPI_INT,MPI_MIN,mpi_comm_level1);return all?mode:0;
}
typedef struct {double *****d;double *force;size_t pairs;int projectors;} StreamForce;
/* VNA is a spin-independent scalar operator. Native SOC Force4B contracts
 * its derivative with D_up,up + D_down,down; the off-diagonal spin blocks do
 * not contribute. Half the trace retains the scalar path's native factors. */
static double stream_density(StreamForce *c,int local,int slot,int i,int j){
 return SpinP_switch==0?c->d[0][local][slot][i][j]:
   .5*(c->d[0][local][slot][i][j]+c->d[1][local][slot][i][j]);
}
/* Each ordered offsite AO pair contributes its basis derivative, and the
 * opposite projector-center derivative. Native on-center replacements are
 * handled separately below. Periodic images retain their neighbor slots. */
static void consume_pair(void *opaque,int local,int qslot,float *gradient,int ni,int p){
 StreamForce *c=opaque;int a=M2G[local],q=natn[a][qslot],species=WhatSpecies[q];
 if(p!=c->projectors||ni!=Spe_Total_CNO[WhatSpecies[a]])stream_fail("pair shape differs from native basis");
 size_t plane=(size_t)ni*p;int l=0;
 for(int angular=0;angular<=List_YOUSO[35];angular++)for(int radial=0;radial<List_YOUSO[34];radial++){
  double e=VNA_proj_ene[species][angular][radial];for(int m=-angular;m<=angular;m++,l++)for(int i=0;i<ni;i++)for(int xyz=0;xyz<3;xyz++)gradient[(size_t)xyz*plane+(size_t)i*p+l]=(float)(e*gradient[(size_t)xyz*plane+(size_t)i*p+l]);
 }
 if(l!=p)stream_fail("projector enumeration mismatch");
 for(int jslot=1;jslot<=FNAN[a];jslot++){
  int rev=RMI1[local][jslot][qslot];if(rev<0)continue;int b=natn[a][jslot],nj=Spe_Total_CNO[WhatSpecies[b]];
  double damp=dampingF(Spe_Atom_Cut1[WhatSpecies[a]]+Spe_Atom_Cut1[WhatSpecies[b]],Dis[a][jslot]);double sum[3]={0.,0.,0.};
  for(int i=0;i<ni;i++)for(int j=0;j<nj;j++){
   const float *raw=CpuVnaShared_Row(0,b,rev,j);const float *gx=gradient+(size_t)i*p,*gy=gx+plane,*gz=gy+plane;double x=0.,y=0.,z=0.;
   for(int k=0;k<p;k++){x+=(double)(gx[k]*raw[k]);y+=(double)(gy[k]*raw[k]);z+=(double)(gz[k]*raw[k]);}
   double scale=4.*stream_density(c,local,jslot,i,j);sum[0]+=scale*x;sum[1]+=scale*y;sum[2]+=scale*z;
  }
  for(int xyz=0;xyz<3;xyz++){double v=damp*sum[xyz];c->force[3*(size_t)a+xyz]+=v;c->force[3*(size_t)q+xyz]-=v;}
 }
 c->pairs++;
}
static void native_extras(StreamForce *c){
 for(int local=1;local<=Matomnum;local++){
  int a=M2G[local],na=Spe_Total_CNO[WhatSpecies[a]];
  /* Native h=q=0: exact one-center-orbital/two-center VNA derivative. */
  double self_damp=dampingF(2.*Spe_Atom_Cut1[WhatSpecies[a]],Dis[a][0]);
  for(int i=0;i<na;i++)for(int j=0;j<na;j++)for(int xyz=0;xyz<3;xyz++){
   double h=0.;for(int q=1;q<=FNAN[a];q++)h+=HVNA2[xyz+1][local][q][i][j];c->force[3*(size_t)a+xyz]+=2.*stream_density(c,local,0,i,j)*(self_damp*h);
  }
  for(int q=1;q<=FNAN[a];q++){
   int b=natn[a][q],lb=F_G2M[b],nb=Spe_Total_CNO[WhatSpecies[b]],pair=RMI1[local][0][q];
   /* Native second case h=q!=0: the matching exact HVNA3 term. */
   double onsite_damp=dampingF(2.*Spe_Atom_Cut1[WhatSpecies[b]],Dis[b][0]);
   for(int i=0;i<nb;i++)for(int j=0;j<nb;j++)for(int xyz=0;xyz<3;xyz++)c->force[3*(size_t)a+xyz]+=2.*stream_density(c,lb,0,i,j)*(-onsite_damp*HVNA3[xyz+1][local][q][i][j]);
   /* Both native cases include the offsite damping derivative. Preserve
    * those two additions explicitly rather than changing that convention. */
   double cut=Spe_Atom_Cut1[WhatSpecies[a]]+Spe_Atom_Cut1[WhatSpecies[b]],r=Dis[a][q],scale=0.;
   if(cut>r){double d=dampingF(cut,r);scale=deri_dampingF(cut,r)/d;}if(r<1.e-10)r=1.e-10;
   for(int pass=0;pass<2;pass++)for(int i=0;i<na;i++)for(int j=0;j<nb;j++)for(int xyz=0;xyz<3;xyz++){
    double direction=(Gxyz[a][xyz+1]+atv[ncn[a][0]][xyz+1]-Gxyz[b][xyz+1]-atv[ncn[a][q]][xyz+1])/r;
    double h=HVNA[local][q][i][j]*(scale*direction);c->force[3*(size_t)a+xyz]+=4.*stream_density(c,local,pair,i,j)*h;
   }
  }
 }
}
int CpuVnaForceStream_Run(double *****density){
 int mode=CpuVnaForceStream_Mode();if(!mode)return 0;if(!density||!F_VNA_flag)stream_fail("native VNA force/density unavailable");
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);size_t count=3*((size_t)atomnum+1);if(count>INT32_MAX)stream_fail("force collective count exceeds MPI int");
 StreamForce c={density,calloc(count,sizeof(double)),0,(List_YOUSO[35]+1)*(List_YOUSO[35]+1)*List_YOUSO[34]};double *sum=calloc(count,sizeof(double));if(!c.force||!sum)stream_fail("force workspace allocation failed");
 double start=MPI_Wtime();native_extras(&c);double dt=CpuVnaDeferred_Stream(DS_VNA,consume_pair,&c);(void)dt;
 MPI_Allreduce(c.force,sum,(int)count,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 for(int local=1;local<=Matomnum;local++){int a=M2G[local];for(int xyz=0;xyz<3;xyz++){Gxyz[a][41+xyz]=sum[3*(size_t)a+xyz];Gxyz[a][17+xyz]+=sum[3*(size_t)a+xyz];}}
 /* This fixed-geometry path has no remaining projector-value readers. Retain
  * the shared window and ownership metadata so native cleanup still knows
  * that its row views must not be passed to free(). */
 size_t advised=CpuVnaShared_DropRaw0();
 if(CpuVnaShared_OwnedBytes()){printf("CPU_VNA_FORCE_CACHE_RELEASE stage=after_stream_force4 rank=%d raw_payload_bytes=%zu advised_range_bytes=%zu shared_window_retained=1\n",rank,CpuVnaShared_OwnedBytes()/4,advised);fflush(stdout);}
 unsigned long long pairs=c.pairs,total=0;MPI_Reduce(&pairs,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,mpi_comm_level1);double elapsed=MPI_Wtime()-start,critical=0.;MPI_Reduce(&elapsed,&critical,1,MPI_DOUBLE,MPI_MAX,0,mpi_comm_level1);
 if(!rank){printf("CPU_VNA_FORCE_STREAM mode=%d pairs=%llu derivative_components_resident=private_pair gradient_publish=0 critical_seconds=%.9f\n",mode,total,critical);fflush(stdout);}free(c.force);free(sum);return 1;
}
