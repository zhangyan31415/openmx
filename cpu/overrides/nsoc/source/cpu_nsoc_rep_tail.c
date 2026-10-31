#include "cpu_nsoc_rep_tail.h"
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
int CpuRepTailAdmit(MPI_Comm comm,int nao,int first,int last,const double*c,const double*f,const double*energy,int count,const int*ids,double amplification,double tolerance,CpuRepTailEvidence*out){
 int np;MPI_Comm_size(comm,&np);int64_t span=last>=first?(int64_t)last-first+1:0;int nk=span<=INT_MAX?(int)span:0,bad=span>INT_MAX||nao<1||first<1||last<first-1||count<0||count==INT_MAX||(count&&!ids)||(nk&&(!c||!f||!energy))||!isfinite(amplification)||amplification<1.||!isfinite(tolerance)||tolerance<=0.||tolerance>1e-13;
 if(!bad&&nk&&(size_t)nao>SIZE_MAX/sizeof(double)/(size_t)nk)bad=1;
 if(!bad)for(int q=0;q<count;q++)if(ids[q]<1||(q&&ids[q]<=ids[q-1]))bad=1;
 uint64_t hash=UINT64_C(1469598103934665603),bits;
 hash=(hash^(uint32_t)nao)*UINT64_C(1099511628211);hash=(hash^(uint32_t)count)*UINT64_C(1099511628211);
 memcpy(&bits,&amplification,sizeof(bits));hash=(hash^bits)*UINT64_C(1099511628211);memcpy(&bits,&tolerance,sizeof(bits));hash=(hash^bits)*UINT64_C(1099511628211);
 if(!bad)for(int q=0;q<count;q++)hash=(hash^(uint32_t)ids[q])*UINT64_C(1099511628211);
 unsigned long long protocol[3]={hash,UINT64_MAX-hash,bad?0:1},agreed[3];MPI_Allreduce(protocol,agreed,3,MPI_UNSIGNED_LONG_LONG,MPI_MIN,comm);
 if(!agreed[2]||agreed[0]!=UINT64_MAX-agreed[1]){if(out)*out=(CpuRepTailEvidence){0,count,0,1,INFINITY,INFINITY,amplification,tolerance};return 0;}
 int*coverage=calloc((size_t)count+1,sizeof(int)),allocation_bad=!coverage,any_allocation_bad;MPI_Allreduce(&allocation_bad,&any_allocation_bad,1,MPI_INT,MPI_MAX,comm);
 if(any_allocation_bad){free(coverage);if(out)*out=(CpuRepTailEvidence){0,count,0,1,INFINITY,INFINITY,amplification,tolerance};return 0;}
 long double d=0.,e=0.;int covered=0;
 if(!bad)for(int q=0;q<count;q++){
  int state=ids[q];if(state<first||state>last)continue;int local=state-first;covered++;coverage[q]=1;
  if(!isfinite(f[local])||!isfinite(energy[local])){bad=1;continue;}
  long double norm=0.;for(int row=0;row<nao;row++){double x=c[(size_t)row*nk+local];if(!isfinite(x)){bad=1;break;}norm+=(long double)x*(long double)x;}
  long double weight=fabsl((long double)f[local]);d+=weight*norm;e+=weight*fabsl((long double)energy[local])*norm;
  if(!isfinite(norm)||!isfinite(d)||!isfinite(e))bad=1;
 }
 double local[2]={0.,0.},global[2];
 if(!bad){local[0]=nextafter((double)d,INFINITY);local[1]=nextafter((double)e,INFINITY);if(!isfinite(local[0])||!isfinite(local[1]))bad=1;}
 /* No NaN enters MPI_SUM. A separate integer bad flag rejects globally. */
 if(bad)local[0]=local[1]=0.;
 coverage[count]=bad;MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,comm);MPI_Allreduce(MPI_IN_PLACE,coverage,count+1,MPI_INT,MPI_SUM,comm);
 bad=coverage[count]!=0||!isfinite(global[0])||!isfinite(global[1]);covered=0;for(int q=0;q<count;q++){if(coverage[q]!=1)bad=1;covered+=coverage[q];}free(coverage);
 /* Positive sums, products and conversion are conservatively inflated to
  * cover FP64 roundoff, including arbitrary MPI reduction association. Local
  * norms use long double; this FP64 gamma envelope is deliberately larger. */
 double ku=4.*DBL_EPSILON*((double)nao+(double)np+(double)count+16.);
 if(!isfinite(ku)||ku>=.125)bad=1;
 double db=INFINITY,eb=INFINITY;
 if(!bad){double inflate=1./(1.-ku);db=nextafter(global[0]*amplification*inflate,INFINITY);eb=nextafter(global[1]*amplification*inflate,INFINITY);if(!isfinite(db)||!isfinite(eb))bad=1;}
 int admitted=!bad&&db<=tolerance&&eb<=tolerance;
 if(out)*out=(CpuRepTailEvidence){admitted,count,covered,bad,db,eb,amplification,tolerance};return admitted;
}

double CpuRepTailAmplification(int operations,int species,const int*norb,const double*const*rotations,double*row_l1){
 if(operations<1||species<1||!norb||!rotations)return NAN;
 int width_max=0;for(int s=0;s<species;s++){if(norb[s]<1||norb[s]>128)return NAN;if(norb[s]>width_max)width_max=norb[s];}
 long double maximum=0.;for(int g=0;g<operations;g++){
  if(!rotations[g])return NAN;size_t offset=0;
  for(int s=0;s<species;s++){int n=norb[s];for(int i=0;i<n;i++){long double sum=0.;for(int j=0;j<n;j++){double x=rotations[g][offset+(size_t)i*n+j];if(!isfinite(x))return NAN;sum+=fabsl((long double)x);}if(sum>maximum)maximum=sum;}offset+=(size_t)n*n;}
 }
 double l=nextafter((double)maximum,INFINITY),inflate=1.+8.*DBL_EPSILON*((double)width_max+2.);
 double bound=nextafter(1.+l*l*inflate,INFINITY);if(!isfinite(l)||!isfinite(bound)||maximum<=0.)return NAN;if(row_l1)*row_l1=l;return bound;
}
