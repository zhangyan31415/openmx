#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "state_coefficients.h"
#include "cpu_nsoc_owner.h"
#include "cpu_nsoc_rep_tail.h"
#ifndef CPU_NO_STAR_CACHE_BYTES
#define CPU_NO_STAR_CACHE_BYTES (64u*1024u*1024u)
#endif
static int own(int g,int b,int p){return(g/b)%p;}
static int local(int g,int b,int p){return(g/(b*p))*b+g%b;}
static double sign4(int a,int b){return((a&b)==1||(a&b)==2)?-1.:1.;}
static double entry(int row,int column){return .03125*(1+row+3*column);}
typedef struct {CpuStateCoefficients*c;int calls;} Tracked;
static int rows(void*context,int first,int count,double*out,int ld){Tracked*t=context;t->calls++;return CpuStateRows(t->c,first,count,out,ld);}
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
 int n=48,pc=(int)sqrt((double)np);while(np%pc)pc--;int pr=np/pc,physical=(rank+1)%np,row=physical/pc,col=physical%pc;
 int64_t cp[49];int ri[192],inv[48],ends[4]={12,24,36,48};double cv[192];
 for(int q=0;q<n;q++){cp[q]=4*q;inv[q]=q*5%n;for(int a=0;a<4;a++){ri[4*q+a]=4*(q%12)+a;cv[4*q+a]=.5*sign4(a,q/12);}}cp[n]=192;
 CpuNativeSector sectors[4];size_t source_used=0;
 for(int s=0;s<4;s++){int nr=0,nc=0;for(int g=0;g<12;g++){nr+=own(g,3,pr)==row;nc+=own(g,3,pc)==col;}sectors[s]=(CpuNativeSector){12,s*12,nr,nc,3,source_used};source_used+=(size_t)nr*nc;}
 double*source=malloc((source_used?source_used:1)*sizeof(double));
 for(int s=0;s<4;s++)for(int q=0;q<12;q++)if(own(q,3,pc)==col)for(int i=0;i<12;i++)if(own(i,3,pr)==row)
  source[sectors[s].offset+(size_t)local(q,3,pc)*sectors[s].nr+local(i,3,pr)]=entry(i+s*12,q+s*12);
 int neighbors[6]={0,1,2,3,1,0};CpuNOBlock blocks[72];for(int a=0;a<12;a++)for(int h=0;h<6;h++)blocks[6*a+h]=(CpuNOBlock){4*a,4,4*((a+neighbors[h])%12),4,a%np};
 CpuNOPlan*plan=CpuNOCreate(MPI_COMM_WORLD,n,72,blocks,1024);if(!plan)MPI_Abort(MPI_COMM_WORLD,1);
 CpuStateCoefficients state={0};int bad=CpuStateInit(&state,n,cp,ri,cv);size_t entries=CpuNOEntries(plan);
 for(int trial=0;trial<3;trial++){
  int wanted=trial==0?n:trial==1?26:3;
  int first=(int)((double)wanted/np*rank)+1,last=(int)((double)wanted/np*(rank+1));if(rank==np-1)last=wanted;
  if(wanted<np){first=rank<wanted?rank+1:1;last=rank<wanted?rank+1:0;}
  int nk=last>=first?last-first+1:0;size_t capacity=(size_t)n*nk;
  double*q=malloc((capacity+8)*sizeof(double));if(CpuStateBind(&state,4,ends,inv,wanted,first,last,q,capacity))bad=1;
  for(size_t k=0;k<capacity+8;k++)q[k]=1234567.;
  size_t temporary=0;if(CpuStatePack(MPI_COMM_WORLD,pr,pc,row,col,sectors,source,inv,wanted,&state,7,&temporary))bad=1;
  for(size_t k=state.used;k<capacity+8;k++)if(q[k]!=1234567.)bad=1;
  double*reference=malloc((capacity?capacity:1)*sizeof(double));
  for(int original=0;original<n;original++)if(inv[original]>=first-1&&inv[original]<last)for(int i=0;i<n;i++)
   reference[(size_t)i*nk+inv[original]-(first-1)]=.5*sign4(i%4,original/12)*entry((original/12)*12+i/4,original);
  state.reference=reference;
  double*ao=malloc((capacity+8)*sizeof(double));for(size_t k=0;k<capacity+8;k++)ao[k]=98765.;
  if(CpuStateRows(&state,0,n,ao,nk))bad=1;
  for(size_t k=0;k<capacity;k++)if(ao[k]!=reference[k])bad=1;
  for(size_t k=capacity;k<capacity+8;k++)if(ao[k]!=98765.)bad=1;
  double*f=malloc((size_t)(nk?nk:1)*sizeof(double)),*energy=malloc((size_t)(nk?nk:1)*sizeof(double));
  for(int k=0;k<nk;k++){f[k]=.125*(1+(k%3));energy[k]=.0625*(first+k)-1.75;}
  if(CpuNOCompute(plan,nk,reference,f,energy))bad=1;
  double*saved=malloc((entries?2*entries:1)*sizeof(double));memcpy(saved,CpuNODensity(plan),entries*sizeof(double));memcpy(saved+entries,CpuNOEnergyDensity(plan),entries*sizeof(double));
  Tracked tracked={&state,0};CpuNORowProvider provider={&tracked,rows};
  if(CpuNOComputeProvided(plan,nk,&provider,f,energy))bad=1;
  for(size_t k=0;k<entries;k++)if(saved[k]!=CpuNODensity(plan)[k]||saved[entries+k]!=CpuNOEnergyDensity(plan)[k])bad=1;
  int expected_calls=CPU_NO_STAR_CACHE_BYTES<(size_t)4*4*nk*sizeof(double)?60:48;
  if(nk&&tracked.calls!=expected_calls)bad=1;
  int ids[2]={1,wanted},cuts=wanted>1?2:1;long double*norms=malloc((size_t)(nk?nk:1)*sizeof(long double));
  if(CpuStateCutNorms(&state,cuts,ids,norms))bad=1;
  CpuRepTailEvidence a,b;int old=CpuRepTailAdmit(MPI_COMM_WORLD,n,first,last,reference,f,energy,cuts,ids,2.,1e-13,&a);
  int provided=CpuRepTailAdmitNorms(MPI_COMM_WORLD,n,first,last,norms,f,energy,cuts,ids,2.,1e-13,&b);
  if(old!=provided||a.density_bound!=b.density_bound||a.energy_density_bound!=b.energy_density_bound||a.invalid!=b.invalid||a.covered_states!=b.covered_states)bad=1;
  for(int k=0;k<nk;k++)f[k]*=1e-30;
  old=CpuRepTailAdmit(MPI_COMM_WORLD,n,first,last,reference,f,energy,cuts,ids,2.,1e-13,&a);
  provided=CpuRepTailAdmitNorms(MPI_COMM_WORLD,n,first,last,norms,f,energy,cuts,ids,2.,1e-13,&b);
  if(!old||!provided||a.density_bound!=b.density_bound||a.energy_density_bound!=b.energy_density_bound)bad=1;
  if(first<=ids[0]&&ids[0]<=last)norms[ids[0]-first]=NAN;
  if(CpuRepTailAdmitNorms(MPI_COMM_WORLD,n,first,last,norms,f,energy,cuts,ids,2.,1e-13,&b)||!b.invalid)bad=1;
  if(state.audit_error!=0.)bad=1;
  /* Retirement keeps sparse row metadata for the next Bind but must reject
   * an already-issued provider after its coefficient pages become dead. */
  state.data=NULL;if(!CpuStateRows(&state,0,1,ao,nk))bad=1;
  free(saved);free(norms);free(f);free(energy);free(ao);free(q);
 }
 CpuStateDestroy(&state);CpuNODestroy(plan);free(source);
 int any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
 if(!rank)printf("state_owner_density ranks=%d Q_full_partial_empty=PASS AO_coefficients=PASS star_cache=PASS signed_DM_EDM=PASS tail_bound=%s\n",np,any?"FAIL":"PASS");
 MPI_Finalize();return any?1:0;
}
