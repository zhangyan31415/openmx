#include "soc_recv_capacity.h"
#include "cpu_memory_phase.h"
/* Private Gamma-point C3 spinor generalized eigensolver. Original H/S and all
 * SOC components retained; modes2/3 explicitly project onto C3 symmetry.
 * Dense storage reuses the caller's three existing spinor buffers. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include "mpi.h"
#include "openmx_common.h"
#include "cpu_c3.h"
#include "cpu_soc_band_exchange.h"
extern int numroc_(int*,int*,int*,int*,int*);
extern void descinit_(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void pzpotrf_(char*,int*,dcomplex*,int*,int*,int*,int*);
extern void pzhegst_(int*,char*,int*,dcomplex*,int*,int*,int*,dcomplex*,int*,int*,int*,double*,int*);
extern void pztranc_(int*,int*,dcomplex*,dcomplex*,int*,int*,int*,dcomplex*,dcomplex*,int*,int*,int*);
extern void pztrsm_(char*,char*,char*,char*,int*,int*,dcomplex*,dcomplex*,int*,int*,int*,dcomplex*,int*,int*,int*);
extern void elpa_solve_evp_complex_2stage_double_impl_(int*,int*,dcomplex*,int*,double*,dcomplex*,int*,int*,int*,int*,int*,int*);
extern int CpuC3_GeometryToleranceAng(const char*,double*);
extern double CpuC3_PeriodicDistanceBohr(const double*,const double*);
extern int flag_export_DM;

typedef struct {int n,ofs,nr,nc,nb,desc[9];size_t soff;} Sector;
typedef struct {int band,row;dcomplex z;} Coeff;
static int active_config=0,force_full=0,full_steps=0,last_step=0;
static struct {int loaded,mode,N,ns,natoms,nsym,d3,sz[4],nnz,*perm,*c2perm,*ri,*rp,*ci;int64_t *cp;dcomplex *cv,*rv;double shift[3],c2rot[9],c2shift[3];Sector sec[4];int nsec;size_t sfsize;} B;
static dcomplex zz(double r,double i){dcomplex z;z.r=r;z.i=i;return z;}
static dcomplex add(dcomplex a,dcomplex b){return zz(a.r+b.r,a.i+b.i);}
static dcomplex mul(dcomplex a,dcomplex b){return zz(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);}
static dcomplex cj(dcomplex a){return zz(a.r,-a.i);}
static void die(const char *s){int r;MPI_Comm_rank(mpi_comm_level1,&r);fprintf(stderr,"CPU_C3 ERROR rank=%d %s\n",r,s);fflush(stderr);MPI_Abort(mpi_comm_level1,95);abort();}
static void *mem(size_t n,size_t s){if(s && n>SIZE_MAX/s)die("allocation overflow");void *p=calloc(n?n:1,s);if(!p)die("allocation failed");return p;}
static void rd(void *p,size_t s,size_t n,FILE*f){if(fread(p,s,n,f)!=n)die("short basis file");}
int CpuC3_Mode(void){const char*p=getenv("OPENMX_CPU_C3");int m=p?atoi(p):0;if(m<0||m>6)die("unsupported mode");return m;}
int CpuC3_WarmupActive(void){return active_config==6&&B.loaded&&B.mode==3;}
void CpuC3_RequestFull(void){force_full=1;}
int CpuC3_ConvergenceBlocked(void){return active_config==6&&(!B.loaded||B.mode==3||full_steps<2);}
int CpuC3_FullSteps(void){return full_steps;}
int CpuC3_MixingAge(void){
 const char*reset=getenv("OPENMX_CPU_C3_RESET_DIIS");
 if(!reset||!atoi(reset)||active_config!=6||B.mode!=4||Mixing_switch!=1||last_step-full_steps+1<Pulay_SCF)return 0;
 return full_steps;
}
static void bounds(int n,int rank,int size,int *lo,int *hi){*lo=(int)((double)n/size*rank);*hi=rank==size-1?n:(int)((double)n/size*(rank+1));}
static int owner(int n,int p,int g){return (int)(((int64_t)(g+1)*p-1)/n);}
static int sec_of(int col){return col<B.sz[0]?0:col<B.sz[0]+B.sz[1]?1:2;}
static void geometry(void){
 if(B.mode==1||B.mode==4)return;
 double geometry_tol=5e-5,cell_tol=5e-5; /* experimental v18: allow projected-all smoke on measured symmetric cells */
 const char *geometry_input=getenv("OPENMX_CPU_C3_GEOMETRY_TOL_ANG");
 if(geometry_input){
  double tolerance_ang;
  if(!CpuC3_GeometryToleranceAng(geometry_input,&tolerance_ang))die("geometry tolerance must be finite, positive and <=0.02 Ang");
  geometry_tol=tolerance_ang/BohrR;cell_tol=geometry_tol;
 }
 double mx=0,c=-.5,s=sqrt(3.)/2,cell[9];
 for(int i=0;i<3;i++)for(int j=0;j<3;j++)cell[3*i+j]=tv[i+1][j+1];
 double flat=fabs(tv[1][3])+fabs(tv[2][3])+fabs(tv[3][1])+fabs(tv[3][2]);
 if(flat>(geometry_input?4*cell_tol:1e-8))die("requires a slab cell within the declared geometry tolerance");
 double drcell[3]={c*tv[1][1]-s*tv[1][2]-tv[2][1],s*tv[1][1]+c*tv[1][2]-tv[2][2],tv[1][3]-tv[2][3]};
 if(sqrt(drcell[0]*drcell[0]+drcell[1]*drcell[1]+drcell[2]*drcell[2])>cell_tol)die("cell is not C3 hexagonal within declared tolerance");
 for(int a=0;a<B.natoms;a++){
  int b=B.perm[a];if(b<0||b>=B.natoms||WhatSpecies[a+1]!=WhatSpecies[b+1])die("invalid atom permutation");
  double target[3]={c*Gxyz[a+1][1]-s*Gxyz[a+1][2],s*Gxyz[a+1][1]+c*Gxyz[a+1][2],Gxyz[a+1][3]},dr[3];
  for(int i=0;i<3;i++){for(int j=0;j<3;j++)target[i]+=B.shift[j]*tv[j+1][i+1];dr[i]=Gxyz[b+1][i+1]-target[i];}
  double er=CpuC3_PeriodicDistanceBohr(dr,cell);if(!isfinite(er))die("invalid full-cell periodic geometry metric");if(er>mx)mx=er;
 }
 if(mx>geometry_tol)die("geometry breaks C3; use separately validated symmetric input");
 if(B.d3){
  mx=0;
  for(int a=0;a<B.natoms;a++){
   int b=B.c2perm[a];if(b<0||b>=B.natoms||WhatSpecies[a+1]!=WhatSpecies[b+1])die("invalid C2 atom permutation");
   double dr[3];for(int i=0;i<3;i++){double x=0;for(int j=0;j<3;j++)x+=B.c2rot[3*i+j]*Gxyz[a+1][j+1]+B.c2shift[j]*tv[j+1][i+1];dr[i]=Gxyz[b+1][i+1]-x;}
   double er=CpuC3_PeriodicDistanceBohr(dr,cell);if(!isfinite(er))die("invalid C2 full-cell periodic geometry metric");if(er>mx)mx=er;
  }
  if(mx>geometry_tol)die("geometry breaks D3 C2 operation");
 }
}
/* Internal sparse gauges; raw physical operators are checked by MarkSolve. */
static struct {int tr,real,prepared,check,pick[4],*sector_by_state;} AU;
extern int CpuSOCNative_TRPassed(void);
extern int CpuSOCNative_TRRequested(void);
void CpuC3_AntiunitaryConfigure(int tr,int real){
 if((tr!=0&&tr!=1)||(real!=0&&real!=1)||real>tr)die("antiunitary requires0/1 and REAL requires TR");
 if(B.loaded)die("antiunitary configuration after basis load");
 if((tr||real)&&!CpuC3_CompactEnabled())die("antiunitary requires compact coefficients");
 AU.tr=tr;AU.real=real;AU.prepared=0;
}
int CpuC3_AntiunitaryReal(void){return AU.real;}
int CpuC3_AntiunitaryTR(void){return AU.tr;}
static void au_record_state_sectors(const int*inverse,int one_based){
 if(!AU.check)return;free(AU.sector_by_state);AU.sector_by_state=mem(B.ns,sizeof(int));
 for(int k=0;k<B.nsec;k++){
  AU.pick[k]=inverse[B.sec[k].ofs]+(one_based?0:1);
  for(int b=0;b<B.sec[k].n;b++){int sorted=inverse[B.sec[k].ofs+b]-(one_based?1:0);if(sorted<0||sorted>=B.ns)die("bad verifier sorted sector map");AU.sector_by_state[sorted]=k;}
 }
}
static int au_validation_picks(int maxn,int*pick){
 int count=0,initial[3]={1,(maxn+1)/2,maxn};
 if(!AU.sector_by_state||maxn<1)die("missing full-sector verifier map");
 for(int q=0;q<3+B.nsec;q++){
  int id=q<3?initial[q]:AU.pick[q-3];if(id<1||id>maxn)die("CHECK sector sample outside returned state range");
  int found=0;for(int j=0;j<count;j++)if(pick[j]==id)found=1;if(!found)pick[count++]=id;
 }return count;
}
static void au_monomial(int source,int target,int sign,int*perm,dcomplex*phase){
 int n=B.sz[source],so=0,to=0;for(int k=0;k<source;k++)so+=B.sz[k];for(int k=0;k<target;k++)to+=B.sz[k];
 int*seen=mem(n,sizeof(int));
 for(int j=0;j<n;j++){
  int ids[256],count=0;dcomplex vals[256];
  for(int64_t a=B.cp[so+j];a<B.cp[so+j+1];a++){
   int old=B.ri[a],row=old<B.N?old+B.N:old-B.N;dcomplex z=cj(B.cv[a]);if(old<B.N){z.r=-z.r;z.i=-z.i;}
   for(int b=B.rp[row];b<B.rp[row+1];b++){
    int i=B.ci[b]-to;if(i<0||i>=n)continue;int q=0;while(q<count&&ids[q]!=i)q++;
    if(q==count){if(count==256)die("antiunitary sparse overlap degree");ids[count]=i;vals[count++]=zz(0,0);}
    vals[q]=add(vals[q],mul(cj(B.rv[b]),z));
   }
  }
  int found=-1;dcomplex z=zz(0,0);
  for(int q=0;q<count;q++)if(hypot(vals[q].r,vals[q].i)>1e-12){if(found>=0)die("antiunitary map is not monomial");found=ids[q];z=vals[q];}
  if(found<0||fabs(hypot(z.r,z.i)-1.)>1e-11||seen[found]++)die("antiunitary monomial is not unitary");
  perm[j]=found;phase[j]=zz(sign*z.r,sign*z.i);
 }
 free(seen);
}
static void au_rebuild_rows(void){
 free(B.rp);free(B.ci);free(B.rv);B.rp=mem(B.ns+1,sizeof(int));B.ci=mem(B.nnz,sizeof(int));B.rv=mem(B.nnz,sizeof(dcomplex));
 for(int a=0;a<B.nnz;a++)B.rp[B.ri[a]+1]++;
 for(int i=0;i<B.ns;i++)B.rp[i+1]+=B.rp[i];int*pos=mem(B.ns,sizeof(int));memcpy(pos,B.rp,(size_t)B.ns*sizeof(int));
 for(int j=0;j<B.ns;j++)for(int64_t a=B.cp[j];a<B.cp[j+1];a++){int q=pos[B.ri[a]]++;B.ci[q]=j;B.rv[q]=B.cv[a];}free(pos);
}
static double au_column_t_error(int source,int target,int sign){
 int capacity=(int)(B.cp[source+1]-B.cp[source]+B.cp[target+1]-B.cp[target]);
 int*rows=mem(capacity,sizeof(int));dcomplex*vals=mem(capacity,sizeof(dcomplex));int used=0;
 for(int part=0;part<2;part++){
  int col=part?target:source;
  for(int64_t a=B.cp[col];a<B.cp[col+1];a++){
   int row=B.ri[a];dcomplex z=B.cv[a];
   if(!part){z=cj(z);if(row<B.N){row+=B.N;z.r=-z.r;z.i=-z.i;}else row-=B.N;}
   else{z.r*=-sign;z.i*=-sign;}
   int q=0;while(q<used&&rows[q]!=row)q++;if(q==used){rows[used]=row;vals[used++]=zz(0,0);}vals[q]=add(vals[q],z);
  }
 }
 double norm=0;for(int i=0;i<used;i++)norm+=vals[i].r*vals[i].r+vals[i].i*vals[i].i;free(rows);free(vals);return sqrt(norm);
}
static void antiunitary_prepare(void){
 if(!AU.tr)return;if(AU.prepared)die("repeated antiunitary basis transform");
 if(!B.d3||B.mode!=3||B.nsym!=4)die("antiunitary supports projected D3 only");
 int m=B.sz[0],r=B.sz[2];int*bp=mem(m,sizeof(int)),*tp=mem(r,sizeof(int));dcomplex*ba=mem(m,sizeof(dcomplex)),*ta=mem(r,sizeof(dcomplex));
 au_monomial(0,1,-1,bp,ba);au_monomial(2,3,1,tp,ta);
 for(int j=0;j<m;j++)if(bp[bp[j]]!=j||hypot(ba[j].r-ba[bp[j]].r,ba[j].i-ba[bp[j]].i)>1e-11)die("C2T map is not a symmetric involution");
 int maxcol=0;for(int j=0;j<B.ns;j++){int n=(int)(B.cp[j+1]-B.cp[j]);if(n>maxcol)maxcol=n;}
 if(B.nnz>INT_MAX/2||maxcol>INT_MAX/2)die("antiunitary sparse gauge overflow");
 int64_t*cp=mem(B.ns+1,sizeof(int64_t));int*ri=mem((size_t)B.nnz*2,sizeof(int));dcomplex*cv=mem((size_t)B.nnz*2,sizeof(dcomplex));
 int*rows=mem((size_t)maxcol*2,sizeof(int));dcomplex*vals=mem((size_t)maxcol*2,sizeof(dcomplex));int nnz=0,paired=0;
 for(int col=0;col<B.ns;col++){
  int cols[2]={col,0},nc=1;dcomplex weights[2]={zz(1,0),zz(0,0)};
  if(AU.real&&col<2*m){
   int sec=col/m,j=col%m,i=bp[j],low=j<i?j:i,high=j>i?j:i;
   double angle=.5*atan2(ba[low].i,ba[low].r);dcomplex phase=zz(cos(angle),sin(angle));
   if(i==j){weights[0]=phase;}
   else{nc=2;cols[0]=sec*m+low;cols[1]=sec*m+high;
    if(j==low){weights[0]=weights[1]=zz(phase.r/sqrt(2.),phase.i/sqrt(2.));if(!sec)paired++;}
    else{weights[0]=zz(-phase.i/sqrt(2.),phase.r/sqrt(2.));weights[1]=zz(-weights[0].r,-weights[0].i);}
   }
  }else if(col>=2*m+r){int j=col-(2*m+r);cols[0]=2*m+r+tp[j];weights[0]=ta[j];}
  int used=0;cp[col]=nnz;
  for(int c=0;c<nc;c++)for(int64_t a=B.cp[cols[c]];a<B.cp[cols[c]+1];a++){
   int row=B.ri[a],q=0;while(q<used&&rows[q]!=row)q++;
   if(q==used){rows[used]=row;vals[used++]=zz(0,0);}vals[q]=add(vals[q],mul(B.cv[a],weights[c]));
  }
  for(int q=0;q<used;q++){ri[nnz]=rows[q];cv[nnz++]=vals[q];}
 }
 cp[B.ns]=nnz;free(rows);free(vals);free(B.cp);free(B.ri);free(B.cv);B.cp=cp;B.ri=ri;B.cv=cv;B.nnz=nnz;au_rebuild_rows();
 double leakage=0.;for(int j=0;j<r;j++){double err=au_column_t_error(2*m+j,2*m+r+j,1);if(err>leakage)leakage=err;}
 if(AU.real)for(int j=0;j<m;j++){double err=au_column_t_error(j,m+j,-1);if(err>leakage)leakage=err;}
 if(!isfinite(leakage)||leakage>1e-11)die("antiunitary sparse gauge has physical-sector leakage");
 free(bp);free(tp);free(ba);free(ta);AU.prepared=1;
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);if(!rank){printf("CPU_SOC_ANTIUNITARY_BASIS tr=%d real=%d real_gauge_pairs=%d sparse_nnz=%d physical_map_error=%.12e\n",AU.tr,AU.real,paired,B.nnz,leakage);fflush(stdout);}
}

static void c3_layout(int mode){
 int p,r,zero=0,info;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);B.mode=mode;
 B.nsec=mode==3?B.nsym:1;B.sfsize=0;int ofs=0;
 for(int k=0;k<B.nsec;k++){
  Sector*t=B.sec+k;t->n=mode==3?B.sz[k]:B.ns;t->ofs=mode==3?ofs:0;ofs+=t->n;t->nb=nblk2;while(t->nb>1&&(t->nb*np_rows2>t->n||t->nb*np_cols2>t->n))t->nb/=2;
  t->nr=numroc_(&t->n,&t->nb,&my_prow2,&zero,&np_rows2);t->nc=numroc_(&t->n,&t->nb,&my_pcol2,&zero,&np_cols2);t->soff=B.sfsize;if(!(B.d3&&mode==3&&(k==1||(AU.tr&&k==3))))B.sfsize+=(size_t)t->nr*t->nc;
  int ld=t->nr?t->nr:1;descinit_(t->desc,&t->n,&t->n,&t->nb,&t->nb,&zero,&zero,&ictxt1_2,&ld,&info);if(info)die("descriptor initialization");
 }
 if(B.sfsize>(size_t)na_rows2*na_cols2)die("existing S buffer too small");geometry();
}
static void c3_init(int mode){
 int p,r,zero=0,info;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);
 if(B.loaded){if(B.mode!=mode){if(B.mode!=3||mode!=4)die("unsupported mode transition");c3_layout(mode);}geometry();return;}
 B.mode=mode;const char*path=getenv("OPENMX_CPU_C3_BASIS");if(!path)die("missing OPENMX_CPU_C3_BASIS");FILE*f=fopen(path,"rb");if(!f)die("cannot read basis");
 int32_t hd[9]={0};int64_t nz;rd(hd,4,8,f);if(hd[1]==2)rd(hd+8,4,1,f);rd(&nz,8,1,f);
 if(hd[0]!=0x43334253||!((hd[1]==1&&hd[4]==3)||(hd[1]==2&&hd[4]==4))||hd[2]!=atomnum||hd[3]<1||hd[3]>INT_MAX/2||nz<1||nz>32LL*hd[3]*2)die("invalid basis header");
 B.natoms=hd[2];B.N=hd[3];B.ns=2*B.N;B.nsym=hd[4];B.d3=B.nsym==4;B.nnz=(int)nz;int64_t sumsz=0;for(int k=0;k<B.nsym;k++){B.sz[k]=hd[5+k];if(B.sz[k]<1)die("empty symmetry sector");sumsz+=B.sz[k];}
 if(sumsz!=B.ns||(B.d3&&(B.sz[0]!=B.sz[1]||B.sz[2]!=B.sz[3])))die("invalid sector sizes");
 const char*group=getenv("OPENMX_CPU_D3");if((group&&atoi(group)!=0)!=B.d3)die("D3 selector and basis header disagree");
 if(B.d3&&(mode==1||mode==2))die("D3 supports adaptive/projected-block/raw-full modes only");
 int count=0;for(int a=1;a<=atomnum;a++)count+=Spe_Total_CNO[WhatSpecies[a]];if(count!=B.N||SpinP_switch!=3||SO_switch!=1||Cnt_switch)die("unsupported basis or spin mode");
 B.perm=mem(B.natoms,sizeof(int));rd(B.perm,4,B.natoms,f);rd(B.shift,8,3,f);
 if(B.d3){B.c2perm=mem(B.natoms,sizeof(int));rd(B.c2perm,4,B.natoms,f);rd(B.c2rot,8,9,f);rd(B.c2shift,8,3,f);for(int i=0;i<9;i++)if(!isfinite(B.c2rot[i]))die("nonfinite C2 rotation");for(int i=0;i<3;i++){if(!isfinite(B.c2shift[i]))die("nonfinite C2 shift");for(int j=0;j<3;j++){double v=0;for(int k=0;k<3;k++)v+=B.c2rot[3*k+i]*B.c2rot[3*k+j];if(fabs(v-(i==j))>1e-12)die("nonorthogonal C2 rotation");}}}
 B.cp=mem(B.ns+1,8);B.ri=mem(nz,4);B.cv=mem(nz,sizeof(dcomplex));rd(B.cp,8,B.ns+1,f);rd(B.ri,4,nz,f);rd(B.cv,sizeof(dcomplex),nz,f);if(fgetc(f)!=EOF)die("trailing basis data");fclose(f);
 if(B.cp[0]!=0||B.cp[B.ns]!=nz)die("invalid CSC pointer");
 if(mode==1){free(B.cp);free(B.ri);free(B.cv);B.nnz=B.ns;B.cp=mem(B.ns+1,8);B.ri=mem(B.ns,4);B.cv=mem(B.ns,sizeof(dcomplex));for(int i=0;i<B.ns;i++){B.cp[i]=i;B.ri[i]=i;B.cv[i]=zz(1,0);}B.cp[B.ns]=B.ns;}
 B.rp=mem(B.ns+1,4);B.ci=mem(B.nnz,4);B.rv=mem(B.nnz,sizeof(dcomplex));
 for(int j=0;j<B.ns;j++){
  if(B.cp[j]>B.cp[j+1]||B.cp[j+1]>B.nnz)die("invalid CSC range");double nrm=0;
  for(int64_t k=B.cp[j];k<B.cp[j+1];k++){int i=B.ri[k];dcomplex z=B.cv[k];if(i<0||i>=B.ns||!isfinite(z.r)||!isfinite(z.i))die("invalid basis value");B.rp[i+1]++;nrm+=z.r*z.r+z.i*z.i;}
  if(fabs(nrm-1)>1e-12)die("basis column not normalized");
 }
 for(int i=0;i<B.ns;i++)B.rp[i+1]+=B.rp[i];int*pos=mem(B.ns,4);memcpy(pos,B.rp,B.ns*4);
 for(int j=0;j<B.ns;j++)for(int64_t k=B.cp[j];k<B.cp[j+1];k++){int i=B.ri[k],t=pos[i]++;B.ci[t]=j;B.rv[t]=B.cv[k];}free(pos);
 if(AU.tr)antiunitary_prepare();
 /* Sparse Gram check includes cross-column overlaps, not just norms. */
 for(int j=0;j<B.ns;j++){
  int ids[256],cnt=0;dcomplex val[256];
  for(int64_t a=B.cp[j];a<B.cp[j+1];a++){
   int row=B.ri[a];for(int k=B.rp[row];k<B.rp[row+1];k++){
    int col=B.ci[k],q=0;while(q<cnt&&ids[q]!=col)q++;
    if(q==cnt){if(cnt==256)die("unsupported basis overlap degree");ids[cnt]=col;val[cnt++]=zz(0,0);}
    val[q]=add(val[q],mul(cj(B.cv[a]),B.rv[k]));
   }
  }
  for(int q=0;q<cnt;q++){if(ids[q]==j)val[q].r-=1;if(hypot(val[q].r,val[q].i)>1e-12)die("basis is not unitary");}
 }
 c3_layout(mode);B.loaded=1;
 if(!r){printf("CPU_C3 basis mode=%d ns=%d sectors=%d,%d,%d,%d d3=%d nnz=%d extra_dense_arrays=0\n",mode,B.ns,B.sz[0],B.sz[1],B.sz[2],B.sz[3],B.d3,B.nnz);fflush(stdout);}
}
static void project_entry(double*buf,int fields,int tile0,int tile1,Sector*t,int i,int j,dcomplex h,double s){
 int paired=B.d3&&B.mode==3&&t->ofs==0,m=B.sz[0];
 for(int b=B.rp[j];b<B.rp[j+1];b++){
  int gb=B.ci[b],tagb=paired?gb/m:0;if(paired&&tagb>1)continue;
  int cb=paired?gb-tagb*m:gb-t->ofs;if(cb<tile0||cb>=tile1)continue;
  for(int a=B.rp[i];a<B.rp[i+1];a++){
   int ga=B.ci[a];if(paired&&ga/m!=tagb)continue; /* Keep AA/BB, never fold AB/BA. */
   int ca=paired?ga-tagb*m:ga-t->ofs;if(ca<0||ca>=t->n)continue;
   if(B.mode==2&&sec_of(B.ci[a])!=sec_of(B.ci[b]))continue;
   dcomplex z=mul(cj(B.rv[a]),B.rv[b]);if(paired){z.r*=.5;z.i*=.5;}dcomplex v=mul(z,h);size_t off=((size_t)(cb-tile0)*t->n+ca)*fields;buf[off]+=v.r;buf[off+1]+=v.i;if(fields==4){buf[off+2]+=z.r*s;buf[off+3]+=z.i*s;}
  }
 }
}
static int soc_recv_mode(void){
 static int cached=-1;if(cached>=0)return cached;const char*text=getenv("OPENMX_CPU_SOC_TILE_RECV");int mode=0,valid=1;if(text){if(!strcmp(text,"0"))mode=0;else if(!strcmp(text,"1"))mode=1;else valid=0;}int good,lo,hi;MPI_Allreduce(&valid,&good,1,MPI_INT,MPI_MIN,mpi_comm_level1);if(!good)die("SOC_TILE_RECV must be0 or1");MPI_Allreduce(&mode,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(&mode,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);if(lo!=hi)die("SOC_TILE_RECV differs across ranks");cached=mode;return cached;
}
static size_t soc_recv_maximum(Sector*t,int tile,int fields,int rank,int world){
 size_t maximum;if(SocRecvMaximum(t->n,t->nb,np_rows2,np_cols2,my_prow2,my_pcol2,rank,world,tile,fields,&maximum))die("SOC recv maximum count/layout invalid");size_t rows;if(SocRecvPrefix(t->n,t->nb,my_prow2,np_rows2,&rows)||rows!=(size_t)t->nr)die("SOC recv actual descriptor rows inconsistent");return maximum;
}
static void soc_recv_alloc_record(Sector*t,int rank,int fields,int compact,size_t required,size_t allocated,size_t legacy){
 printf("CPU_SOC_RECV_ALLOC rank=%d sector=%d n=%d fields=%d compact=%d required_max_elements=%zu allocated_elements=%zu legacy_elements=%zu allocated_bytes=%zu legacy_bytes=%zu source=actual_mem_argument dummy_zero_owner=%d\n",rank,(int)(t-B.sec),t->n,fields,compact,required,allocated,legacy,allocated*sizeof(double),legacy*sizeof(double),required==0&&allocated==1);fflush(stdout);
}
static void soc_recv_used_record(Sector*t,int rank,int fields,int compact,size_t required,size_t allocated,size_t seen){
 if(required!=seen)die("SOC recv planned maximum disagrees with actual producer");printf("CPU_SOC_RECV_USED rank=%d sector=%d fields=%d compact=%d planned_max_elements=%zu observed_max_elements=%zu allocated_elements=%zu all_counts_bounded=1 all_reads_bounded=1\n",rank,(int)(t-B.sec),fields,compact,required,seen,allocated);fflush(stdout);
}
static void assemble(Sector*t,int first,double*****nh,double*****im,double****olp,dcomplex*H,dcomplex*S){
 int r,p;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);int*mp=mem(atomnum+1,4);int v=0;for(int a=1;a<=atomnum;a++){mp[a]=v;v+=Spe_Total_CNO[WhatSpecies[a]];}
 int fields=first?4:2,tile=128,compact_recv=soc_recv_mode();size_t cap=(size_t)t->n*tile*fields;
 size_t recv_max=soc_recv_maximum(t,tile,fields,r,p),recv_alloc=compact_recv?(recv_max?recv_max:1):cap,recv_seen=0;
 double*buf=mem(cap,8),*sum=r?NULL:mem(cap,8),*packed=r?NULL:mem(cap,8),*recv=mem(recv_alloc,8);int*ct=mem(p,4),*dp=mem(p,4),*pos=mem(p,4);
 soc_recv_alloc_record(t,r,fields,compact_recv,recv_max,recv_alloc,cap);CpuMemoryPhase("soc_complex_recv_allocated");
 memset(H,0,(size_t)t->nr*t->nc*sizeof(dcomplex));if(first)memset(S,0,(size_t)t->nr*t->nc*sizeof(dcomplex));
 for(int j0=0;j0<t->n;j0+=tile){
  int j1=j0+tile;if(j1>t->n)j1=t->n;size_t count=(size_t)t->n*(j1-j0)*fields;if(count>INT_MAX)die("tile count overflow");memset(buf,0,count*8);
  for(int a=1;a<=Matomnum;a++){
   int ga=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga]];
   for(int nb=0;nb<=FNAN[ga];nb++){
    int gb=natn[ga][nb],nbo=Spe_Total_CNO[WhatSpecies[gb]];
    for(int i=0;i<na;i++)for(int j=0;j<nbo;j++){
     int I=mp[ga]+i,J=mp[gb]+j;double s=first?olp[a][nb][i][j]:0;
     project_entry(buf,fields,j0,j1,t,I,J,zz(nh[0][a][nb][i][j],im[0][a][nb][i][j]),s);
     project_entry(buf,fields,j0,j1,t,I+B.N,J+B.N,zz(nh[1][a][nb][i][j],im[1][a][nb][i][j]),s);
     dcomplex ud=zz(nh[2][a][nb][i][j],nh[3][a][nb][i][j]+im[2][a][nb][i][j]);
     project_entry(buf,fields,j0,j1,t,I,J+B.N,ud,0);project_entry(buf,fields,j0,j1,t,J+B.N,I,cj(ud),0);
    }
   }
  }
  MPI_Reduce(buf,sum,(int)count,MPI_DOUBLE,MPI_SUM,0,mpi_comm_level1);memset(ct,0,p*4);
  for(int j=j0;j<j1;j++)for(int i=0;i<t->n;i++){int dest=((i/t->nb)%np_rows2)*np_cols2+(j/t->nb)%np_cols2;ct[dest]+=fields;}
  dp[0]=0;for(int q=1;q<p;q++)dp[q]=dp[q-1]+ct[q-1];memcpy(pos,dp,p*4);
  if(!r)for(int j=j0;j<j1;j++)for(int i=0;i<t->n;i++){int dest=((i/t->nb)%np_rows2)*np_cols2+(j/t->nb)%np_cols2;size_t from=((size_t)(j-j0)*t->n+i)*fields;if(from>count||(size_t)fields>count-from||pos[dest]<dp[dest]||pos[dest]>dp[dest]+ct[dest]-fields||(size_t)pos[dest]>cap-(size_t)fields)die("SOC complex root pack bounds");for(int q=0;q<fields;q++)packed[pos[dest]++]=sum[from+q];}
  if(!SocRecvCountValid(recv_alloc,recv_max,ct[r]))die("SOC Scatterv receive count exceeds actual allocation");if((size_t)ct[r]>recv_seen)recv_seen=(size_t)ct[r];
  MPI_Scatterv(packed,ct,dp,MPI_DOUBLE,recv,ct[r],MPI_DOUBLE,0,mpi_comm_level1);
  int off=0;for(int j=j0;j<j1;j++)if((j/t->nb)%np_cols2==my_pcol2)for(int i=0;i<t->n;i++)if((i/t->nb)%np_rows2==my_prow2){int il=(i/(t->nb*np_rows2))*t->nb+i%t->nb,jl=(j/(t->nb*np_cols2))*t->nb+j%t->nb;size_t ix=(size_t)jl*t->nr+il;if(ix>=(size_t)t->nr*t->nc||!SocRecvReadValid(recv_alloc,ct[r],(size_t)off,(size_t)fields))die("SOC complex unpack bounds");H[ix]=zz(recv[off],recv[off+1]);if(first)S[ix]=zz(recv[off+2],recv[off+3]);off+=fields;}
  if(off!=ct[r])die("tile distribution mismatch");
 }
 soc_recv_used_record(t,r,fields,compact_recv,recv_max,recv_alloc,recv_seen);
 free(mp);free(buf);free(sum);free(packed);free(recv);free(ct);free(dp);free(pos);CpuMemoryPhase("soc_complex_recv_freed");
}
static int global_index(int local,int nb,int coord,int nprocs){return (local/nb)*nb*nprocs+coord*nb+local%nb;}
static void dump_small(Sector*t,dcomplex*H,dcomplex*S,int step){
 if(t->n>512||B.mode!=1||step!=1||!getenv("OPENMX_CPU_C3_DUMP"))return;int r;MPI_Comm_rank(mpi_comm_level1,&r);size_t nn=(size_t)t->n*t->n;dcomplex*a=mem(2*nn,sizeof(dcomplex)),*sum=mem(2*nn,sizeof(dcomplex));
 for(int j=0;j<t->nc;j++){int gj=global_index(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++){int gi=global_index(i,t->nb,my_prow2,np_rows2);if(gi<t->n&&gj<t->n){a[(size_t)gi*t->n+gj]=H[(size_t)j*t->nr+i];a[nn+(size_t)gi*t->n+gj]=S[(size_t)j*t->nr+i];}}}
 MPI_Reduce(a,sum,(int)(nn*4),MPI_DOUBLE,MPI_SUM,0,mpi_comm_level1);if(!r){FILE*f=fopen("native_raw_HS.bin","wb");if(!f)die("dump file");int hd[4]={t->n,step,B.mode,t->nb};fwrite(hd,4,4,f);fwrite(sum,sizeof(dcomplex),nn*2,f);fclose(f);}free(a);free(sum);
}
static void expand_vectors(Sector*t,dcomplex*Q,int nev,dcomplex*ev){
 int r,p,lo,hi;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);bounds(B.ns,r,p,&lo,&hi);
 int*ct=mem(p,4),*rc=mem(p,4),*dp=mem(p,4),*rp=mem(p,4),*pos=mem(p,4);
 for(int b0=0;b0<nev;b0+=256){
  int b1=b0+256;if(b1>nev)b1=nev;memset(ct,0,p*4);
  for(int j=0;j<t->nc;j++){int b=global_index(j,t->nb,my_pcol2,np_cols2);if(b<b0||b>=b1)continue;int dst=owner(B.ns,p,t->ofs+b);for(int i=0;i<t->nr;i++)if(global_index(i,t->nb,my_prow2,np_rows2)<t->n)ct[dst]++;}
  MPI_Alltoall(ct,1,MPI_INT,rc,1,MPI_INT,mpi_comm_level1);dp[0]=rp[0]=0;for(int k=1;k<p;k++){dp[k]=dp[k-1]+ct[k-1];rp[k]=rp[k-1]+rc[k-1];}int ns=dp[p-1]+ct[p-1],nr=rp[p-1]+rc[p-1];Coeff*sb=mem(ns,sizeof(Coeff)),*rb=mem(nr,sizeof(Coeff));memcpy(pos,dp,p*4);
  for(int j=0;j<t->nc;j++){
   int b=global_index(j,t->nb,my_pcol2,np_cols2);if(b<b0||b>=b1)continue;int dst=owner(B.ns,p,t->ofs+b);
   for(int i=0;i<t->nr;i++){int a=global_index(i,t->nb,my_prow2,np_rows2);if(a>=t->n)continue;Coeff*z=sb+pos[dst]++;z->band=t->ofs+b;z->row=t->ofs+a;z->z=Q[(size_t)j*t->nr+i];}
  }
  for(int k=0;k<p;k++){if((int64_t)(ct[k]+dp[k]+rc[k]+rp[k])*sizeof(Coeff)>INT_MAX)die("coefficient MPI count overflow");ct[k]*=sizeof(Coeff);dp[k]*=sizeof(Coeff);rc[k]*=sizeof(Coeff);rp[k]*=sizeof(Coeff);}
  MPI_Alltoallv(sb,ct,dp,MPI_BYTE,rb,rc,rp,MPI_BYTE,mpi_comm_level1);
  for(int k=0;k<nr;k++){
   Coeff*z=rb+k;if(z->band<lo||z->band>=hi||z->row<0||z->row>=B.ns)die("coefficient distribution mismatch");size_t base=(size_t)(z->band-lo)*B.ns;
   for(int64_t a=B.cp[z->row];a<B.cp[z->row+1];a++)ev[base+B.ri[a]]=add(ev[base+B.ri[a]],mul(B.cv[a],z->z));
  }
  free(sb);free(rb);
 }
 free(ct);free(rc);free(dp);free(rp);free(pos);
}
typedef struct {double e;int source;} Eigen;
static int eigcmp(const void*pa,const void*pb){const Eigen*a=pa,*b=pb;if(a->e<b->e)return -1;if(a->e>b->e)return 1;return (a->source>b->source)-(a->source<b->source);}
static int final_owner(int p,int band,int*is2,int*ie2){for(int r=0;r<p;r++)if(is2[r]<=band+1&&band+1<=ie2[r])return r;die("missing final band owner");return 0;}
static void reorder(double*eval,double*ko,int maxn,int*is2,int*ie2,dcomplex*ev,dcomplex*scratch){
 int r,p,lo,hi;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);bounds(B.ns,r,p,&lo,&hi);Eigen*sort=mem(B.ns,sizeof(Eigen));int*inv=mem(B.ns,4);
 for(int i=0;i<B.ns;i++){if(!isfinite(eval[i]))die("nonfinite eigenvalue");sort[i].e=eval[i];sort[i].source=i;}qsort(sort,B.ns,sizeof(Eigen),eigcmp);for(int i=0;i<B.ns;i++){ko[i+1]=sort[i].e;inv[sort[i].source]=i;}
 au_record_state_sectors(inv,0);
 int *ct=mem(p,4),*rc=mem(p,4),*dp=mem(p,4),*rdp=mem(p,4),*pos=mem(p,4);
 for(int a=lo;a<hi;a++)if(inv[a]<maxn)ct[final_owner(p,inv[a],is2,ie2)]++;
 MPI_Alltoall(ct,1,MPI_INT,rc,1,MPI_INT,mpi_comm_level1);for(int k=1;k<p;k++){dp[k]=dp[k-1]+ct[k-1];rdp[k]=rdp[k-1]+rc[k-1];}int nsend=dp[p-1]+ct[p-1],nrecv=rdp[p-1]+rc[p-1];int*ids=mem(nsend,4),*rid=mem(nrecv,4),*local=mem(nsend,4);memcpy(pos,dp,p*4);
 for(int a=lo;a<hi;a++)if(inv[a]<maxn){int dst=final_owner(p,inv[a],is2,ie2),q=pos[dst]++;ids[q]=inv[a];local[q]=a-lo;}
 MPI_Alltoallv(ids,ct,dp,MPI_INT,rid,rc,rdp,MPI_INT,mpi_comm_level1);
 size_t capacity=(size_t)na_rows_max2*na_cols_max2,need=(size_t)nrecv*B.ns;int allocated=need>capacity;dcomplex*recv=allocated?mem(need,sizeof(dcomplex)):scratch;
 int transport=CpuSocExchangeBands(B.ns,ev,local,ct,dp,rc,rdp,recv,64u*1024u*1024u,mpi_comm_level1);
 if(transport)die("bounded band redistribution failed");
 if(nrecv!=ie2[r]-is2[r]+1)die("returned eigenstate count mismatch");
 for(int a=0;a<nrecv;a++){int dest=rid[a]+1-is2[r];if(dest<0||dest>=nrecv)die("sorted band owner mismatch");memcpy(ev+(size_t)dest*B.ns,recv+(size_t)a*B.ns,B.ns*sizeof(dcomplex));}
 if(allocated)free(recv);free(sort);free(inv);free(ct);free(rc);free(dp);free(rdp);free(pos);free(ids);free(rid);free(local);
 if(!r){printf("CPU_C3_BAND_REORDER max_packet_bytes=67108864 wide_pointer_offsets=1 transport=Sendrecv_ring\n");fflush(stdout);}
}
/* Symmetry-coordinate storage, with a physical band slot independent of the
 * logical state's sorted position. D3 partner1 aliases partner0's storage. */
static struct {
 int ready,local_bands,rank,np,first[4],count[4],stride[4];
 int *logical_sector,*logical_band,*sorted;
 int *rowptr[4],*column[4];dcomplex *value[4],*q[4];double *qr[4];
} C;
static int compact_sort_local=0,compact_cyclic=0;
void CpuC3_CompactBalanceConfigure(int sort,int cyclic){
 if((sort!=0&&sort!=1)||(cyclic!=0&&cyclic!=1))die("compact balance configuration requires0/1");
 if(C.ready||C.logical_sector||C.q[0]||C.qr[0])die("compact balance must be configured before storage initialization");
 if((sort||cyclic)&&!CpuC3_CompactEnabled())die("compact balance requires compact storage");
 compact_sort_local=sort;compact_cyclic=cyclic;
}
int CpuC3_CompactEnabled(void){
 const char*s=getenv("OPENMX_CPU_SOC_COMPACT");
 if(!s||!strcmp(s,"0"))return 0;
 if(strcmp(s,"1"))die("OPENMX_CPU_SOC_COMPACT must be exactly0 or1");
 return 1;
}
void CpuC3_CompactRelease(void){
 for(int k=0;k<4;k++){free(C.rowptr[k]);free(C.column[k]);free(C.value[k]);if(!(B.d3&&(k==1||(AU.tr&&k==3)))){free(C.q[k]);free(C.qr[k]);}}
 free(C.logical_sector);free(C.logical_band);free(C.sorted);memset(&C,0,sizeof(C));
}
int CpuC3_CompactLocalBands(void){return C.ready?C.local_bands:0;}
size_t CpuC3_CompactBytes(void){
 size_t n=0;for(int k=0;k<B.nsec;k++)if(!(B.d3&&(k==1||(AU.tr&&k==3)))){if(C.q[k])n+=(size_t)C.count[k]*B.sec[k].n*sizeof(dcomplex);if(C.qr[k])n+=(size_t)C.count[k]*B.sec[k].n*sizeof(double);}
 return n;
}
int CpuC3_CompactSortedBand(int local){return C.ready&&local>=0&&local<C.local_bands?C.sorted[local]:-1;}
static int compact_slot(int sector,int band){
 if(sector<0||sector>=B.nsec||C.stride[sector]<=0)return -1;
 int64_t d=(int64_t)band-C.first[sector];if(d<0||d%C.stride[sector])return -1;
 int64_t slot=d/C.stride[sector];return slot<C.count[sector]?(int)slot:-1;
}
static int compact_destination(int n,int p,int band){return compact_cyclic?band%p:owner(n,p,band);}
/* Half-open physical slot range whose global bands lie in [b0,b1). */
static void compact_packet_slots(int k,int b0,int b1,int*lo,int*hi){
 int64_t first=C.first[k],stride=C.stride[k];
 int64_t a=b0>first?((int64_t)b0-first+stride-1)/stride:0;
 int64_t b=b1>first?((int64_t)b1-first+stride-1)/stride:0;
 if(a>C.count[k])a=C.count[k];if(b>C.count[k])b=C.count[k];
 if(b<a)die("compact packet slot inversion");*lo=(int)a;*hi=(int)b;
}
int CpuC3_CompactPack(int row0,int rows,int band0,int bands,void*out,int ld){
 if(!C.ready||row0<0||rows<0||row0>B.ns||rows>B.ns-row0||band0<0||bands<0||
    band0>C.local_bands||bands>C.local_bands-band0||ld<rows||ld<0)return -1;
 if(!rows||!bands)return 0;
 if(!out||(size_t)ld>SIZE_MAX/sizeof(dcomplex)/(size_t)bands)return -1;
 dcomplex*dest=(dcomplex*)out;
 for(int b=0;b<bands;b++){
  int at=band0+b,k=C.logical_sector[at],src=(B.d3&&k==1)?0:(AU.tr&&k==3?2:k);
  int slot=compact_slot(k,C.logical_band[at]);if(slot<0)return -1;
  const dcomplex*q=C.q[src]?C.q[src]+(size_t)slot*B.sec[src].n:NULL;
  const double*qr=C.qr[src]?C.qr[src]+(size_t)slot*B.sec[src].n:NULL;
  if(!q&&!qr)return -1;
  for(int i=0;i<rows;i++){
   int row=row0+i;dcomplex z=zz(0,0);
   for(int j=C.rowptr[k][row];j<C.rowptr[k][row+1];j++){int col=C.column[k][j];dcomplex v=qr?zz(qr[col],0):q[col];if(AU.tr&&k==3)v=cj(v);z=add(z,mul(C.value[k][j],v));}
   dest[(size_t)b*ld+i]=z;
  }
 }
 return 0;
}
static void compact_begin(void){
 int r,p;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);CpuC3_CompactRelease();
 C.rank=r;C.np=p;if(B.mode!=3)die("compact wavefunctions require projected mode3");
 for(int k=0;k<B.nsec;k++){
  Sector*t=B.sec+k;
  if(compact_cyclic){C.first[k]=r;C.stride[k]=p;C.count[k]=t->n>r?1+(t->n-1-r)/p:0;}
  else{int lo,hi;bounds(t->n,r,p,&lo,&hi);C.first[k]=lo;C.stride[k]=1;C.count[k]=hi-lo;}
  if(C.count[k]>INT_MAX-C.local_bands)die("compact logical state count overflow");C.local_bands+=C.count[k];
  C.rowptr[k]=mem(B.ns+1,sizeof(int));int nz=0;
  for(int i=0;i<B.ns;i++){for(int a=B.rp[i];a<B.rp[i+1];a++)if(B.ci[a]>=t->ofs&&B.ci[a]<t->ofs+t->n)nz++;C.rowptr[k][i+1]=nz;}
  C.column[k]=mem(nz,sizeof(int));C.value[k]=mem(nz,sizeof(dcomplex));int at=0;
  for(int i=0;i<B.ns;i++)for(int a=B.rp[i];a<B.rp[i+1];a++)if(B.ci[a]>=t->ofs&&B.ci[a]<t->ofs+t->n){C.column[k][at]=B.ci[a]-t->ofs;C.value[k][at++]=B.rv[a];}
 }
 if(B.d3&&(C.first[0]!=C.first[1]||C.stride[0]!=C.stride[1]||C.count[0]!=C.count[1]))die("D3 partner layout differs");
 if(AU.tr&&(C.first[2]!=C.first[3]||C.stride[2]!=C.stride[3]||C.count[2]!=C.count[3]))die("TR partner layout differs");
 C.logical_sector=mem(C.local_bands,sizeof(int));C.logical_band=mem(C.local_bands,sizeof(int));C.sorted=mem(C.local_bands,sizeof(int));
 int at=0;for(int k=0;k<B.nsec;k++)for(int slot=0;slot<C.count[k];slot++){C.logical_sector[at]=k;C.logical_band[at++]=C.first[k]+slot*C.stride[k];}
}
typedef struct {int band,row;double value;} RealCoeff;
static void compact_store_typed(int sector,void*Q,int nev,int isreal){
 const char*profile_env=getenv("OPENMX_CPU_DM_BALANCE_PROFILE");int profile=0;
 if(profile_env&&strcmp(profile_env,"0")){if(strcmp(profile_env,"1"))die("balance profile requires0/1");profile=1;}
 double store_start=profile?MPI_Wtime():0.;
 int r,p;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);Sector*t=B.sec+sector;
 if(nev!=t->n||(B.d3&&sector==1)||(AU.tr&&sector==3)||isreal!=(AU.real&&sector==0))die("compact store requires one complete independent typed sector");
 if(isreal)C.qr[sector]=mem((size_t)C.count[sector]*t->n,sizeof(double));else C.q[sector]=mem((size_t)C.count[sector]*t->n,sizeof(dcomplex));
 size_t recsize=isreal?sizeof(RealCoeff):sizeof(Coeff);
 int*ct=mem(p,4),*rc=mem(p,4),*dp=mem(p,4),*rp=mem(p,4),*pos=mem(p,4);
 for(int b0=0;b0<nev;){
  int b1=nev-b0<256?nev:b0+256;memset(ct,0,p*sizeof(int));
  for(int j=0;j<t->nc;j++){int b=global_index(j,t->nb,my_pcol2,np_cols2);if(b<b0||b>=b1)continue;int dest=compact_destination(t->n,p,b);for(int i=0;i<t->nr;i++)if(global_index(i,t->nb,my_prow2,np_rows2)<t->n){if(ct[dest]==INT_MAX)die("compact count overflow");ct[dest]++;}}
  MPI_Alltoall(ct,1,MPI_INT,rc,1,MPI_INT,mpi_comm_level1);dp[0]=rp[0]=0;
  for(int i=1;i<p;i++){if(ct[i-1]>INT_MAX-dp[i-1]||rc[i-1]>INT_MAX-rp[i-1])die("compact displacement overflow");dp[i]=dp[i-1]+ct[i-1];rp[i]=rp[i-1]+rc[i-1];}
  if(ct[p-1]>INT_MAX-dp[p-1]||rc[p-1]>INT_MAX-rp[p-1])die("compact packet overflow");
  int ns=dp[p-1]+ct[p-1],nr=rp[p-1]+rc[p-1];void*send=mem(ns,recsize),*recv=mem(nr,recsize);memcpy(pos,dp,p*4);
  for(int j=0;j<t->nc;j++){int b=global_index(j,t->nb,my_pcol2,np_cols2);if(b<b0||b>=b1)continue;int dest=compact_destination(t->n,p,b);for(int i=0;i<t->nr;i++){int row=global_index(i,t->nb,my_prow2,np_rows2);if(row>=t->n)continue;int at=pos[dest]++;if(isreal){RealCoeff*z=(RealCoeff*)send+at;z->band=b;z->row=row;z->value=((double*)Q)[(size_t)j*t->nr+i];}else{Coeff*z=(Coeff*)send+at;z->band=b;z->row=row;z->z=((dcomplex*)Q)[(size_t)j*t->nr+i];}}}
  for(int i=0;i<p;i++){if(ct[i]>INT_MAX/(int)recsize||rc[i]>INT_MAX/(int)recsize||dp[i]>INT_MAX/(int)recsize||rp[i]>INT_MAX/(int)recsize)die("compact byte count overflow");ct[i]*=recsize;rc[i]*=recsize;dp[i]*=recsize;rp[i]*=recsize;}
  MPI_Alltoallv(send,ct,dp,MPI_BYTE,recv,rc,rp,MPI_BYTE,mpi_comm_level1);
  int lo,hi;compact_packet_slots(sector,b0,b1,&lo,&hi);size_t need=(size_t)(hi-lo)*t->n;unsigned char*seen=mem(need,1);
  if((size_t)nr!=need)die("compact coefficient count mismatch");
  for(int i=0;i<nr;i++){
   int band=isreal?((RealCoeff*)recv)[i].band:((Coeff*)recv)[i].band,row=isreal?((RealCoeff*)recv)[i].row:((Coeff*)recv)[i].row;
   int slot=compact_slot(sector,band);
   if(band<b0||band>=b1||slot<lo||slot>=hi||row<0||row>=t->n)die("compact coefficient owner mismatch");
   size_t ix=(size_t)(slot-lo)*t->n+row;if(seen[ix]++)die("duplicate compact coefficient");
   if(isreal)C.qr[sector][(size_t)slot*t->n+row]=((RealCoeff*)recv)[i].value;else C.q[sector][(size_t)slot*t->n+row]=((Coeff*)recv)[i].z;
  }
  for(size_t i=0;i<need;i++)if(!seen[i])die("missing compact coefficient");
  free(seen);free(send);free(recv);b0=b1;
 }
 free(ct);free(rc);free(dp);free(rp);free(pos);if(B.d3&&sector==0){C.q[1]=C.q[0];C.qr[1]=C.qr[0];}if(AU.tr&&sector==2)C.q[3]=C.q[2];
 if(profile){printf("COMPACT_STORE_PROFILE rank=%d sector=%d seconds=%.9f\n",r,sector,MPI_Wtime()-store_start);fflush(stdout);}
}
static void compact_store_sector(int sector,dcomplex*Q,int nev){compact_store_typed(sector,Q,nev,0);}
static void compact_store_real(int sector,double*Q,int nev){compact_store_typed(sector,Q,nev,1);}
typedef struct {int sector,band,sorted;} CompactLogical;
static int compact_logical_compare(const void*aa,const void*bb){const CompactLogical*a=aa,*b=bb;return (a->sorted>b->sorted)-(a->sorted<b->sorted);}
static void compact_sort(double*eval,double*ko){
 Eigen*sort=mem(B.ns,sizeof(Eigen));int*inv=mem(B.ns,sizeof(int));
 for(int i=0;i<B.ns;i++){if(!isfinite(eval[i]))die("nonfinite compact eigenvalue");sort[i].e=eval[i];sort[i].source=i;}qsort(sort,B.ns,sizeof(Eigen),eigcmp);
 for(int i=0;i<B.ns;i++){ko[i+1]=sort[i].e;inv[sort[i].source]=i+1;}
 au_record_state_sectors(inv,1);
 for(int i=0;i<C.local_bands;i++)C.sorted[i]=inv[B.sec[C.logical_sector[i]].ofs+C.logical_band[i]];
 if(compact_sort_local){
  CompactLogical*logical=mem(C.local_bands,sizeof(*logical));
  for(int i=0;i<C.local_bands;i++){logical[i].sector=C.logical_sector[i];logical[i].band=C.logical_band[i];logical[i].sorted=C.sorted[i];}
  qsort(logical,C.local_bands,sizeof(*logical),compact_logical_compare);
  for(int i=0;i<C.local_bands;i++){C.logical_sector[i]=logical[i].sector;C.logical_band[i]=logical[i].band;C.sorted[i]=logical[i].sorted;}
  free(logical);
 }
 C.ready=1;free(sort);free(inv);
 unsigned long long local=CpuC3_CompactBytes(),total=0;MPI_Reduce(&local,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,mpi_comm_level1);
 if(!C.rank){printf("CPU_SOC_COMPACT vectors_bytes=%llu logical_states=%d dense_AO_expansion=0 full_band_reorder=0\n",total,B.ns);printf("CPU_COMPACT_BALANCE sort=%d cyclic=%d sector_local_owner=1 D3_partner_alias_retained=1\n",compact_sort_local,compact_cyclic);fflush(stdout);}
}

extern void pdpotrf_(char*,int*,double*,int*,int*,int*,int*);
extern void pdsygst_(int*,char*,int*,double*,int*,int*,int*,double*,int*,int*,int*,double*,int*);
extern void pdtran_(int*,int*,double*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
extern void pdtrsm_(char*,char*,char*,char*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*);
extern void elpa_solve_evp_real_2stage_double_impl_(int*,int*,double*,int*,double*,double*,int*,int*,int*,int*,int*,int*);
/* Stats: max|Im H|, max|Im S|, sum(ReH²), sum(ImH²), sum(ReS²), sum(ImS²). */
static int au_imag_valid(const double*x){
 for(int i=0;i<6;i++)if(!isfinite(x[i])||x[i]<0)return 0;
 double h=sqrt(x[3])/fmax(sqrt(x[2]),DBL_MIN),s=sqrt(x[5])/fmax(sqrt(x[4]),DBL_MIN);
 return isfinite(h)&&isfinite(s)&&x[0]<=1e-10&&x[1]<=1e-12&&h<=1e-10&&s<=1e-10;
}
static void assemble_real(Sector*t,double*****nh,double*****im,double****olp,double*H,double*S){
 int r,p;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);int*mp=mem(atomnum+1,4);int off=0;
 for(int a=1;a<=atomnum;a++){mp[a]=off;off+=Spe_Total_CNO[WhatSpecies[a]];}
 const int tile=128,fields=4;size_t cap=(size_t)t->n*tile*fields;int compact_recv=soc_recv_mode();
 size_t recv_max=soc_recv_maximum(t,tile,2,r,p),recv_alloc=compact_recv?(recv_max?recv_max:1):cap/2,recv_seen=0;
 double*buf=mem(cap,8),*sum=r?NULL:mem(cap,8),*packed=r?NULL:mem(cap/2,8),*recv=mem(recv_alloc,8);
 soc_recv_alloc_record(t,r,2,compact_recv,recv_max,recv_alloc,cap/2);CpuMemoryPhase("soc_real_recv_allocated");
 int*ct=mem(p,4),*dp=mem(p,4),*pos=mem(p,4);double stats[6]={0};
 memset(H,0,(size_t)t->nr*t->nc*sizeof(double));memset(S,0,(size_t)t->nr*t->nc*sizeof(double));
 for(int j0=0;j0<t->n;j0+=tile){
  int j1=j0+tile;if(j1>t->n)j1=t->n;size_t count=(size_t)t->n*(j1-j0)*fields;if(count>INT_MAX)die("real tile count overflow");memset(buf,0,count*sizeof(double));
  for(int a=1;a<=Matomnum;a++){
   int ga=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga]];
   for(int b=0;b<=FNAN[ga];b++){int gb=natn[ga][b],nb=Spe_Total_CNO[WhatSpecies[gb]];
    for(int i=0;i<na;i++)for(int j=0;j<nb;j++){
     int I=mp[ga]+i,J=mp[gb]+j;double s=olp[a][b][i][j];
     project_entry(buf,fields,j0,j1,t,I,J,zz(nh[0][a][b][i][j],im[0][a][b][i][j]),s);
     project_entry(buf,fields,j0,j1,t,I+B.N,J+B.N,zz(nh[1][a][b][i][j],im[1][a][b][i][j]),s);
     dcomplex ud=zz(nh[2][a][b][i][j],nh[3][a][b][i][j]+im[2][a][b][i][j]);
     project_entry(buf,fields,j0,j1,t,I,J+B.N,ud,0);project_entry(buf,fields,j0,j1,t,J+B.N,I,cj(ud),0);
    }
   }
  }
  MPI_Reduce(buf,sum,(int)count,MPI_DOUBLE,MPI_SUM,0,mpi_comm_level1);memset(ct,0,p*4);
  for(int j=j0;j<j1;j++)for(int i=0;i<t->n;i++){int dest=((i/t->nb)%np_rows2)*np_cols2+(j/t->nb)%np_cols2;ct[dest]+=2;}
  dp[0]=0;for(int q=1;q<p;q++)dp[q]=dp[q-1]+ct[q-1];memcpy(pos,dp,p*4);
  if(!r)for(int j=j0;j<j1;j++)for(int i=0;i<t->n;i++){
   size_t from=((size_t)(j-j0)*t->n+i)*4;if(from>count||4>count-from)die("SOC real reduction read bounds");double hr=sum[from],hi=sum[from+1],sr=sum[from+2],si=sum[from+3];
   if(fabs(hi)>stats[0])stats[0]=fabs(hi);if(fabs(si)>stats[1])stats[1]=fabs(si);
   stats[2]+=hr*hr;stats[3]+=hi*hi;stats[4]+=sr*sr;stats[5]+=si*si;
   int dest=((i/t->nb)%np_rows2)*np_cols2+(j/t->nb)%np_cols2;if(pos[dest]<dp[dest]||pos[dest]>dp[dest]+ct[dest]-2||(size_t)pos[dest]>cap/2-2)die("SOC real root pack bounds");packed[pos[dest]++]=hr;packed[pos[dest]++]=sr;
  }
  if(!SocRecvCountValid(recv_alloc,recv_max,ct[r]))die("SOC Scatterv receive count exceeds actual allocation");if((size_t)ct[r]>recv_seen)recv_seen=(size_t)ct[r];
  MPI_Scatterv(packed,ct,dp,MPI_DOUBLE,recv,ct[r],MPI_DOUBLE,0,mpi_comm_level1);
  int at=0;for(int j=j0;j<j1;j++)if((j/t->nb)%np_cols2==my_pcol2)for(int i=0;i<t->n;i++)if((i/t->nb)%np_rows2==my_prow2){
   int il=(i/(t->nb*np_rows2))*t->nb+i%t->nb,jl=(j/(t->nb*np_cols2))*t->nb+j%t->nb;size_t ix=(size_t)jl*t->nr+il;if(ix>=(size_t)t->nr*t->nc||!SocRecvReadValid(recv_alloc,ct[r],(size_t)at,2))die("SOC real unpack bounds");H[ix]=recv[at++];S[ix]=recv[at++];
  }if(at!=ct[r])die("real tile distribution mismatch");
 }
 MPI_Bcast(stats,6,MPI_DOUBLE,0,mpi_comm_level1);
 if(!au_imag_valid(stats))die("assembled C2T gauge H/S imaginary guard failed");
 if(!r){printf("CPU_SOC_REAL_IMAG sector=0 passed=1 H_max=%.12e H_rel=%.12e S_max=%.12e S_rel=%.12e\n",stats[0],sqrt(stats[3])/fmax(sqrt(stats[2]),DBL_MIN),stats[1],sqrt(stats[5])/fmax(sqrt(stats[4]),DBL_MIN));fflush(stdout);}
 soc_recv_used_record(t,r,2,compact_recv,recv_max,recv_alloc,recv_seen);
 free(mp);free(buf);free(sum);free(packed);free(recv);free(ct);free(dp);free(pos);CpuMemoryPhase("soc_real_recv_freed");
}
static void antiunitary_solve_real(Sector*t,int step,double*eval,double*****nh,double*****im,double****olp,int rows,int cols,int world){
 if(t->nr<1||t->nc<1)die("real SOC requires nonempty local ELPA layout");
 int rank,one=1,info;MPI_Comm_rank(mpi_comm_level1,&rank);double a0=MPI_Wtime();size_t count=(size_t)t->nr*t->nc;
 double*S=mem(count,sizeof(double)),*H=mem(count,sizeof(double)),*Q=mem(count,sizeof(double));
 unsigned long long local=count*3*sizeof(double),total=0;MPI_Reduce(&local,&total,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,mpi_comm_level1);
 if(!rank){printf("CPU_SOC_TYPED_ALLOC sector=0 dtype=real64 work_S_H_Q_bytes=%llu caller_complex_buffers_used=0\n",total);fflush(stdout);}
 CpuMemoryPhase("real_sector_0_allocated");assemble_real(t,nh,im,olp,H,S);double a1=MPI_Wtime();
 pdpotrf_("L",&t->n,S,&one,&one,t->desc,&info);if(info)die("real overlap Cholesky failed");
 double scale=1.;pdsygst_(&one,"L",&t->n,H,&one,&one,t->desc,S,&one,&one,t->desc,&scale,&info);if(info||!isfinite(scale)||scale<=0)die("real generalized transform failed");
 for(int j=0;j<t->nc;j++){int gj=global_index(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++)if(global_index(i,t->nb,my_prow2,np_rows2)<gj)H[(size_t)j*t->nr+i]=0.;}
 double alpha=1.,beta=0.;pdtran_(&t->n,&t->n,&alpha,H,&one,&one,t->desc,&beta,Q,&one,&one,t->desc);
 for(int j=0;j<t->nc;j++){int gj=global_index(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++)if(global_index(i,t->nb,my_prow2,np_rows2)!=gj)H[(size_t)j*t->nr+i]+=Q[(size_t)j*t->nr+i];}
 double a2=MPI_Wtime();CpuMemoryPhase("real_sector_0_before_elpa");int nev=t->n,wanted=nev;
 elpa_solve_evp_real_2stage_double_impl_(&t->n,&nev,H,&t->nr,eval+t->ofs,Q,&t->nr,&t->nb,&t->nc,&rows,&cols,&world);if(nev!=wanted)die("real ELPA changed state count");
 MPI_Bcast(eval+t->ofs,t->n,MPI_DOUBLE,0,mpi_comm_level1);for(int i=0;i<t->n;i++)eval[t->ofs+i]*=scale;
 double a3=MPI_Wtime();CpuMemoryPhase("real_sector_0_after_elpa");
 pdtrsm_("L","L","T","N",&t->n,&wanted,&alpha,S,&one,&one,t->desc,Q,&one,&one,t->desc);compact_store_real(0,Q,wanted);
 free(S);free(H);free(Q);CpuMemoryPhase("real_sector_0_work_freed");double a4=MPI_Wtime();
 printf("CPU_C3_PHASE rank=%d step=%d sector=0 n=%d nev=%d assemble=%.6f transform=%.6f eigen=%.6f reconstruct=%.6f\n",rank,step,t->n,wanted,a1-a0,a2-a1,a3-a2,a4-a3);
 if(!rank){printf("CPU_SOC_REAL_SOLVE step=%d sector=0 n=%d dtype=real64\n",step,t->n);printf("CPU_D3_REUSE step=%d sector=1 source=0 states=%d distinct_eigenstates_retained=1\n",step,t->n);fflush(stdout);}
}

static void c2_apply(const dcomplex*v,dcomplex*out){
 /* U=[A,SA,P,Q]: S maps A->B, B->-A, P->iP, Q->-iQ. */
 if(!B.d3)die("C2 action requires D3 basis");int m=B.sz[0],p0=2*m,q0=p0+B.sz[2];dcomplex*c=mem(B.ns,sizeof(dcomplex));
 for(int j=0;j<B.ns;j++)for(int64_t k=B.cp[j];k<B.cp[j+1];k++)c[j]=add(c[j],mul(cj(B.cv[k]),v[B.ri[k]]));
 memset(out,0,(size_t)B.ns*sizeof(dcomplex));
 for(int j=0;j<B.ns;j++){
  dcomplex z=j<m?zz(-c[j+m].r,-c[j+m].i):j<p0?c[j-m]:j<q0?zz(-c[j].i,c[j].r):zz(c[j].i,-c[j].r);
  for(int64_t k=B.cp[j];k<B.cp[j+1];k++)out[B.ri[k]]=add(out[B.ri[k]],mul(B.cv[k],z));
 }
 free(c);
}
static void time_reverse(const dcomplex*v,dcomplex*out){
 for(int i=0;i<B.N;i++){out[i]=cj(v[i+B.N]);out[i+B.N]=zz(-v[i].r,v[i].i);}
}
static void project_vector(const dcomplex*v,dcomplex*out,int sector){
 /* The third C3 eigenspace includes both +/-i C2 eigenspaces in D3 mode. */
 memset(out,0,(size_t)B.ns*sizeof(dcomplex));int lo=sector?B.sz[0]:0;if(sector==2)lo+=B.sz[1];int hi=sector==2?B.ns:lo+B.sz[sector];
 for(int j=lo;j<hi;j++){dcomplex z=zz(0,0);for(int64_t k=B.cp[j];k<B.cp[j+1];k++)z=add(z,mul(cj(B.cv[k]),v[B.ri[k]]));for(int64_t k=B.cp[j];k<B.cp[j+1];k++)out[B.ri[k]]=add(out[B.ri[k]],mul(B.cv[k],z));}
}
static void raw_apply(const dcomplex*v,dcomplex*hv,dcomplex*sv,double*****nh,double*****im,double****olp){
 dcomplex*work=mem((size_t)B.ns*2,sizeof(dcomplex));int*mp=mem(atomnum+1,4);int off=0;for(int a=1;a<=atomnum;a++){mp[a]=off;off+=Spe_Total_CNO[WhatSpecies[a]];}
 for(int a=1;a<=Matomnum;a++){int ga=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga]];for(int nb=0;nb<=FNAN[ga];nb++){int gb=natn[ga][nb],nbo=Spe_Total_CNO[WhatSpecies[gb]];for(int i=0;i<na;i++)for(int j=0;j<nbo;j++){
  int I=mp[ga]+i,J=mp[gb]+j;dcomplex u=zz(nh[0][a][nb][i][j],im[0][a][nb][i][j]),d=zz(nh[1][a][nb][i][j],im[1][a][nb][i][j]),ud=zz(nh[2][a][nb][i][j],nh[3][a][nb][i][j]+im[2][a][nb][i][j]);double s=olp[a][nb][i][j];
  work[I]=add(work[I],add(mul(u,v[J]),mul(ud,v[J+B.N])));work[I+B.N]=add(work[I+B.N],mul(d,v[J+B.N]));work[J+B.N]=add(work[J+B.N],mul(cj(ud),v[I]));
  work[B.ns+I]=add(work[B.ns+I],zz(s*v[J].r,s*v[J].i));work[B.ns+I+B.N]=add(work[B.ns+I+B.N],zz(s*v[J+B.N].r,s*v[J+B.N].i));
 }}}
 dcomplex*all=mem((size_t)B.ns*2,sizeof(dcomplex));MPI_Allreduce(work,all,B.ns*4,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);memcpy(hv,all,B.ns*sizeof(dcomplex));memcpy(sv,all+B.ns,B.ns*sizeof(dcomplex));free(all);free(work);free(mp);
}
static double norm2(dcomplex z){return z.r*z.r+z.i*z.i;}
static void c3_projected_apply(const dcomplex*v,dcomplex*hv,dcomplex*sv,double*****nh,double*****im,double****olp){
 dcomplex*pv=mem(B.ns,sizeof(dcomplex)),*h=mem(B.ns,sizeof(dcomplex)),*s=mem(B.ns,sizeof(dcomplex)),*ph=mem(B.ns,sizeof(dcomplex)),*ps=mem(B.ns,sizeof(dcomplex));
 memset(hv,0,(size_t)B.ns*sizeof(dcomplex));memset(sv,0,(size_t)B.ns*sizeof(dcomplex));
 for(int sec=0;sec<3;sec++){
  project_vector(v,pv,sec);raw_apply(pv,h,s,nh,im,olp);project_vector(h,ph,sec);project_vector(s,ps,sec);
  for(int i=0;i<B.ns;i++){hv[i]=add(hv[i],ph[i]);sv[i]=add(sv[i],ps[i]);}
 }
 free(pv);free(h);free(s);free(ph);free(ps);
}
static void probe_symmetry(int step,double*****nh,double*****im,double****olp){
 if(!B.d3)return;int r;MPI_Comm_rank(mpi_comm_level1,&r);
 dcomplex*v=mem(B.ns,sizeof(dcomplex)),*h=mem(B.ns,sizeof(dcomplex)),*s=mem(B.ns,sizeof(dcomplex)),*pv=mem(B.ns,sizeof(dcomplex)),*hp=mem(B.ns,sizeof(dcomplex)),*sp=mem(B.ns,sizeof(dcomplex)),*ph=mem(B.ns,sizeof(dcomplex)),*ps=mem(B.ns,sizeof(dcomplex));
 for(int i=0;i<B.ns;i++)v[i]=zz(sin((i+1)*.317)/sqrt((double)B.ns),cos((i+1)*.173)/sqrt((double)B.ns));
 raw_apply(v,h,s,nh,im,olp);
 for(int kind=0;kind<2;kind++){
  if(kind==0){c2_apply(v,pv);c2_apply(h,ph);c2_apply(s,ps);}else{time_reverse(v,pv);time_reverse(h,ph);time_reverse(s,ps);}
  raw_apply(pv,hp,sp,nh,im,olp);double eh=0,es=0,dh0=0,dh1=0,ds0=0,ds1=0,mh=0,ms=0;
  for(int i=0;i<B.ns;i++){dcomplex a=add(hp[i],zz(-ph[i].r,-ph[i].i)),b=add(sp[i],zz(-ps[i].r,-ps[i].i));eh+=norm2(a);es+=norm2(b);dh0+=norm2(hp[i]);dh1+=norm2(ph[i]);ds0+=norm2(sp[i]);ds1+=norm2(ps[i]);if(hypot(a.r,a.i)>mh)mh=hypot(a.r,a.i);if(hypot(b.r,b.i)>ms)ms=hypot(b.r,b.i);}
  if(!r){printf("CPU_SYMMETRY_PROBE step=%d mode=%d kind=%s H_rel=%.12e H_max=%.12e S_rel=%.12e S_max=%.12e\n",step,B.mode,kind==0?"C2":"T",sqrt(eh)/fmax(sqrt(dh0)+sqrt(dh1),1e-300),mh,sqrt(es)/fmax(sqrt(ds0)+sqrt(ds1),1e-300),ms);fflush(stdout);}
 }
 free(v);free(h);free(s);free(pv);free(hp);free(sp);free(ph);free(ps);
}
static void operator_bounds(double*****nh,double*****im,double****olp,double*out){
 double*w=mem((size_t)B.ns*2,8),*all=mem((size_t)B.ns*2,8);int*mp=mem(atomnum+1,4),off=0;
 for(int a=1;a<=atomnum;a++){mp[a]=off;off+=Spe_Total_CNO[WhatSpecies[a]];}
 for(int a=1;a<=Matomnum;a++){int ga=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga]];for(int k=0;k<=FNAN[ga];k++){int gb=natn[ga][k],nb=Spe_Total_CNO[WhatSpecies[gb]];for(int i=0;i<na;i++)for(int j=0;j<nb;j++){
  int I=mp[ga]+i,J=mp[gb]+j;double ud=hypot(nh[2][a][k][i][j],nh[3][a][k][i][j]+im[2][a][k][i][j]);w[I]+=hypot(nh[0][a][k][i][j],im[0][a][k][i][j])+ud;w[I+B.N]+=hypot(nh[1][a][k][i][j],im[1][a][k][i][j]);w[J+B.N]+=ud;w[B.ns+I]+=fabs(olp[a][k][i][j]);w[B.ns+I+B.N]+=fabs(olp[a][k][i][j]);
 }}}
 MPI_Allreduce(w,all,B.ns*2,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);out[0]=out[1]=0;for(int i=0;i<B.ns;i++){if(all[i]>out[0])out[0]=all[i];if(all[B.ns+i]>out[1])out[1]=all[B.ns+i];}free(w);free(all);free(mp);
}
static void validate_native(int step,int maxn,double*ko,int*is2,int*ie2,dcomplex*ev,double*****nh,double*****im,double****olp){
 int r,p;MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);double opnorm[2];operator_bounds(nh,im,olp,opnorm);int pick[7],npick=au_validation_picks(maxn,pick);dcomplex*vs=mem((size_t)B.ns*npick,sizeof(dcomplex)),*ss=mem((size_t)B.ns*npick,sizeof(dcomplex));
 dcomplex*pv=mem(B.ns,sizeof(dcomplex)),*h=mem(B.ns,sizeof(dcomplex)),*s=mem(B.ns,sizeof(dcomplex)),*ph=mem(B.ns,sizeof(dcomplex)),*ps=mem(B.ns,sizeof(dcomplex)),*rh=mem(B.ns,sizeof(dcomplex)),*rs=mem(B.ns,sizeof(dcomplex)),*sh=mem(B.ns,sizeof(dcomplex)),*sS=mem(B.ns,sizeof(dcomplex));
 for(int k=0;k<npick;k++){
  int band=pick[k],root=final_owner(p,band-1,is2,ie2);dcomplex*v=vs+(size_t)k*B.ns;memset(rh,0,B.ns*sizeof(dcomplex));memset(rs,0,B.ns*sizeof(dcomplex));memset(sh,0,B.ns*sizeof(dcomplex));memset(sS,0,B.ns*sizeof(dcomplex));
  if(C.ready){
   /* Diagnostic-only full vector: exactly one compact owner contributes. */
   int found=0,all=0;memset(v,0,(size_t)B.ns*sizeof(dcomplex));
   for(int b=0;b<C.local_bands;b++)if(C.sorted[b]==band){if(CpuC3_CompactPack(0,B.ns,b,1,v,B.ns))die("compact verification pack failed");found++;}
   MPI_Allreduce(&found,&all,1,MPI_INT,MPI_SUM,mpi_comm_level1);if(all!=1)die("compact sampled state lacks unique owner");
   if(B.ns>INT_MAX/2)die("compact verifier MPI count overflow");MPI_Allreduce(MPI_IN_PLACE,v,B.ns*2,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
  }else{if(r==root)memcpy(v,ev+(size_t)(band-is2[r])*B.ns,B.ns*sizeof(dcomplex));MPI_Bcast(v,B.ns*2,MPI_DOUBLE,root,mpi_comm_level1);}
  if(B.mode==1||B.mode==4){raw_apply(v,rh,rs,nh,im,olp);memcpy(sh,rh,B.ns*sizeof(dcomplex));memcpy(sS,rs,B.ns*sizeof(dcomplex));}
  else for(int sec=0;sec<3;sec++){project_vector(v,pv,sec);raw_apply(pv,h,s,nh,im,olp);project_vector(h,ph,sec);project_vector(s,ps,sec);for(int i=0;i<B.ns;i++){rh[i]=add(rh[i],h[i]);rs[i]=add(rs[i],s[i]);sh[i]=add(sh[i],ph[i]);sS[i]=add(sS[i],ps[i]);}}
  if(B.d3&&B.mode==3){
   /* P_D3(X)=1/2[P_C3(X)+S^dagger P_C3(X)S], with S^dagger=-S. */
   c2_apply(v,pv);c3_projected_apply(pv,ph,ps,nh,im,olp);c2_apply(ph,h);c2_apply(ps,s);
   for(int i=0;i<B.ns;i++){sh[i]=zz(.5*(sh[i].r-h[i].r),.5*(sh[i].i-h[i].i));sS[i]=zz(.5*(sS[i].r-s[i].r),.5*(sS[i].i-s[i].i));}
  }
  double raw=0,sym=0,hn=0,sn=0;for(int i=0;i<B.ns;i++){double e=ko[band];raw+=norm2(add(rh[i],zz(-e*rs[i].r,-e*rs[i].i)));sym+=norm2(add(sh[i],zz(-e*sS[i].r,-e*sS[i].i)));hn+=norm2(sh[i]);sn+=norm2(sS[i]);}
  double actionden=sqrt(hn)+fabs(ko[band])*sqrt(sn),vn=0;for(int i=0;i<B.ns;i++)vn+=norm2(v[i]);double den=(opnorm[0]+fabs(ko[band])*opnorm[1])*sqrt(vn);
  if(!isfinite(den)||den<=0||!isfinite(actionden)||actionden<=0)die("bad residual denominator");double action=sqrt(sym)/actionden;raw=sqrt(raw)/den;sym=sqrt(sym)/den;
  if(!r){printf("CPU_C3_VERIFY step=%d band=%d eigen=%.12e projected_residual=%.12e raw_residual=%.12e action_residual=%.12e H_norm_bound=%.12e S_norm_bound=%.12e vector_norm=%.12e sector=%d\n",step,band,ko[band],sym,raw,action,opnorm[0],opnorm[1],sqrt(vn),AU.sector_by_state[band-1]);fflush(stdout);}
  if(!isfinite(sym)||sym>1e-10)die("projected generalized eigenpair residual failed");memcpy(ss+(size_t)k*B.ns,sS,B.ns*sizeof(dcomplex));
 }
 double orth=0;for(int a=0;a<npick;a++)for(int b=0;b<npick;b++){dcomplex z=zz(0,0);for(int i=0;i<B.ns;i++)z=add(z,mul(cj(vs[(size_t)a*B.ns+i]),ss[(size_t)b*B.ns+i]));if(a==b)z.r-=1;double err=hypot(z.r,z.i);if(err>orth)orth=err;}
 if(!isfinite(orth)||orth>1e-8)die("sampled S-orthonormality failed");if(!r){printf("CPU_C3_VERIFY step=%d sampled_S_orth=%.12e\n",step,orth);fflush(stdout);}
 free(vs);free(ss);free(pv);free(h);free(s);free(ph);free(ps);free(rh);free(rs);free(sh);free(sS);
}
void CpuC3_Solve(int step,int maxn,double*ko,double*****nh,double*****im,double****overlap,int*is2,int*ie2,dcomplex*sf,dcomplex*H,dcomplex*Q,dcomplex*ev){
 int config=CpuC3_Mode(),mode=config,r,p,one=1,info,lo,hi,compact=CpuC3_CompactEnabled();
 if(AU.tr&&(CpuSOCNative_TRRequested()!=1||CpuSOCNative_TRPassed()!=1))die("antiunitary requires actual raw H/S time-reversal guard");
 if(compact){
  const char*dm=getenv("OPENMX_CPU_DM_BLAS3");
  if(config!=3||step!=1||Solver!=2||SpinP_switch!=3||SO_switch!=1||DFTSCF_loop!=1||MD_switch!=0||Cnt_switch||
     level_fileout>=2||CDDF_on||LNAO_calc_flag||LNBO_calc_flag||NBO_switch||!dm||strcmp(dm,"1")||Band_disp_switch||MO_fileout||Dos_fileout||DosGauss_fileout||
     Wannier_Func_Calc||Wannier90_fileout||CWF_fileout_flag||flag_export_DM||pop_anal_aow_flag||
     unfold_electronic_band||Polarization_ABC_flag||CWF_Calc||xanes_gs_fileout||xanes_calc||cal_partial_charge)
   die("compact vectors require fixed one-step projected SOC, BLAS3 owner DM and no state-based exports");
 }
 if(config==5){const char*from=getenv("OPENMX_CPU_C3_FULL_FROM");int release=from?atoi(from):3;if(release<1)die("invalid full-solver release step");mode=step>=release?4:3;}
 if(config==6){
  if(step==1 && (MD_switch!=0 || SCF_Control_Temp!=0 || Scf_RestartFromFile!=0))
   die("adaptive C3 requires Nomd, temperature-controlled mixing off, and restart reading off");
  if(step==1){if(B.loaded)die("adaptive prototype requires a fresh process for each geometry");force_full=0;full_steps=0;}
  const char*tol=getenv("OPENMX_CPU_C3_RELEASE_NORM");double release=tol?atof(tol):1e-3;
  if(!isfinite(release)||release<=0)die("invalid SCF release norm");
  if(step>=DFTSCF_loop||(step>1&&NormRD[0]<release))force_full=1;
  mode=force_full?4:3;
 }
 const char*au_check=getenv("OPENMX_CPU_C3_CHECK");AU.check=au_check&&atoi(au_check);
 active_config=config;last_step=step;MPI_Bcast(&mode,1,MPI_INT,0,mpi_comm_level1);if(mode==4)full_steps++;else full_steps=0;
 int rebuild=(step==1||!B.loaded||B.mode!=mode);MPI_Comm_rank(mpi_comm_level1,&r);MPI_Comm_size(mpi_comm_level1,&p);double start=MPI_Wtime();c3_init(mode);bounds(B.ns,r,p,&lo,&hi);CpuMemoryPhase("before_vector_init");if(compact)compact_begin();else memset(ev,0,(size_t)(hi-lo)*B.ns*sizeof(dcomplex));CpuMemoryPhase("after_vector_init");double*eval=mem(B.ns,sizeof(double));int rows,cols,world;
 MPI_Comm rowcomm,colcomm;
 MPI_Comm_split(mpi_comm_level1,my_pcol2,my_prow2,&rowcomm);MPI_Comm_split(mpi_comm_level1,my_prow2,my_pcol2,&colcomm);rows=MPI_Comm_c2f(rowcomm);cols=MPI_Comm_c2f(colcomm);world=MPI_Comm_c2f(mpi_comm_level1);
 for(int k=0;k<B.nsec;k++){
  if(B.d3&&mode==3&&(k==1||(AU.tr&&k==3)))continue; /* Logical partners retain distinct states. */
  if(AU.real&&k==0){antiunitary_solve_real(B.sec,step,eval,nh,im,overlap,rows,cols,world);memcpy(eval+B.sec[1].ofs,eval,(size_t)B.sec[0].n*sizeof(double));continue;}
  Sector*t=B.sec+k;dcomplex*savedH=H,*savedQ=Q;
  if(AU.real){if(t->nr<1||t->nc<1)die("typed SOC requires nonempty local ELPA layout");H=mem((size_t)t->nr*t->nc,sizeof(dcomplex));Q=mem((size_t)t->nr*t->nc,sizeof(dcomplex));}
  dcomplex*S=AU.real?mem((size_t)t->nr*t->nc,sizeof(dcomplex)):sf+t->soff;double a0=MPI_Wtime();assemble(t,rebuild,nh,im,overlap,H,S);double a1=MPI_Wtime();dump_small(t,H,S,step);
  if(rebuild){pzpotrf_("L",&t->n,S,&one,&one,t->desc,&info);if(info)die("overlap Cholesky failed");}
  double scale=1.;pzhegst_(&one,"L",&t->n,H,&one,&one,t->desc,S,&one,&one,t->desc,&scale,&info);if(info||!isfinite(scale)||scale<=0)die("generalized transform failed");
  /* PZHEGST defines the selected lower triangle only. Fill the other half
   * explicitly before calling ELPA; never pass stale upper-triangle entries. */
  for(int j=0;j<t->nc;j++){int gj=global_index(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++){int gi=global_index(i,t->nb,my_prow2,np_rows2);if(gi<gj)H[(size_t)j*t->nr+i]=zz(0,0);}}
  dcomplex alpha=zz(1,0),beta=zz(0,0);pztranc_(&t->n,&t->n,&alpha,H,&one,&one,t->desc,&beta,Q,&one,&one,t->desc);
  for(int j=0;j<t->nc;j++){int gj=global_index(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++){int gi=global_index(i,t->nb,my_prow2,np_rows2);size_t at=(size_t)j*t->nr+i;if(gi!=gj)H[at]=add(H[at],Q[at]);else H[at].i=0;}}
  double a2=MPI_Wtime();int nev=(mode==4?maxn:t->n),wanted=nev;elpa_solve_evp_complex_2stage_double_impl_(&t->n,&nev,H,&t->nr,eval+t->ofs,Q,&t->nr,&t->nb,&t->nc,&rows,&cols,&world);if(nev!=wanted)die("ELPA changed state count");MPI_Bcast(eval+t->ofs,t->n,MPI_DOUBLE,0,mpi_comm_level1);for(int a=0;a<wanted;a++)eval[t->ofs+a]*=scale;for(int a=wanted;a<t->n;a++)eval[t->ofs+a]=DBL_MAX;double a3=MPI_Wtime();
  pztrsm_("L","L","C","N",&t->n,&wanted,&alpha,S,&one,&one,t->desc,Q,&one,&one,t->desc);if(compact)compact_store_sector(k,Q,wanted);else expand_vectors(t,Q,wanted,ev);
  if(B.d3&&mode==3&&k==0){
   Sector*partner=B.sec+1;memcpy(eval+partner->ofs,eval+t->ofs,(size_t)t->n*sizeof(double));if(!compact)expand_vectors(partner,Q,wanted,ev);
   if(!r){printf("CPU_D3_REUSE step=%d sector=1 source=0 states=%d distinct_eigenstates_retained=1\n",step,t->n);fflush(stdout);}
  }
  if(AU.tr&&k==2){memcpy(eval+B.sec[3].ofs,eval+t->ofs,(size_t)t->n*sizeof(double));C.q[3]=C.q[2];if(!r){printf("CPU_TR_REUSE step=%d sector=3 source=2 distinct_eigenstates_retained=1 coefficients=conjugate_alias\n",step);fflush(stdout);}}
  if(AU.real){free(S);free(H);free(Q);H=savedH;Q=savedQ;}
  char memory_stage[80];snprintf(memory_stage,sizeof(memory_stage),"after_sector_%d",k);CpuMemoryPhase(memory_stage);
  double a4=MPI_Wtime();
  printf("CPU_C3_PHASE rank=%d step=%d sector=%d n=%d nev=%d assemble=%.6f transform=%.6f eigen=%.6f reconstruct=%.6f\n",r,step,k,t->n,wanted,a1-a0,a2-a1,a3-a2,a4-a3);fflush(stdout);
 }
 CpuMemoryPhase("before_vector_order");
 if(compact)compact_sort(eval,ko);else reorder(eval,ko,maxn,is2,ie2,ev,Q);
 CpuMemoryPhase("after_vector_order");
 const char*check=getenv("OPENMX_CPU_C3_CHECK");if(check&&atoi(check))validate_native(step,maxn,ko,is2,ie2,ev,nh,im,overlap);
 const char*probe=getenv("OPENMX_CPU_SYMMETRY_PROBE");if(probe&&atoi(probe))probe_symmetry(step,nh,im,overlap);
 free(AU.sector_by_state);AU.sector_by_state=NULL;free(eval);MPI_Comm_free(&rowcomm);MPI_Comm_free(&colcomm);double elapsed=MPI_Wtime()-start;printf("CPU_C3 step=%d rank=%d mode=%d seconds=%.6f maxn=%d\n",step,r,mode,elapsed,maxn);fflush(stdout);
}
