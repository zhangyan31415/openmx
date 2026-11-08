#include "cpu_reclaim.h"
/* Optional fixed-geometry Gamma spinless D3 warmup. Stock raw solver resumes
 * before convergence. One additional distributed real matrix, no replicated
 * dense U/H/S. Original occupations, density, EDM and force routines retained. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mpi.h"
#include "openmx_common.h"
#include "cpu_nsoc_sym.h"
#include "cpu_nsoc_scratch.h"
#include "cpu_nsoc_sparse.h"
#include "cpu_nsoc_owner_openmx.h"
#include "native_panel_solver.h"
extern int flag_export_DM;
static int direct_active=0;
static int direct_mode=-1,direct_step=0;
static double direct_reference_s_defect=0.;
static void direct_assemble(int step,int kind,double*full,int*desc,double*work);
static void direct_release(void);
static void native_state_release(void);

extern int numroc_(int*,int*,int*,int*,int*);
extern void descinit_(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void pdgemm_(char*,char*,int*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
extern void pdgemr2d_(int*,int*,double*,int*,int*,int*,double*,int*,int*,int*,int*);
extern void pdpotrf_(char*,int*,double*,int*,int*,int*,int*);
extern void pdtrsm_(char*,char*,char*,char*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*);
extern void pdtran_(int*,int*,double*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
extern void elpa_solve_evp_real_2stage_double_impl_(int*,int*,double*,int*,double*,double*,int*,int*,int*,int*,int*,int*);

extern void pzpotrf_(char*,int*,dcomplex*,int*,int*,int*,int*);
extern void pztrsm_(char*,char*,char*,char*,int*,int*,dcomplex*,dcomplex*,int*,int*,int*,dcomplex*,int*,int*,int*);
extern void pztranc_(int*,int*,dcomplex*,dcomplex*,int*,int*,int*,dcomplex*,dcomplex*,int*,int*,int*);
extern void elpa_solve_evp_complex_2stage_double_impl_(int*,int*,dcomplex*,int*,double*,dcomplex*,int*,int*,int*,int*,int*,int*);
typedef struct {int n,ofs,nr,nc,nb,desc[9];size_t soff;} Sector;
typedef struct {int row,col;double value;} Entry;
typedef struct {double value;int source;} Eigen;
#define NSYM_MAX_SECTORS 48
static struct {int loaded,version,N,natoms,nnz,nsec,sz[NSYM_MAX_SECTORS],*perm,*c2perm,*ri,*rankmap;int64_t *cp;double *cv,*U,*partial_eval_cache;size_t cap,active,scap,hcap,hoff[NSYM_MAX_SECTORS];Sector sec[NSYM_MAX_SECTORS];double shift[3],rot[9],c2shift[3];int partial_cache[NSYM_MAX_SECTORS],partial_cache_valid;} B;
/* SECTOR_PARALLEL_PATCH_BEGIN */
static void sector_parallel_release(void);
/* SECTOR_PARALLEL_PATCH_END */
static CpuNSocEigenScratch eigen_scratch={0};
static int borrowed_eigen_scratch=0;
static int representative_certificate_ready=0,representative_certificate_n=0,representative_certificate_maxn=0;
static int representative_tail_ready=0,representative_cut_count=0,*representative_cut_ids=NULL,representative_cut_cap=0;
static int prepared=0,last_step=0,warm=0,force_full=0,raw_rebuild=0,full_steps=0,first_raw_step=0;

static void die(const char*s){int r;MPI_Comm_rank(mpi_comm_level1,&r);fprintf(stderr,"CPU_NSOC_SYM ERROR rank=%d %s\n",r,s);fflush(stderr);MPI_Abort(mpi_comm_level1,95);abort();}
static void *mem(size_t n,size_t s){if(s&&n>SIZE_MAX/s)die("allocation overflow");void*p=calloc(n?n:1,s);if(!p)die("allocation failed");return p;}
static void rd(void*p,size_t s,size_t n,FILE*f){if(fread(p,s,n,f)!=n)die("short real basis file");}
static double option(const char*k,double fallback){const char*p=getenv(k);if(!p||!*p)return fallback;char*end;double q=strtod(p,&end);if(end==p||*end||!isfinite(q))die("invalid numeric environment setting");return q;}
static struct {double*coeff;size_t capacity;int step,first,last,returned;} native_return;
static int native_return_mode(void){
 static int cached=-1;if(cached>=0)return cached;
 double q=option("OPENMX_CPU_NSOC_NATIVE_RETURN",0);int v=(int)q,lo,hi;
 if(q!=0&&q!=1&&q!=2&&q!=3&&q!=4)die("NATIVE_RETURN must be 0, 1, 2(audit), 3(provider) or 4(provider audit)");
 MPI_Allreduce(&v,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);MPI_Allreduce(&v,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);
 if(lo!=hi)die("NATIVE_RETURN differs across ranks");cached=v;return cached;
}
int CpuNSym_NativeReturnEnabled(void){return native_return_mode()!=0;}
void CpuNSym_SetNativeReturn(int step,double*coeff,size_t capacity,int first,int last){
 if(!CpuNSym_NativeReturnEnabled())return;
 if(step<1||!coeff||last<first-1||(last>=first&&first<1))die("invalid native coefficient registration");
 native_return.coeff=coeff;native_return.capacity=capacity;native_return.step=step;
 native_return.first=first;native_return.last=last;native_return.returned=0;
}
int CpuNSym_NativeReturned(int step){return CpuNSym_NativeReturnEnabled()&&native_return.returned==step;}
static int rank0(void){int r;MPI_Comm_rank(mpi_comm_level1,&r);return r==0;}
static int raw_reset_at(void){
 double q=option("OPENMX_CPU_NSOC_RAW_RESET_AT",0);
 if(q<0||q>INT_MAX||q!=(int)q)die("RAW_RESET_AT must be a nonnegative integer");
 return (int)q;
}
int CpuNSym_FullProjected(void){
 double q=option("OPENMX_CPU_NSOC_FULL_PROJECT",0);
 if(q!=0&&q!=1)die("FULL_PROJECT must be0 or1");
 return q!=0;
}
int CpuNSym_Enabled(void){
 double q=option("OPENMX_CPU_NSOC_SYM",0);int control=raw_reset_at();
 if(q!=0&&q!=1)die("OPENMX_CPU_NSOC_SYM must be0 or1");
 if(q&&control)die("RAW_RESET_AT attribution control cannot be combined with symmetry");
 if(CpuNSym_FullProjected()){
  double cap=option("OPENMX_CPU_NSOC_WARM_STEPS",0),reset=option("OPENMX_CPU_NSOC_RESET_DIIS",0);
  if(!q)die("FULL_PROJECT requires OPENMX_CPU_NSOC_SYM=1");
  if(control||cap!=0||reset!=0)die("FULL_PROJECT requires RAW_RESET_AT=0 WARM_STEPS=0 RESET_DIIS=0");
 }
 if((q||control) && (SpinP_switch!=0||SO_switch!=0||Solver!=2||Cnt_switch||MD_switch!=0||SCF_Control_Temp||Scf_RestartFromFile))die("requires nonmagnetic NSOC Gamma Cluster, Nomd, uncontracted basis and restart off");
 return q!=0;
}
/* Fixed-geometry, fully projected NSOC solves do not read the previous
 * iteration's EVec1. Reuse that caller-owned storage for the U/sector scratch
 * only inside CpuNSym_Solve. Native Hs-to-EVec1 redistribution follows its last
 * scratch read. Distinct H/S/Cs buffers and every eigensolver are unchanged. */
static int compact_eigen_scratch(void){
 double q=option("OPENMX_CPU_NSOC_EIGEN_SCRATCH",0);
 if(q!=0&&q!=1)die("EIGEN_SCRATCH must be 0 or 1");
 if(q&&(!CpuNSym_Enabled()||!CpuNSym_FullProjected()))die("EIGEN_SCRATCH requires full projected NSOC");
 return q!=0;
}
void CpuNSym_ReserveEigenScratch(double **coeff,size_t native_capacity,size_t matrix_capacity){
 if(!compact_eigen_scratch())return;
 if(B.loaded||eigen_scratch.ptr)die("EIGEN_SCRATCH must be reserved once before the first solve");
 if(CpuNSocScratch_Reserve(&eigen_scratch,coeff,native_capacity,matrix_capacity))die("EIGEN_SCRATCH capacity allocation failed");
 int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
 printf("CPU_NSOC_EIGEN_SCRATCH_RESERVE rank=%d native_doubles=%zu matrix_doubles=%zu reserved_doubles=%zu extra_doubles=%zu owner=EVec1\n",rank,native_capacity,matrix_capacity,eigen_scratch.capacity,eigen_scratch.capacity-native_capacity);fflush(stdout);
}
size_t CpuNSym_EigenAllocationCapacity(size_t native_capacity){return eigen_scratch.ptr?eigen_scratch.capacity:native_capacity;}
void CpuNSym_ReleaseEigenScratch(void){
 native_state_release();
 direct_release();
/* SECTOR_PARALLEL_PATCH_BEGIN */
 sector_parallel_release();
/* SECTOR_PARALLEL_PATCH_END */

 if(borrowed_eigen_scratch){B.U=NULL;borrowed_eigen_scratch=0;}
 CpuNSocScratch_Release(&eigen_scratch);
}
int CpuNSym_Prepare(int step){
 int enabled=CpuNSym_Enabled(),control=raw_reset_at();
 if(!enabled&&!control)return 0;
 if(prepared&&last_step==step)return warm;
 if(step<1||step<last_step||(enabled&&!CpuNSym_FullProjected()&&DFTSCF_loop<3))die("requires increasing SCF steps and at least3 allowed symmetry iterations");
 if(step==1){
  if(prepared)die("requires a fresh process for each geometry");
  warm=enabled;force_full=0;full_steps=0;first_raw_step=0;
  if(control&&rank0()){
   int effective=control>=(Pulay_SCF/2+1)?control:Pulay_SCF/2+1;
   printf("CPU_NSOC_RAW_RESET_CONTROL requested=%d effective=%d\n",control,effective);fflush(stdout);
  }
 }
 if(CpuNSym_FullProjected()){
  /* This explicit experimental mode solves the projected operator through
   * native convergence. Residual-based warm release is intentionally disabled. */
  if(step==1&&rank0()){
   printf("CPU_NSOC_FULL_PROJECT representation=PROJECTED release_policy=disabled_full_projected raw_restoration=disabled\n");fflush(stdout);
  }
  prepared=1;last_step=step;warm=1;force_full=0;raw_rebuild=0;full_steps=0;return 1;
 }
 if(control){
  prepared=1;last_step=step;warm=0;raw_rebuild=0;full_steps++;return 0;
 }
 double limit=option("OPENMX_CPU_NSOC_RELEASE_NORM",1e-3),cap=option("OPENMX_CPU_NSOC_WARM_STEPS",0);
 if(limit<0||cap<0||cap>INT_MAX||cap!=(int)cap)die("invalid release norm/warm-step limit");
 if(step>=DFTSCF_loop-1||(step>1&&limit>0&&NormRD[0]<limit)||(cap>0&&step>cap))force_full=1;
 raw_rebuild=warm&&force_full;warm=!force_full;prepared=1;last_step=step;
 if(warm)full_steps=0;else {full_steps++;if(!first_raw_step)first_raw_step=step;}
 if(raw_rebuild&&rank0()){printf("CPU_NSOC_SYM_RESTORE step=%d rebuild_raw_S=1\n",step);fflush(stdout);}
 if(!warm&&rank0()){printf("CPU_NSOC_SYM_RAW step=%d full_steps=%d\n",step,full_steps);fflush(stdout);}
 return warm;
}
int CpuNSym_RawSRebuild(void){return CpuNSym_Enabled()&&raw_rebuild;}
int CpuNSym_ConvergenceBlocked(void){return CpuNSym_Enabled()&&!CpuNSym_FullProjected()&&(!prepared||warm||full_steps<2);}
void CpuNSym_RequestFull(void){if(CpuNSym_Enabled()&&!CpuNSym_FullProjected())force_full=1;}
int CpuNSym_FullSteps(void){return CpuNSym_Enabled()?full_steps:0;}
int CpuNSym_MixingAge(void){
 int enabled=CpuNSym_Enabled(),start=raw_reset_at();
 if(enabled){if(warm||option("OPENMX_CPU_NSOC_RESET_DIIS",0)==0)return 0;start=first_raw_step;}
 if(!start||!prepared)return 0;
 if(start<Pulay_SCF/2+1)start=Pulay_SCF/2+1;
 return last_step>=start?last_step-start+1:0;
}
static int global(int i,int nb,int coord,int np){return(i/nb)*nb*np+coord*nb+i%nb;}
static int local(int i,int nb,int np){return(i/(nb*np))*nb+i%nb;}
static int sector(int i){for(int k=0;k<B.nsec;k++)if(i<B.sec[k].ofs+B.sec[k].n)return (B.version==4&&k==1)?0:k;die("sector index");return 0;}
static int d3_basis(void){return (B.version==1||B.version==3)&&B.nsec==4;}
static int c3_basis(void){return B.version==4&&B.nsec==3;}
static int paired_basis(void){return d3_basis()||c3_basis();}
static void geometry(void){
 if(B.nsec==1||B.version>=2)return; /* v2 carries a generic validated sector basis; geometry is checked by its generator. */
 double reciprocal[3][3],det=0,mx=0;
 for(int i=0;i<3;i++){int j=(i+1)%3,k=(i+2)%3;
  for(int q=0;q<3;q++)reciprocal[i][q]=tv[j+1][(q+1)%3+1]*tv[k+1][(q+2)%3+1]-tv[j+1][(q+2)%3+1]*tv[k+1][(q+1)%3+1];
 }
 for(int q=0;q<3;q++)det+=tv[1][q+1]*reciprocal[0][q];
 if(!isfinite(det)||fabs(det)<1e-12)die("singular geometry cell");
 for(int op=0;op<2;op++)for(int a=0;a<B.natoms;a++){
  int b=op?B.c2perm[a]:B.perm[a];if(b<0||b>=B.natoms||WhatSpecies[a+1]!=WhatSpecies[b+1])die("invalid atom permutation");
  double v[3],d[3];
  if(!op){v[0]=-.5*Gxyz[a+1][1]-sqrt(3.)*.5*Gxyz[a+1][2];v[1]=sqrt(3.)*.5*Gxyz[a+1][1]-.5*Gxyz[a+1][2];v[2]=Gxyz[a+1][3];}
  else for(int i=0;i<3;i++){v[i]=0;for(int j=0;j<3;j++)v[i]+=B.rot[3*i+j]*Gxyz[a+1][j+1];}
  for(int i=0;i<3;i++){for(int j=0;j<3;j++)v[i]+=(op?B.c2shift[j]:B.shift[j])*tv[j+1][i+1];d[i]=Gxyz[b+1][i+1]-v[i];}
  double frac[3],best=HUGE_VAL;
  for(int j=0;j<3;j++){frac[j]=0;for(int i=0;i<3;i++)frac[j]+=d[i]*reciprocal[j][i]/det;frac[j]-=nearbyint(frac[j]);}
  for(int u=-1;u<=1;u++)for(int v=-1;v<=1;v++)for(int w=-1;w<=1;w++){
   double shift[3]={frac[0]+u,frac[1]+v,frac[2]+w},norm=0;
   for(int i=0;i<3;i++){double x=0;for(int j=0;j<3;j++)x+=shift[j]*tv[j+1][i+1];norm+=x*x;}
   if(norm<best)best=norm;
  }
  double e=sqrt(best)*BohrR;if(!isfinite(e))die("nonfinite geometry residual");if(e>mx)mx=e;
 }
 if(mx>5e-5)die("geometry differs from verified D3 mapping");
}
static void basis_init(int n){
 if(B.loaded){if(B.N!=n)die("basis changed");geometry();return;}
 const char*path=getenv("OPENMX_CPU_NSOC_BASIS");if(!path)die("missing OPENMX_CPU_NSOC_BASIS");FILE*f=fopen(path,"rb");if(!f)die("cannot open real basis");
 int32_t h[9];int64_t nz;rd(h,4,9,f);
 if(h[0]!=0x4e534433||(h[1]!=1&&h[1]!=2&&h[1]!=3&&h[1]!=4)||h[2]!=atomnum||h[3]!=n||h[4]<1||h[4]>NSYM_MAX_SECTORS)die("unsupported NSYM basis header");
 B.version=h[1];B.natoms=atomnum;B.N=n;B.nsec=h[4];int sum=0;
 for(int k=0;k<B.nsec;k++){
  int32_t width=(k<4)?h[5+k]:0;
  if(k>=4)rd(&width,4,1,f);
  B.sz[k]=width;if(B.sz[k]<1)die("empty symmetry sector unsupported");
  if(sum>n-B.sz[k])die("sector size overflow");sum+=B.sz[k];
 }
 if(sum!=n)die("sector sizes do not span basis");
 if((B.version==3&&(B.nsec!=4||B.sz[0]!=B.sz[1])) ||
    (B.version==4&&(B.nsec!=3||B.sz[0]!=B.sz[1])))die("paired basis requires aligned equal-sized partners");
 rd(&nz,8,1,f);
 if(nz<1||nz>256LL*n||nz>INT_MAX)die("invalid NSYM sparse basis size");
 B.nnz=(int)nz;
 B.perm=mem(atomnum,4);B.c2perm=mem(atomnum,4);rd(B.perm,4,atomnum,f);rd(B.shift,8,3,f);rd(B.c2perm,4,atomnum,f);rd(B.rot,8,9,f);rd(B.c2shift,8,3,f);
 B.cp=mem(n+1,8);B.ri=mem(nz,4);B.cv=mem(nz,8);rd(B.cp,8,n+1,f);rd(B.ri,4,nz,f);rd(B.cv,8,nz,f);if(fgetc(f)!=EOF)die("trailing real basis data");fclose(f);
 if(B.cp[0]!=0||B.cp[n]!=nz)die("invalid CSC pointers");
 int*rp=mem(n+1,4),*ci=mem(nz,4);double*rv=mem(nz,8);
 for(int j=0;j<n;j++){
  if(B.cp[j]>B.cp[j+1]||B.cp[j+1]>nz)die("CSC pointer range");
  for(int64_t z=B.cp[j];z<B.cp[j+1];z++){int i=B.ri[z];if(i<0||i>=n||!isfinite(B.cv[z]))die("invalid sparse basis value");rp[i+1]++;}
 }
 for(int i=0;i<n;i++)rp[i+1]+=rp[i];int*pos=mem(n,4);memcpy(pos,rp,n*4);
 for(int j=0;j<n;j++)for(int64_t z=B.cp[j];z<B.cp[j+1];z++){int t=pos[B.ri[z]]++;ci[t]=j;rv[t]=B.cv[z];}free(pos);
 for(int j=0;j<n;j++){
  int ids[256],count=0;double val[256];
  for(int64_t z=B.cp[j];z<B.cp[j+1];z++)for(int a=rp[B.ri[z]];a<rp[B.ri[z]+1];a++){
   int c=ci[a],q=0;while(q<count&&ids[q]!=c)q++;if(q==count){if(count==256)die("basis overlap degree too large");ids[count]=c;val[count++]=0;}val[q]+=B.cv[z]*rv[a];
  }
  for(int q=0;q<count;q++)if(fabs(val[q]-(ids[q]==j))>1e-10)die("basis is not orthogonal");
 }
 free(rp);free(ci);free(rv);
 for(int i=0;i<3;i++)for(int j=0;j<3;j++){double v=0;for(int k=0;k<3;k++)v+=B.rot[3*k+i]*B.rot[3*k+j];if(!isfinite(v)||fabs(v-(i==j))>1e-10)die("invalid C2 rotation");}
 for(int i=0;i<3;i++)if(!isfinite(B.shift[i])||!isfinite(B.c2shift[i]))die("invalid symmetry translation");
 geometry();B.cap=(size_t)na_rows_max*na_cols_max;B.active=(size_t)na_rows*na_cols;
 if(compact_eigen_scratch()){
  B.U=CpuNSocScratch_Borrow(&eigen_scratch,B.cap);
  if(!B.U)die("EIGEN_SCRATCH missing or insufficient caller-owned EVec1 capacity");
  borrowed_eigen_scratch=1;
  int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
  printf("CPU_NSOC_EIGEN_SCRATCH_ACTIVE rank=%d reused_doubles=%zu removed_separate_bytes=%zu extra_distributed_matrices=0 ownership=borrowed_EVec1\n",rank,B.cap,B.cap*sizeof(double));fflush(stdout);
 }else B.U=mem(B.cap,8);
 B.scap=0;B.hcap=0;int ofs=0,zero=0,info;
 for(int k=0;k<B.nsec;k++){
  Sector*t=B.sec+k;t->n=B.sz[k];t->ofs=ofs;ofs+=t->n;t->nb=nblk;
  while(t->nb>1&&(t->nb*np_rows>t->n||t->nb*np_cols>t->n))t->nb/=2;
  if(t->n<np_rows||t->n<np_cols)die("too many processes for small symmetry sector");
  t->nr=numroc_(&t->n,&t->nb,&my_prow,&zero,&np_rows);t->nc=numroc_(&t->n,&t->nb,&my_pcol,&zero,&np_cols);int ld=t->nr?t->nr:1;
  descinit_(t->desc,&t->n,&t->n,&t->nb,&t->nb,&zero,&zero,&ictxt1,&ld,&info);if(info)die("sector descriptor initialization");
  t->soff=(paired_basis()&&k==1)?0:B.scap;if(!(paired_basis()&&k==1))B.scap+=(size_t)t->nr*t->nc*(c3_basis()&&k==0?2:1);
  B.hoff[k]=B.hcap;B.hcap+=(size_t)t->nr*t->nc;
  if((size_t)t->nr*t->nc>B.cap)die("sector workspace exceeds stock allocation");
 }
 if(B.scap>B.active)die("overlap factor cache exceeds actual local stock S allocation");
 if(CpuNSym_NativeReturnEnabled()&&B.hcap>B.active)die("packed sector H exceeds actual local stock H allocation");
 int p;MPI_Comm_size(mpi_comm_level1,&p);if(np_rows*np_cols!=p)die("requires one full real BLACS world");
 int me[2]={my_prow,my_pcol},*all=mem(2*p,4);MPI_Allgather(me,2,MPI_INT,all,2,MPI_INT,mpi_comm_level1);B.rankmap=mem(p,4);for(int k=0;k<p;k++)B.rankmap[k]=-1;
 for(int k=0;k<p;k++){int q=all[2*k]*np_cols+all[2*k+1];if(q<0||q>=p||B.rankmap[q]!=-1)die("invalid BLACS-to-MPI mapping");B.rankmap[q]=k;}free(all);
 B.loaded=1;
 if(rank0()){printf("CPU_NSOC_SYM_BASIS n=%d sector_count=%d sectors=%d,%d,%d,%d nnz=%d extra_distributed_matrices=%d local_workspace_bytes=%zu version=%d paired_reuse=%s\n",n,B.nsec,B.sz[0],B.sz[1],B.sz[2],B.sz[3],B.nnz,borrowed_eigen_scratch?0:1,B.cap*8,B.version,d3_basis()?"D3_real":(c3_basis()?"C3_complex":"off"));fflush(stdout);}
}
static void buildU(void){
 memset(B.U,0,B.cap*8);
 for(int j=0;j<B.N;j++)if((j/nblk)%np_cols==my_pcol){int jl=local(j,nblk,np_cols);for(int64_t z=B.cp[j];z<B.cp[j+1];z++){int i=B.ri[z];if((i/nblk)%np_rows==my_prow)B.U[(size_t)jl*na_rows+local(i,nblk,np_rows)]=B.cv[z];}}
}
static void copyblock(double*full,int*desc,Sector*t,double*small){int one=1,a=t->ofs+1;pdgemr2d_(&t->n,&t->n,full,&a,&a,desc,small,&one,&one,t->desc,&ictxt1);}
static void putblock(double*small,Sector*t,double*full,int*desc,int ofs){int one=1,a=ofs+1;pdgemr2d_(&t->n,&t->n,small,&one,&one,t->desc,full,&a,&a,desc,&ictxt1);}
static void putblock_rect(double*small,Sector*t,int cols,double*full,int*desc,int ofs){int one=1,a=ofs+1;pdgemr2d_(&t->n,&cols,small,&one,&one,t->desc,full,&a,&a,desc,&ictxt1);}
static CpuNSparsePlan *sparse_plan_get(double*plan_time){
 static CpuNSparsePlan *sparse_plan=NULL;
 if(plan_time)*plan_time=0;
 if(!sparse_plan){
  char error[256];double start=MPI_Wtime();
  sparse_plan=CpuNSparse_Create(B.N,nblk,na_rows,na_cols,np_rows,np_cols,my_prow,my_pcol,B.cp,B.ri,B.cv,mpi_comm_level1,error,sizeof(error));
  if(!sparse_plan)die(error);double elapsed=MPI_Wtime()-start;
  if(plan_time)*plan_time=elapsed;
  int rank;MPI_Comm_rank(mpi_comm_level1,&rank);
  printf("CPU_NSOC_SPARSE_PLAN rank=%d n=%d components=%d maxblock=%d panel_bytes=%zu metadata_bytes=%zu local_rows=%d local_cols=%d plan_seconds=%.9f\n",rank,B.N,CpuNSparse_Components(sparse_plan),CpuNSparse_MaxBlock(sparse_plan),CpuNSparse_PanelBytes(sparse_plan),CpuNSparse_MetadataBytes(sparse_plan),na_rows,na_cols,elapsed);fflush(stdout);
 }
 return sparse_plan;
}
static void congruence(double*A,int*descA,double*work,int*descW){
 static unsigned long profile_call=0;
 double backend=option("OPENMX_CPU_NSOC_SPARSE_U",0),profiling=option("OPENMX_CPU_NSOC_SPARSE_PROFILE",0);
 if((backend!=0&&backend!=1)||(profiling!=0&&profiling!=1))die("SPARSE_U and SPARSE_PROFILE must be0 or1");
 int profile=profiling!=0;double begin=profile?MPI_Wtime():0,plan_time=0,right=0,left=0;
 CpuNSparseTimings t={0};
 if(backend!=0){
  if(descA[2]!=B.N||descA[3]!=B.N||descW[2]!=B.N||descW[3]!=B.N||descA[4]!=nblk||descA[5]!=nblk||descW[4]!=nblk||descW[5]!=nblk||descA[6]||descA[7]||descW[6]||descW[7]||descA[8]!=na_rows||descW[8]!=na_rows)die("sparse congruence requires matching stock block-cyclic descriptors");
  char error[256];
  CpuNSparsePlan*sparse_plan=sparse_plan_get(&plan_time);
  if(CpuNSparse_Congruence(sparse_plan,A,descA[8],work,descW[8],profile,&t,error,sizeof(error)))die(error);
 }
 else{
  /* Default backend is the unchanged fp1 pair of full-N dense products. */
  int n=B.N,one=1;double a=1,b=0;buildU();
  double p0=profile?MPI_Wtime():0;
  pdgemm_("N","N",&n,&n,&n,&a,A,&one,&one,descA,B.U,&one,&one,descC,&b,work,&one,&one,descW);
  if(profile){right=MPI_Wtime()-p0;p0=MPI_Wtime();}
  pdgemm_("T","N",&n,&n,&n,&a,B.U,&one,&one,descC,work,&one,&one,descW,&b,A,&one,&one,descA);
  if(profile)left=MPI_Wtime()-p0;
 }
 if(profile){
  double total=MPI_Wtime()-begin;int rank;MPI_Comm_rank(mpi_comm_level1,&rank);profile_call++;
  printf("CPU_NSOC_SPARSE_PHASE rank=%d step=%d matrix=%s call=%lu backend=%s plan=%.9f right_pack=%.9f right_collective=%.9f right_multiply=%.9f left_pack=%.9f left_collective=%.9f left_multiply=%.9f dense_right=%.9f dense_left=%.9f total=%.9f\n",rank,last_step,descA==descS?"S":"H",profile_call,backend!=0?"sparse":"dense",plan_time,t.right_pack,t.right_collective,t.right_multiply,t.left_pack,t.left_collective,t.left_multiply,right,left,total);fflush(stdout);
 }
}
static void backtransform(int step,int n,double*H,double*W,int phase){
 double setting=option("OPENMX_CPU_NSOC_SPARSE_BACK",0);
 if(setting!=0&&setting!=1&&setting!=2)die("SPARSE_BACK must be0,1 or2");
 int mode=(int)setting,one=1,rank;double a=1,b=0;
 double batch_setting=option("OPENMX_CPU_NSOC_BACK_BATCH_ROWS",0);
 if(batch_setting<0||batch_setting>1024||batch_setting!=(int)batch_setting)die("BACK_BATCH_ROWS must be an integer0..1024");
 int batch_rows=(int)batch_setting,batches=0;size_t batch_bytes=0;
 int profile=phase||mode==2||option("OPENMX_CPU_NSOC_SPARSE_PROFILE",0)!=0;
 double start=MPI_Wtime(),t0=start,mask=0,dense=0,plan=0,transpose=0,compare=0,copy=0;
 double maxabs=0,scale=0,relative=0;CpuNSparseTimings timings={0};
 if(mode){
  /* H still contains the transformed Hamiltonian outside recovered diagonal
   * blocks. The dense reference never reads those entries; U*E must not. */
  for(int j=0;j<na_cols;j++){
   int gj=global(j,nblk,my_pcol,np_cols);if(gj>=n)continue;
   for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi<n&&sector(gi)!=sector(gj))H[(size_t)j*na_rows+i]=0;}
  }
  mask=MPI_Wtime()-t0;t0=MPI_Wtime();
 }
 if(mode!=1){
  buildU();memset(W,0,B.cap*sizeof(double));
  for(int k=0;k<B.nsec;k++){if(c3_basis()&&k==1)continue;Sector*t=B.sec+k;int off=t->ofs+1,width=t->n*(c3_basis()&&k==0?2:1);pdgemm_("T","T",&width,&n,&width,&a,H,&off,&off,descH,B.U,&one,&off,descC,&b,W,&off,&one,descC);}
  dense=MPI_Wtime()-t0;
 }
 if(mode){
  char error[256];CpuNSparsePlan*p=sparse_plan_get(&plan);
  /* Sparse U is borrowed from the plan, so B.U can hold U*E after the
   * dense verification product. Existing three distributed buffers suffice. */
  if(CpuNSparse_Left(p,H,na_rows,B.U,na_rows,batch_rows,profile,&timings,error,sizeof(error)))die(error);
  batch_bytes=CpuNSparse_BackPanelBytes(p);batches=CpuNSparse_BackBatches(p);
  t0=MPI_Wtime();pdtran_(&n,&n,&a,B.U,&one,&one,descC,&b,H,&one,&one,descH);transpose=MPI_Wtime()-t0;
  if(mode==2){
   t0=MPI_Wtime();double local_values[2]={0,0},global_values[2];int bad=0,global_bad;
   for(size_t i=0;i<B.active;i++){
    double x=W[i],y=H[i];if(!isfinite(x)||!isfinite(y)){bad=1;continue;}
    double diff=fabs(x-y);if(!isfinite(diff)){bad=1;continue;}
    if(diff>local_values[0])local_values[0]=diff;if(fabs(x)>local_values[1])local_values[1]=fabs(x);
   }
   MPI_Allreduce(&bad,&global_bad,1,MPI_INT,MPI_MAX,mpi_comm_level1);
   MPI_Allreduce(local_values,global_values,2,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);
   maxabs=global_values[0];scale=global_values[1];relative=maxabs/fmax(scale,1e-300);
   MPI_Comm_rank(mpi_comm_level1,&rank);
   if(rank==0){printf("CPU_NSOC_SPARSE_BACK_VERIFY step=%d max_abs=%.17g scale=%.17g relative=%.17g nonfinite=%d\n",step,maxabs,scale,relative,global_bad);fflush(stdout);}
   if(global_bad||!isfinite(relative)||relative>1e-10)die("sparse eigenvector backtransform differs from dense reference");
   compare=MPI_Wtime()-t0;
  }
  t0=MPI_Wtime();memcpy(W,H,B.active*sizeof(double));copy=MPI_Wtime()-t0;
 }
 if(profile){
  double total=MPI_Wtime()-start;MPI_Comm_rank(mpi_comm_level1,&rank);
  printf("CPU_NSOC_SPARSE_BACK_PHASE rank=%d step=%d mode=%d batch_rows=%d batches=%d back_panel_bytes=%zu mask=%.9f dense=%.9f plan=%.9f pack=%.9f collective=%.9f multiply=%.9f transpose=%.9f compare=%.9f copy=%.9f total=%.9f\n",rank,step,mode,batch_rows,batches,batch_bytes,mask,dense,plan,timings.left_pack,timings.left_collective,timings.left_multiply,transpose,compare,copy,total);fflush(stdout);
 }
}
static double offblock_norm(double*A){
 double a[2]={0,0},s[2];for(int j=0;j<na_cols;j++){int gj=global(j,nblk,my_pcol,np_cols);if(gj>=B.N)continue;for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi>=B.N)continue;double v=A[(size_t)j*na_rows+i];a[1]+=v*v;if(sector(gi)!=sector(gj))a[0]+=v*v;}}
 MPI_Allreduce(a,s,2,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);return sqrt(s[0]/fmax(s[1],1e-300));
}
static double paired_relative(const double*x,const double*y,size_t count,double sign){
 double v[2]={0,0},g[2];for(size_t i=0;i<count;i++){double d=x[i]-sign*y[i];v[0]+=d*d;v[1]+=x[i]*x[i]+y[i]*y[i];}
 MPI_Allreduce(v,g,2,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 return sqrt(g[0]/fmax(g[1],1e-300));
}
static void paired_verify(double defect,const char*what){
 if(!isfinite(defect)||defect>1e-10)die(what);
}
/* Convert a canonical real E block [A,-B;B,A] into A+iB.
 * Two existing distributed real buffers hold the complex problem. No NxN
 * replicated arrays and no additional distributed matrix allocation. */
static void c3_complex_extract(double*full,int*desc,double*packed,double*scratch,int step,const char*field){
 Sector*t=B.sec;size_t z=(size_t)t->nr*t->nc;int one=1,x=t->ofs+1,y=B.sec[1].ofs+1;
 if(2*z>B.cap)die("C3 complex buffer exceeds existing workspace");
 pdgemr2d_(&t->n,&t->n,full,&x,&x,desc,packed,&one,&one,t->desc,&ictxt1);
 pdgemr2d_(&t->n,&t->n,full,&y,&y,desc,scratch,&one,&one,t->desc,&ictxt1);
 double local_defect[3]={0.,0.,0.},global_defect[3];
 for(size_t i=0;i<z;i++){double d=packed[i]-scratch[i];local_defect[0]+=d*d;local_defect[2]+=packed[i]*packed[i]+scratch[i]*scratch[i];packed[i]=.5*(packed[i]+scratch[i]);}
 pdgemr2d_(&t->n,&t->n,full,&y,&x,desc,scratch,&one,&one,t->desc,&ictxt1);
 pdgemr2d_(&t->n,&t->n,full,&x,&y,desc,packed+z,&one,&one,t->desc,&ictxt1);
 for(size_t i=0;i<z;i++){double d=scratch[i]+packed[z+i];local_defect[1]+=d*d;local_defect[2]+=scratch[i]*scratch[i]+packed[z+i]*packed[z+i];}
 MPI_Allreduce(local_defect,global_defect,3,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 double diagonal=sqrt(global_defect[0]/fmax(global_defect[2],1e-300));
 double cross=sqrt(global_defect[1]/fmax(global_defect[2],1e-300));
 paired_verify(diagonal,"C3 E diagonal partners disagree");
 paired_verify(cross,"C3 E cross blocks violate complex structure");
 for(size_t i=0;i<z;i++)scratch[i]=.5*(scratch[i]-packed[z+i]);
 for(size_t i=z;i-->0;){packed[2*i]=packed[i];packed[2*i+1]=scratch[i];}
 if(rank0()){printf("CPU_NSOC_C3_COMPLEX_COV step=%d field=%s diagonal_relative=%.12e cross_relative=%.12e n=%d\n",step,field,diagonal,cross,t->n);fflush(stdout);}
}
static void overlap_factors(double*S,double*work){
 double defect;int packed_direct=direct_active&&direct_mode==1&&d3_basis()&&CpuNSym_NativeReturnEnabled();
 if(direct_active){
  if(direct_mode==2){congruence(S,descS,work,descC);direct_reference_s_defect=offblock_norm(S);}else direct_reference_s_defect=0.;
  direct_assemble(direct_step,2,S,descS,work);defect=direct_reference_s_defect;
 }else{congruence(S,descS,work,descC);defect=offblock_norm(S);}
 double strict=option("OPENMX_CPU_NSOC_GENERIC_COV_STRICT",1);double cov_tol=option("OPENMX_CPU_NSOC_GENERIC_COV_TOL",1e-5);
 if(!isfinite(cov_tol)||cov_tol<=0)die("GENERIC_COV_TOL must be positive");
 if(!isfinite(defect)||(d3_basis()&&defect>1e-5)||(B.version>=2&&strict!=0&&defect>cov_tol))die("raw overlap fails sector covariance; check the validated basis and operation convention");
 if(packed_direct){
  /* Rpa produced one E representative and both distinct A blocks. Their
   * 32-block outputs were redistributed directly into the same packed
   * factor descriptors; no full AO S or duplicate E was ever stored. */
  for(int k=0;k<B.nsec;k++)if(k!=1){Sector*t=B.sec+k;int one=1,info;pdpotrf_("L",&t->n,S+t->soff,&one,&one,t->desc,&info);if(info)die("direct packed overlap Cholesky failed");}
  if(rank0()){printf("CPU_NSOC_DIRECT_OVERLAP packed_factors=1 full_AO_S_intermediate=0 source_blocksize=32 target_sector_descriptors=unchanged equivalent_E=one_certified_representative\n");fflush(stdout);}
 }else{
 if(c3_basis()){
  c3_complex_extract(S,descS,work,B.U,1,"S");Sector*t=B.sec;int one=1,info;
  pzpotrf_("L",&t->n,(dcomplex*)work,&one,&one,t->desc,&info);if(info)die("C3 complex overlap Cholesky failed");
  memcpy(B.U,work,2*(size_t)t->nr*t->nc*sizeof(double));
 }
 for(int k=0;k<B.nsec;k++){
  if((paired_basis()&&k==1)||(c3_basis()&&k==0))continue;
  Sector*t=B.sec+k;double*p=B.U+t->soff;copyblock(S,descS,t,p);
  if(d3_basis()&&k==0){copyblock(S,descS,B.sec+1,work);if(B.version==3)paired_verify(paired_relative(p,work,(size_t)t->nr*t->nc,1.),"D3 overlap equivalent blocks disagree");for(size_t i=0;i<(size_t)t->nr*t->nc;i++)p[i]=.5*(p[i]+work[i]);}
  int one=1,info;pdpotrf_("L",&t->n,p,&one,&one,t->desc,&info);if(info)die("projected overlap Cholesky failed");
 }
 memcpy(S,B.U,B.scap*8);
 }
 size_t reclaimed=CpuReclaim_Pages(S,B.scap*8,B.active*8);
 {int rank;MPI_Comm_rank(mpi_comm_level1,&rank);printf("CPU_NSOC_OVERLAP_RECLAIM rank=%d allocation_bytes=%zu live_factor_bytes=%zu returned_bytes=%zu virtual_capacity_retained=1\n",rank,B.active*8,B.scap*8,reclaimed);fflush(stdout);}
 if(rank0()){
  if(direct_active&&direct_mode==1)printf("CPU_NSOC_SYM_OVERLAP covariance=construction_certificate_actual_D3_AO_intertwiner_distinct_A_fresh_Reynolds raw_offblock_measurement=not_performed cached_factor_doubles=%zu\n",B.scap);
  else printf("CPU_NSOC_SYM_OVERLAP projected_offblock_relative=%.12e cached_factor_doubles=%zu\n",defect,B.scap);
  fflush(stdout);
 }
}
static void dump_small(double*H,double*S,double*e,int maxn,const char*file,int step){
 if(B.N>128||option("OPENMX_CPU_NSOC_DUMP",0)==0)return;
 int n=B.N,ne=e?n:0,fields=e?1:2;size_t nn=(size_t)n*n;double*a=mem(fields*nn,8),*sum=mem(fields*nn,8);
 for(int j=0;j<na_cols;j++){int gj=global(j,nblk,my_pcol,np_cols);if(gj>=n)continue;for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi>=n)continue;size_t at=(size_t)j*na_rows+i;
  if(e)a[(size_t)gj*n+gi]=H[at];else{a[(size_t)gi*n+gj]=H[at];a[nn+(size_t)gi*n+gj]=S[at];}
 }}
 MPI_Reduce(a,sum,(int)(fields*nn),MPI_DOUBLE,MPI_SUM,0,mpi_comm_level1);
 if(rank0()){FILE*f=fopen(file,"wb");if(!f)die("cannot write diagnostic dump");int hd[4]={n,e?maxn:step,e?step:0,nblk};fwrite(hd,4,4,f);if(e){if(maxn!=n)die("small eigen dump requires first full returned spectrum");fwrite(e+1,8,ne,f);}fwrite(sum,8,fields*nn,f);fclose(f);}
 free(a);free(sum);
}
static int eigcmp(const void*a,const void*b){const Eigen*x=a,*y=b;if(x->value<y->value)return-1;if(x->value>y->value)return 1;return(x->source>y->source)-(x->source<y->source);}
static void reorder(double*ev,double*ko,int maxn,double*from,double*to){
 int p;MPI_Comm_size(mpi_comm_level1,&p);Eigen*sort=mem(B.N,sizeof(Eigen));int*inv=mem(B.N,4);
 for(int i=0;i<B.N;i++){if(!isfinite(ev[i]))die("nonfinite projected eigenvalue");sort[i].value=ev[i];sort[i].source=i;}qsort(sort,B.N,sizeof(Eigen),eigcmp);
 for(int i=0;i<B.N;i++){ko[i+1]=sort[i].value;inv[sort[i].source]=i;}
 memset(to,0,B.active*8);int*ct=mem(p,4),*rc=mem(p,4),*sd=mem(p,4),*rdp=mem(p,4),*pos=mem(p,4);
 int tile_size=(int)option("OPENMX_CPU_NSOC_REORDER_TILE",1024);if(tile_size<1||tile_size>B.N)tile_size=B.N;
 for(int tile=0;tile<B.N;tile+=tile_size){int end=tile+tile_size;if(end>B.N)end=B.N;memset(ct,0,p*4);
  for(int j=0;j<na_cols;j++){int gj=global(j,nblk,my_pcol,np_cols);if(gj<tile||gj>=end)continue;for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi>=B.N||inv[gi]>=maxn)continue;int dest=B.rankmap[((inv[gi]/nblk)%np_rows)*np_cols+(gj/nblk)%np_cols];ct[dest]++;}}
  MPI_Alltoall(ct,1,MPI_INT,rc,1,MPI_INT,mpi_comm_level1);sd[0]=rdp[0]=0;for(int q=1;q<p;q++){sd[q]=sd[q-1]+ct[q-1];rdp[q]=rdp[q-1]+rc[q-1];}int ns=sd[p-1]+ct[p-1],nr=rdp[p-1]+rc[p-1];Entry*send=mem(ns,sizeof(Entry)),*recv=mem(nr,sizeof(Entry));memcpy(pos,sd,p*4);
  for(int j=0;j<na_cols;j++){int gj=global(j,nblk,my_pcol,np_cols);if(gj<tile||gj>=end)continue;for(int i=0;i<na_rows;i++){int gi=global(i,nblk,my_prow,np_rows);if(gi>=B.N||inv[gi]>=maxn)continue;int dest=B.rankmap[((inv[gi]/nblk)%np_rows)*np_cols+(gj/nblk)%np_cols];Entry*v=send+pos[dest]++;v->row=inv[gi];v->col=gj;v->value=from[(size_t)j*na_rows+i];}}
  for(int q=0;q<p;q++){if((int64_t)(ct[q]+rc[q]+sd[q]+rdp[q])*sizeof(Entry)>INT_MAX)die("coefficient MPI count overflow");ct[q]*=sizeof(Entry);rc[q]*=sizeof(Entry);sd[q]*=sizeof(Entry);rdp[q]*=sizeof(Entry);}
  MPI_Alltoallv(send,ct,sd,MPI_BYTE,recv,rc,rdp,MPI_BYTE,mpi_comm_level1);
  for(int q=0;q<nr;q++){Entry*v=recv+q;if(v->row<0||v->row>=maxn||v->col<tile||v->col>=end||(v->row/nblk)%np_rows!=my_prow||(v->col/nblk)%np_cols!=my_pcol||!isfinite(v->value))die("coefficient redistribution identity");to[(size_t)local(v->col,nblk,np_cols)*na_rows+local(v->row,nblk,np_rows)]=v->value;}
  free(send);free(recv);
 }
 free(sort);free(inv);free(ct);free(rc);free(sd);free(rdp);free(pos);
}
#include "fused_return_native_v2.inc"
#include "cpu_native_return.inc"
/* CHECK uses only sampled O(N) vectors. Raw matrix blocks remain distributed
 * over owning atoms; the six local group contributions are summed once. */
static void group_vector(const double*v,double*out,double*coeff,int power,int flip,int transpose){
 int n=B.N,m=B.sz[0];
 for(int j=0;j<n;j++){double z=0;for(int64_t a=B.cp[j];a<B.cp[j+1];a++)z+=B.cv[a]*v[B.ri[a]];coeff[j]=z;}
 double c=power?-.5:1,s=power==1?sqrt(3.)*.5:(power==2?-sqrt(3.)*.5:0);
 for(int j=0;j<m;j++){
  double x=coeff[j],y=coeff[m+j],a,b;
  if(transpose){if(flip)y=-y;a=c*x+s*y;b=-s*x+c*y;}
  else{a=c*x-s*y;b=s*x+c*y;if(flip)b=-b;}
  coeff[j]=a;coeff[m+j]=b;
 }
 if(flip)for(int j=B.sec[3].ofs;j<n;j++)coeff[j]=-coeff[j];
 memset(out,0,(size_t)n*8);
 for(int j=0;j<n;j++)for(int64_t a=B.cp[j];a<B.cp[j+1];a++)out[B.ri[a]]+=B.cv[a]*coeff[j];
}
static void raw_block_action(double*****nh,double****olp,const int*offset,const double*v,double*h,double*s){
 memset(h,0,(size_t)B.N*8);memset(s,0,(size_t)B.N*8);
 for(int ma=1;ma<=Matomnum;ma++){
  int ga=M2G[ma],nr=Spe_Total_CNO[WhatSpecies[ga]],row=offset[ga];
  for(int b=0;b<=FNAN[ga];b++){
   int gb=natn[ga][b],nc=Spe_Total_CNO[WhatSpecies[gb]],col=offset[gb];
   for(int i=0;i<nr;i++)for(int j=0;j<nc;j++){
    double x=v[col+j];h[row+i]+=nh[0][ma][b][i][j]*x;s[row+i]+=olp[ma][b][i][j]*x;
   }
  }
 }
}
static void sampled_verify(int step,int maxn,double*ko,double*eval,double*unsorted,double*sorted,double*****nh,double****olp){
 int n=B.N,checked=maxn,probe=0;
 double*probe_ko=NULL;
 if(maxn==n){
  checked=n-1;probe=1;probe_ko=mem(n+1,8);
  /* The same partial-state row redistribution, in an existing workspace,
   * leaves live Hs/ko untouched while tiny fixtures exercise MaxN<N. */
  reorder(eval,probe_ko,checked,unsorted,B.U);sorted=B.U;
 }
 int bands[3]={0,checked/2,checked-1};
 double*C=mem((size_t)3*n,8),*HP=mem((size_t)3*n,8),*SP=mem((size_t)3*n,8);
 for(int k=0;k<3;k++)if((bands[k]/nblk)%np_rows==my_prow){
  int i=local(bands[k],nblk,np_rows);
  for(int j=0;j<na_cols;j++){int gj=global(j,nblk,my_pcol,np_cols);if(gj<n)C[(size_t)k*n+gj]=sorted[(size_t)j*na_rows+i];}
 }
 MPI_Allreduce(MPI_IN_PLACE,C,3*n,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 int*offset=mem(atomnum+1,4);int end=0;
 for(int ga=1;ga<=atomnum;ga++){offset[ga]=end;end+=Spe_Total_CNO[WhatSpecies[ga]];}if(end!=n)die("verification basis ordering");
 double bound[2]={0,0},bounds[2];
 for(int ma=1;ma<=Matomnum;ma++){
  int ga=M2G[ma],nr=Spe_Total_CNO[WhatSpecies[ga]];
  for(int i=0;i<nr;i++){
   double h=0,s=0;
   for(int b=0;b<=FNAN[ga];b++){int gb=natn[ga][b],nc=Spe_Total_CNO[WhatSpecies[gb]];for(int j=0;j<nc;j++){h+=fabs(nh[0][ma][b][i][j]);s+=fabs(olp[ma][b][i][j]);}}
   if(!isfinite(h)||!isfinite(s))die("nonfinite raw verification bound");
   if(h>bound[0])bound[0]=h;if(s>bound[1])bound[1]=s;
  }
 }
 MPI_Allreduce(bound,bounds,2,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);
 double*q=mem(n,8),*coeff=mem(n,8),*h=mem(n,8),*s=mem(n,8),*y=mem(n,8);
 if(B.version>=2){
  /* Generic v2 has already passed the matrix-level covariance gate.  Verify
   * sampled eigenvectors directly against the raw distributed H/S action;
   * the D3 branch below additionally averages the six legacy operations. */
  for(int k=0;k<3;k++){
   raw_block_action(nh,olp,offset,C+(size_t)k*n,h,s);
   memcpy(HP+(size_t)k*n,h,(size_t)n*8);
   memcpy(SP+(size_t)k*n,s,(size_t)n*8);
  }
 }else{
  for(int k=0;k<3;k++)for(int flip=0;flip<2;flip++)for(int power=0;power<3;power++){
   group_vector(C+(size_t)k*n,q,coeff,power,flip,0);
   raw_block_action(nh,olp,offset,q,h,s);
   group_vector(h,y,coeff,power,flip,1);for(int i=0;i<n;i++)HP[(size_t)k*n+i]+=y[i]/6;
   group_vector(s,y,coeff,power,flip,1);for(int i=0;i<n;i++)SP[(size_t)k*n+i]+=y[i]/6;
  }
 }
 MPI_Allreduce(MPI_IN_PLACE,HP,3*n,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 MPI_Allreduce(MPI_IN_PLACE,SP,3*n,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 for(int i=0;i<3*n;i++)if(!isfinite(C[i])||!isfinite(HP[i])||!isfinite(SP[i]))die("nonfinite sampled verification vector");
 double residual=0,scaled=0,orth=0;
 for(int k=0;k<3;k++){
  double rr=0,hh=0,ss=0,rinf=0,cinf=0,e=ko[bands[k]+1];
  for(int i=0;i<n;i++){
   double v=HP[(size_t)k*n+i]-e*SP[(size_t)k*n+i];rr+=v*v;hh+=HP[(size_t)k*n+i]*HP[(size_t)k*n+i];ss+=SP[(size_t)k*n+i]*SP[(size_t)k*n+i];
   if(fabs(v)>rinf)rinf=fabs(v);if(fabs(C[(size_t)k*n+i])>cinf)cinf=fabs(C[(size_t)k*n+i]);
  }
  if(!isfinite(e)||!isfinite(rr)||!isfinite(hh)||!isfinite(ss))die("nonfinite sampled verification norm");
  double r=sqrt(rr)/fmax(sqrt(hh)+fabs(e)*sqrt(ss),1e-300);
  double z=rinf/fmax((bounds[0]+fabs(e)*bounds[1])*cinf,1e-300);if(!isfinite(r)||!isfinite(z))die("nonfinite sampled verification residual");if(r>residual)residual=r;if(z>scaled)scaled=z;
  for(int l=0;l<3;l++){double v=0;for(int i=0;i<n;i++)v+=C[(size_t)k*n+i]*SP[(size_t)l*n+i];v=fabs(v-(k==l));if(!isfinite(v))die("nonfinite sampled verification Sorth");if(v>orth)orth=v;}
 }
 if(rank0()){
  printf("CPU_NSOC_SYM_VERIFY step=%d native_maxn=%d checked_maxn=%d bands=%d,%d,%d partial_reorder_probe=%d action_residual=%.12e scaled_raw_bound_residual=%.12e Sorth=%.12e raw_Hinf_bound=%.12e raw_Sinf_bound=%.12e extra_vector_doubles=%d\n",step,maxn,checked,bands[0]+1,bands[1]+1,bands[2]+1,probe,residual,scaled,orth,bounds[0],bounds[1],14*n+(probe?n+1:0));fflush(stdout);
  if(n<=128&&option("OPENMX_CPU_NSOC_DUMP",0)!=0){
   char name[80];snprintf(name,sizeof(name),"nsym_verify_step%d.bin",step);FILE*f=fopen(name,"wb");if(!f)die("verification dump open");
   int32_t hd[8]={n,maxn,checked,step,bands[0]+1,bands[1]+1,bands[2]+1,probe};double e[3]={ko[bands[0]+1],ko[bands[1]+1],ko[bands[2]+1]};
   fwrite(hd,4,8,f);fwrite(e,8,3,f);fwrite(C,8,3*n,f);fwrite(HP,8,3*n,f);fwrite(SP,8,3*n,f);fclose(f);
  }
 }
 double verify_tol=B.version>=2?option("OPENMX_CPU_NSOC_GENERIC_COV_TOL",1e-5):1e-8;
 if(!isfinite(verify_tol)||verify_tol<=0)die("GENERIC_COV_TOL must be positive");
 if(!isfinite(residual)||!isfinite(orth)||residual>verify_tol||orth>verify_tol)die("sampled projected generalized eigenpair verification failed");
 free(probe_ko);free(C);free(HP);free(SP);free(offset);free(q);free(coeff);free(h);free(s);free(y);
}

static void CpuNSym_Solve_legacy(int step,int n,int maxn,double*ko,double*S,double*H,double*W,double*****nh,double****olp){
 if(!CpuNSym_Prepare(step)||maxn<1||maxn>n)die("invalid projected solver entry");double start=MPI_Wtime();basis_init(n);int phase=option("OPENMX_CPU_NSOC_PHASE",0)!=0;double partial_opt=option("OPENMX_CPU_NSOC_PARTIAL_EVP",0);if(partial_opt!=0&&partial_opt!=1)die("PARTIAL_EVP must be0 or1");int partial=partial_opt!=0;
 if(rank0()){printf("CPU_NSOC_SYM_PROJECT step=%d n=%d sectors=%d,%d,%d,%d maxn=%d evp_mode=%s\n",step,n,B.sz[0],B.sz[1],B.sz[2],B.sz[3],maxn,partial?"partial":"full");fflush(stdout);}
 if(step==1)dump_small(H,S,NULL,maxn,"nsym_raw_HS.bin",step);
 if(step==1)overlap_factors(S,W);double tS=MPI_Wtime();congruence(H,descH,W,descC);double tH=MPI_Wtime();
 double*eval=mem(n,8);int one=1;double a=1,b=0;MPI_Comm rowcomm,colcomm;
 MPI_Comm_split(mpi_comm_level1,my_pcol,my_prow,&rowcomm);MPI_Comm_split(mpi_comm_level1,my_prow,my_pcol,&colcomm);int rows=MPI_Comm_c2f(rowcomm),cols=MPI_Comm_c2f(colcomm),world=MPI_Comm_c2f(mpi_comm_level1);
 for(int k=0;k<B.nsec;k++){
  Sector*t=B.sec+k;double t0=MPI_Wtime();copyblock(H,descH,t,W);
  if(d3_basis()&&k==0){copyblock(H,descH,B.sec+1,B.U);for(size_t i=0;i<(size_t)t->nr*t->nc;i++)W[i]=.5*(W[i]+B.U[i]);}
  double*sf=S+t->soff;
  pdtrsm_("L","L","N","N",&t->n,&t->n,&a,sf,&one,&one,t->desc,W,&one,&one,t->desc);
  pdtrsm_("R","L","T","N",&t->n,&t->n,&a,sf,&one,&one,t->desc,W,&one,&one,t->desc);
  pdtran_(&t->n,&t->n,&a,W,&one,&one,t->desc,&b,B.U,&one,&one,t->desc);for(size_t i=0;i<(size_t)t->nr*t->nc;i++)W[i]=.5*(W[i]+B.U[i]);double t1=MPI_Wtime();
  int wanted=partial?(maxn<t->n?maxn:t->n):t->n,nev=wanted;elpa_solve_evp_real_2stage_double_impl_(&t->n,&nev,W,&t->nr,eval+t->ofs,B.U,&t->nr,&t->nb,&t->nc,&rows,&cols,&world);if(nev!=wanted)die("ELPA changed returned state count");for(int q=wanted;q<t->n;q++)eval[t->ofs+q]=DBL_MAX;MPI_Bcast(eval+t->ofs,t->n,MPI_DOUBLE,0,mpi_comm_level1);double t2=MPI_Wtime();
  pdtrsm_("L","L","T","N",&t->n,&wanted,&a,sf,&one,&one,t->desc,B.U,&one,&one,t->desc);putblock_rect(B.U,t,wanted,H,descH,t->ofs);
  if(d3_basis()&&k==0){putblock_rect(B.U,t,wanted,H,descH,B.sec[1].ofs);memcpy(eval+B.sec[1].ofs,eval+t->ofs,t->n*8);}
  double t3=MPI_Wtime();if(phase){int r;MPI_Comm_rank(mpi_comm_level1,&r);printf("CPU_NSOC_SYM_SECTOR rank=%d step=%d sector=%d n=%d nev=%d transform=%.6f eigen=%.6f recover=%.6f\n",r,step,k,t->n,wanted,t1-t0,t2-t1,t3-t2);}
 }
 double tE=MPI_Wtime();backtransform(step,n,H,W,phase);
 double tBack=MPI_Wtime();reorder(eval,ko,maxn,W,H);double end=MPI_Wtime();
 if(option("OPENMX_CPU_NSOC_CHECK",0)!=0)sampled_verify(step,maxn,ko,eval,W,H,nh,olp);
 if(step==1)dump_small(H,NULL,ko,maxn,"nsym_projected_eigen.bin",step);
 if(phase){int r;MPI_Comm_rank(mpi_comm_level1,&r);printf("CPU_NSOC_SYM_PHASE rank=%d step=%d overlap=%.6f rotate=%.6f sectors=%.6f backtransform=%.6f reorder=%.6f total=%.6f\n",r,step,tS-start,tH-tS,tE-tH,tBack-tE,end-tBack,end-start);fflush(stdout);}

 free(eval);MPI_Comm_free(&rowcomm);MPI_Comm_free(&colcomm);
}

static void solve_sector_once(int k,int wanted,double*Hsrc,double*Hdst,double*S,double*eval,double*work,int rows,int cols,int world,int one,double a,double b,int phase,int step){
 Sector*t=B.sec+k;double t0=MPI_Wtime();int packed=CpuNSym_NativeReturnEnabled();
 if(packed)memcpy(work,Hsrc+B.hoff[k],(size_t)t->nr*t->nc*sizeof(double));else copyblock(Hsrc,descH,t,work);
 if(d3_basis()&&k==0){if(packed)memcpy(B.U,Hsrc+B.hoff[1],(size_t)t->nr*t->nc*sizeof(double));else copyblock(Hsrc,descH,B.sec+1,B.U);if(B.version==3)paired_verify(paired_relative(work,B.U,(size_t)t->nr*t->nc,1.),"D3 Hamiltonian equivalent blocks disagree");for(size_t i=0;i<(size_t)t->nr*t->nc;i++)work[i]=.5*(work[i]+B.U[i]);}
 double*sf=S+t->soff;
 pdtrsm_("L","L","N","N",&t->n,&t->n,&a,sf,&one,&one,t->desc,work,&one,&one,t->desc);
 pdtrsm_("R","L","T","N",&t->n,&t->n,&a,sf,&one,&one,t->desc,work,&one,&one,t->desc);
 pdtran_(&t->n,&t->n,&a,work,&one,&one,t->desc,&b,B.U,&one,&one,t->desc);for(size_t i=0;i<(size_t)t->nr*t->nc;i++)work[i]=.5*(work[i]+B.U[i]);double t1=MPI_Wtime();
 int nev=wanted;elpa_solve_evp_real_2stage_double_impl_(&t->n,&nev,work,&t->nr,eval+t->ofs,B.U,&t->nr,&t->nb,&t->nc,&rows,&cols,&world);if(nev!=wanted)die("ELPA changed returned state count");for(int q=wanted;q<t->n;q++)eval[t->ofs+q]=DBL_MAX;MPI_Bcast(eval+t->ofs,t->n,MPI_DOUBLE,0,mpi_comm_level1);double t2=MPI_Wtime();
 pdtrsm_("L","L","T","N",&t->n,&wanted,&a,sf,&one,&one,t->desc,B.U,&one,&one,t->desc);
 if(packed)memcpy(Hdst+B.hoff[k],B.U,(size_t)t->nr*t->nc*sizeof(double));else putblock_rect(B.U,t,wanted,Hdst,descH,t->ofs);
 if(d3_basis()&&k==0){if(packed)memcpy(Hdst+B.hoff[1],B.U,(size_t)t->nr*t->nc*sizeof(double));else putblock_rect(B.U,t,wanted,Hdst,descH,B.sec[1].ofs);memcpy(eval+B.sec[1].ofs,eval+t->ofs,t->n*8);}
 double t3=MPI_Wtime();if(phase){int r;MPI_Comm_rank(mpi_comm_level1,&r);printf("CPU_NSOC_SYM_SECTOR rank=%d step=%d sector=%d n=%d nev=%d transform=%.6f eigen=%.6f recover=%.6f\n",r,step,k,t->n,wanted,t1-t0,t2-t1,t3-t2);fflush(stdout);}
}

static void solve_c3_complex(int step,double*H,double*S,double*eval,double*work,int rows,int cols,int world,int phase){
 Sector*t=B.sec;int one=1,x=t->ofs+1,y=B.sec[1].ofs+1,m=t->n,nev=m;
 size_t count=(size_t)t->nr*t->nc;dcomplex alpha={1.,0.},zero={0.,0.};double start=MPI_Wtime();
 c3_complex_extract(H,descH,work,B.U,step,"H");dcomplex*factor=(dcomplex*)(S+t->soff);
 pztrsm_("L","L","N","N",&m,&m,&alpha,factor,&one,&one,t->desc,(dcomplex*)work,&one,&one,t->desc);
 pztrsm_("R","L","C","N",&m,&m,&alpha,factor,&one,&one,t->desc,(dcomplex*)work,&one,&one,t->desc);
 pztranc_(&m,&m,&alpha,(dcomplex*)work,&one,&one,t->desc,&zero,(dcomplex*)B.U,&one,&one,t->desc);
 for(size_t i=0;i<2*count;i++)work[i]=.5*(work[i]+B.U[i]);double transformed=MPI_Wtime();
 elpa_solve_evp_complex_2stage_double_impl_(&m,&nev,(dcomplex*)work,&t->nr,eval+t->ofs,(dcomplex*)B.U,&t->nr,&t->nb,&t->nc,&rows,&cols,&world);
 if(nev!=m)die("C3 complex ELPA changed state count");MPI_Bcast(eval+t->ofs,m,MPI_DOUBLE,0,mpi_comm_level1);double eigen=MPI_Wtime();
 pztrsm_("L","L","C","N",&m,&m,&alpha,factor,&one,&one,t->desc,(dcomplex*)B.U,&one,&one,t->desc);
 for(size_t i=0;i<count;i++){work[i]=B.U[2*i];work[count+i]=B.U[2*i+1];}
 pdgemr2d_(&m,&m,work,&one,&one,t->desc,H,&x,&x,descH,&ictxt1);
 pdgemr2d_(&m,&m,work+count,&one,&one,t->desc,H,&y,&x,descH,&ictxt1);
 for(size_t i=0;i<count;i++)work[count+i]=-work[count+i];
 pdgemr2d_(&m,&m,work+count,&one,&one,t->desc,H,&x,&y,descH,&ictxt1);
 pdgemr2d_(&m,&m,work,&one,&one,t->desc,H,&y,&y,descH,&ictxt1);
 memcpy(eval+B.sec[1].ofs,eval+t->ofs,m*sizeof(double));
 if(phase){int rank;MPI_Comm_rank(mpi_comm_level1,&rank);printf("CPU_NSOC_C3_COMPLEX_PHASE rank=%d step=%d n=%d transform=%.9f eigen=%.9f recover=%.9f total=%.9f\n",rank,step,m,transformed-start,eigen-transformed,MPI_Wtime()-eigen,MPI_Wtime()-start);fflush(stdout);}
}

static int partial_grow(int old,int n){if(old>=n)return n;int target=old+1;if(old>1&&old*2>target)target=old*2;if(target>n)target=n;return target;}
static int partial_weighted_count(int*solved){return 2*solved[0]+solved[2]+solved[3];}
static void partial_seed_from_full(double*eval,int maxn){
 int total=2*B.sec[0].n+B.sec[2].n+B.sec[3].n;Eigen*list=mem((size_t)total,sizeof(Eigen));int z=0;
 for(int k=0;k<4;k++)if(k!=1){int w=k==0?2:1;for(int q=0;q<B.sec[k].n;q++)for(int c=0;c<w;c++){list[z].value=eval[B.sec[k].ofs+q];list[z].source=k;z++;}}
 qsort(list,total,sizeof(Eigen),eigcmp);for(int k=0;k<4;k++)B.partial_cache[k]=0;
 int upto=maxn<total?maxn:total;for(int i=0;i<upto;i++)B.partial_cache[list[i].source]++;
 /* A small guard absorbs normal SCF eigenvalue crossings; adaptive_boundary
  * still grows a sector when the guard is insufficient. */
 for(int k=0;k<4;k++)if(k!=1){int g=B.partial_cache[k]/20+8;B.partial_cache[k]+=g;if(B.partial_cache[k]>B.sec[k].n)B.partial_cache[k]=B.sec[k].n;}
 B.partial_cache[1]=B.partial_cache[0];B.partial_cache_valid=1;free(list);
}
static int partial_boundary(double*eval,int*solved,int maxn,double*boundary){
 int count=partial_weighted_count(solved);if(count<maxn)return 0;Eigen*list=mem((size_t)count,sizeof(Eigen));int z=0;
 for(int k=0;k<4;k++)if(k!=1){int w=k==0?2:1;for(int q=0;q<solved[k];q++)for(int c=0;c<w;c++){list[z].value=eval[B.sec[k].ofs+q];list[z].source=k;z++;}}
 qsort(list,count,sizeof(Eigen),eigcmp);*boundary=list[maxn-1].value;free(list);return isfinite(*boundary);
}

static int representative_cut_compare(const void*a,const void*b){int x=*(const int*)a,y=*(const int*)b;return(x>y)-(x<y);}
static void representative_record(int n,int maxn,const double*eval){
 representative_certificate_ready=0;representative_tail_ready=0;representative_cut_count=0;
 const char*flag=getenv("OPENMX_CPU_NSOC_REPRESENTATIVE_DM");if(!flag||!strcmp(flag,"0"))return;
 if(!CpuNSym_FullProjected()||!paired_basis()||option("OPENMX_CPU_NSOC_PARTIAL_EVP",0)!=0)return;
 const char*tail=getenv("OPENMX_CPU_NSOC_REP_TAIL_BOUND");int wanted=tail&&!strcmp(tail,"1");
 if(wanted&&n>representative_cut_cap){free(representative_cut_ids);representative_cut_ids=mem(n,sizeof(int));representative_cut_cap=n;}
 Eigen*sorted=mem(n,sizeof(Eigen));int*rank=mem(n,sizeof(int));
 for(int i=0;i<n;i++){if(!isfinite(eval[i])){free(sorted);free(rank);return;}sorted[i]=(Eigen){eval[i],i};}
 qsort(sorted,n,sizeof(Eigen),eigcmp);for(int i=0;i<n;i++)rank[sorted[i].source]=i;
 int okay=1,pairs_valid=1;for(int i=0;i<B.sec[0].n;i++){
  int a=B.sec[0].ofs+i,b=B.sec[1].ofs+i;
  if(eval[a]!=eval[b]){okay=0;pairs_valid=0;}
  if((rank[a]<maxn)!=(rank[b]<maxn)){okay=0;if(wanted)representative_cut_ids[representative_cut_count++]=(rank[a]<maxn?rank[a]:rank[b])+1;}
 }
 if(wanted&&pairs_valid&&representative_cut_count>0){qsort(representative_cut_ids,representative_cut_count,sizeof(int),representative_cut_compare);representative_tail_ready=1;}
 representative_certificate_ready=okay;representative_certificate_n=n;representative_certificate_maxn=maxn;free(sorted);free(rank);
}
#include "cpu_nsoc_representatives.h"
#include "cpu_nsoc_rep_intertwiner.inc"
static int representative_ao_checked=0,representative_ao_valid=0;
static int representative_context_valid(int n,int maxn,const double*ko,int ng,int nat,const int*const*perms,const int*types,const int*norb,const double*const*rotations){
 if(!prepared||!warm||!CpuNSym_FullProjected()||n!=representative_certificate_n||maxn!=representative_certificate_maxn||nat!=B.natoms||!ko||!perms||!types||!norb||!rotations)return 0;
 if((d3_basis()&&ng!=6)||(c3_basis()&&ng!=3))return 0;
 if(!representative_ao_checked){representative_ao_valid=representative_intertwiner(ng,nat,perms,types,norb,rotations);representative_ao_checked=1;}
 if(!representative_ao_valid)return 0;
 for(int i=1;i<=maxn;i++)if(!isfinite(ko[i]))return 0;
 return 1;
}
int CpuNSym_RepresentativeCertificate(int n,int maxn,const double*ko,int ng,int nat,const int*const*perms,const int*types,const int*norb,const double*const*rotations){
 return representative_certificate_ready&&representative_context_valid(n,maxn,ko,ng,nat,perms,types,norb,rotations);
}
int CpuNSym_RepresentativeTailCandidate(int n,int maxn,const double*ko,int ng,int nat,const int*const*perms,const int*types,const int*norb,const double*const*rotations,int*count,const int**ids){
 if(!representative_tail_ready||!representative_context_valid(n,maxn,ko,ng,nat,perms,types,norb,rotations))return 0;
 *count=representative_cut_count;*ids=representative_cut_ids;return 1;
}

/* SECTOR_PARALLEL_PATCH_BEGIN */
#include "cpu_nsoc_sector_parallel.inc"
/* SECTOR_PARALLEL_PATCH_END */
#include "native_panel_solver.inc"
void CpuNSym_Solve(int step,int n,int maxn,double*ko,double*S,double*H,double*W,double*****nh,double****olp){
 representative_certificate_ready=0;representative_tail_ready=0;representative_cut_count=0;
 if(!CpuNSym_Prepare(step)||maxn<1||maxn>n)die("invalid projected solver entry");double start=MPI_Wtime();basis_init(n);int phase=option("OPENMX_CPU_NSOC_PHASE",0)!=0;double partial_opt=option("OPENMX_CPU_NSOC_PARTIAL_EVP",0);if(partial_opt!=0&&partial_opt!=1&&partial_opt!=2)die("PARTIAL_EVP must be0, 1, or 2");/* Mode 1 keeps the original proportional adaptive seed; mode 2 reuses the previous full-step sector distribution before applying the same exact boundary test. */if(partial_opt!=0)die("CG1 requires PARTIAL_EVP=0");int warm_partial=partial_opt==2&&B.partial_cache_valid&&step>1;if(warm_partial)partial_seed_from_full(B.partial_eval_cache,maxn);int partial=partial_opt==1||warm_partial;
/* SECTOR_PARALLEL_PATCH_BEGIN */
 int parallel=sector_parallel_request(step);if(parallel)sector_parallel_prepare(S,H,W);
 int cfs_mode=CfsNativeRequest();(void)CfsNativeSelectedRequest(cfs_mode);
 if(CpuNSym_NativeReturnEnabled()){
  if(!d3_basis()||B.version!=3||!direct_active||direct_mode!=1||parallel||cfs_mode!=1||
     CfsNativeSelectedRequest(cfs_mode)!=1||!compact_eigen_scratch()||
     option("OPENMX_CPU_NSOC_DUMP",0)!=0||native_return.step!=step||native_return.coeff!=B.U)
   die("NATIVE_RETURN requires registered borrowed EVec1, direct real D3 version3, full-projected full-spectrum, selected/fused return1, no sector parallel or dump");
  if(native_return_mode()>=3&&(!CpuNSOCOwnerMode()||MO_fileout||level_fileout>=2||Dos_fileout||DosGauss_fileout||
      flag_export_DM||CWF_Calc||CWF_fileout_flag||Wannier_Func_Calc||Wannier90_fileout||
      pop_anal_aow_flag||unfold_electronic_band||Polarization_ABC_flag||xanes_calc||xanes_gs_fileout))
   die("state coefficient provider requires owner DM without eigenstate-based output and level.fileout<2");
  /* Previous iteration's native coefficients and packed return input have
   * no readers after its DM/EDM finished. Their current-step first writes
   * follow below; discard stale pages before sector solving, not afterwards. */
  size_t old_u=CpuReclaim_Pages(B.U,0,eigen_scratch.capacity*sizeof(double));
  size_t old_w=CpuReclaim_Pages(W,0,B.cap*sizeof(double));
  if(phase){int rank;MPI_Comm_rank(mpi_comm_level1,&rank);printf("CPU_NSOC_NATIVE_DEAD_SCRATCH rank=%d step=%d old_coeff_pages_returned=%zu old_return_input_pages_returned=%zu previous_DM_complete=1\n",rank,step,old_u,old_w);fflush(stdout);}
 }
/* SECTOR_PARALLEL_PATCH_END */

 if(rank0()){printf("CPU_NSOC_SYM_PROJECT step=%d n=%d sectors=%d,%d,%d,%d maxn=%d evp_mode=%s\n",step,n,B.sz[0],B.sz[1],B.sz[2],B.sz[3],maxn,warm_partial?"partial_warm_adaptive":(partial?"partial_adaptive":"full"));fflush(stdout);}
 if(step==1&&(!direct_active||direct_mode==2))dump_small(H,S,NULL,maxn,"nsym_raw_HS.bin",step);
 if(step==1)overlap_factors(S,W);
/* SECTOR_PARALLEL_PATCH_BEGIN */
 if(parallel){if(step==1)sector_parallel_migrate_factors(S,W);else if(!sector_parallel.migrated)die("SECTOR_PARALLEL factor cache unavailable");}
/* SECTOR_PARALLEL_PATCH_END */
double tS=MPI_Wtime();double hdef;
 if(direct_active){if(direct_mode==2){congruence(H,descH,W,descC);hdef=offblock_norm(H);}else hdef=0.;direct_assemble(step,1,H,descH,W);}
 else{congruence(H,descH,W,descC);hdef=offblock_norm(H);}
double strict_h=option("OPENMX_CPU_NSOC_GENERIC_COV_STRICT",1);double cov_tol=option("OPENMX_CPU_NSOC_GENERIC_COV_TOL",1e-5);if(!isfinite(cov_tol)||cov_tol<=0)die("GENERIC_COV_TOL must be positive");if(!isfinite(hdef)||(B.version>=2&&strict_h!=0&&hdef>cov_tol))die("raw Hamiltonian fails sector covariance; check the validated basis and operation convention");if(rank0()&&step==1){
 if(direct_active&&direct_mode==1)printf("CPU_NSOC_SYM_COVARIANCE overlap_checked=1 h_covariance=construction_certificate_actual_D3_AO_intertwiner_distinct_A_fresh_Reynolds raw_offblock_measurement=not_performed strict=%.0f\n",strict_h);
 else printf("CPU_NSOC_SYM_COVARIANCE overlap_checked=1 h_offblock_relative=%.12e strict=%.0f\n",hdef,strict_h);fflush(stdout);
}double tH=MPI_Wtime();
 double*Hsrc=H;if(partial){Hsrc=mem(B.active,8);memcpy(Hsrc,H,B.active*8);}double*eval=mem(n,8);for(int i=0;i<n;i++)eval[i]=DBL_MAX;
 int one=1;double a=1,b=0;MPI_Comm rowcomm,colcomm;MPI_Comm_split(mpi_comm_level1,my_pcol,my_prow,&rowcomm);MPI_Comm_split(mpi_comm_level1,my_prow,my_pcol,&colcomm);int rows=MPI_Comm_c2f(rowcomm),cols=MPI_Comm_c2f(colcomm),world=MPI_Comm_c2f(mpi_comm_level1);
 int solved[NSYM_MAX_SECTORS]={0},dirty[NSYM_MAX_SECTORS]={0};int denom=2*B.sec[0].n+B.sec[2].n+B.sec[3].n;
 if(!d3_basis() && partial)die("PARTIAL_EVP requires the validated four-sector D3 basis");
/* SECTOR_PARALLEL_PATCH_BEGIN */
 if(parallel){for(int k=0;k<4;k++)solved[k]=B.sec[k].n;sector_parallel_solve(step,H,S,W,eval,phase);}else
/* SECTOR_PARALLEL_PATCH_END */
 if(c3_basis()){
  if(partial)die("C3 complex path requires complete conjugate sectors");
  solve_c3_complex(step,H,S,eval,W,rows,cols,world,phase);
  solve_sector_once(2,B.sec[2].n,Hsrc,H,S,eval,W,rows,cols,world,one,a,b,phase,step);
 }else if(d3_basis()){
  dirty[0]=dirty[2]=dirty[3]=1;
  if(warm_partial){for(int k=0;k<4;k++)solved[k]=B.partial_cache[k];}else for(int k=0;k<4;k++)solved[k]=partial?((maxn*B.sec[k].n+denom-1)/denom):B.sec[k].n;solved[1]=solved[0];
  for(;;){
   for(int k=0;k<4;k++)if(k!=1&&dirty[k]){if(solved[k]<1)die("empty adaptive sector");solve_sector_once(k,solved[k],Hsrc,H,S,eval,W,rows,cols,world,one,a,b,phase,step);dirty[k]=0;}
   if(!partial)break;
   int total=partial_weighted_count(solved),grew=0;double boundary=0;
   if(total<maxn){for(int k=0;k<4;k++)if(k!=1&&solved[k]<B.sec[k].n){solved[k]=partial_grow(solved[k],B.sec[k].n);dirty[k]=1;grew=1;}if(grew)continue;die("adaptive partial EVP cannot provide requested states");}
   if(!partial_boundary(eval,solved,maxn,&boundary))die("adaptive partial EVP boundary unavailable");
   double tol=1e-12*fmax(1.0,fabs(boundary));
   for(int k=0;k<4;k++)if(k!=1&&solved[k]<B.sec[k].n){double last=eval[B.sec[k].ofs+solved[k]-1];if(!isfinite(last)||last<=boundary+tol){solved[k]=partial_grow(solved[k],B.sec[k].n);dirty[k]=1;grew=1;}}
   if(!grew)break;
  }
 }else{
  if(partial)die("PARTIAL_EVP is currently implemented only for the legacy D3 basis");
  for(int k=0;k<B.nsec;k++){
   solved[k]=B.sec[k].n;dirty[k]=1;
   solve_sector_once(k,solved[k],Hsrc,H,S,eval,W,rows,cols,world,one,a,b,phase,step);
  }
 }
 if(d3_basis())solved[1]=solved[0];
 if(rank0()&&partial){printf("CPU_NSOC_SYM_PARTIAL_COUNTS solved=%d,%d,%d,%d weighted=%d requested=%d\n",solved[0],solved[1],solved[2],solved[3],partial_weighted_count(solved),maxn);fflush(stdout);}
 double tE=MPI_Wtime(),tBack;
 if(CpuNSym_NativeReturnEnabled()){CpuNativeReturnTail(step,n,maxn,eval,ko,H,W,phase);tBack=MPI_Wtime();}
 else if(cfs_mode){CfsNativeTail(step,n,maxn,eval,ko,H,W,phase,cfs_mode);tBack=MPI_Wtime();}
 else{backtransform(step,n,H,W,phase);tBack=MPI_Wtime();reorder(eval,ko,maxn,W,H);}
 representative_record(n,maxn,eval);double end=MPI_Wtime();
 if(option("OPENMX_CPU_NSOC_CHECK",0)!=0)sampled_verify(step,maxn,ko,eval,W,H,nh,olp);
 if(step==1&&!CpuNSym_NativeReturnEnabled())dump_small(H,NULL,ko,maxn,"nsym_projected_eigen.bin",step);
 if(phase){int r;MPI_Comm_rank(mpi_comm_level1,&r);printf("CPU_NSOC_SYM_PHASE rank=%d step=%d overlap=%.6f rotate=%.6f sectors=%.6f backtransform=%.6f reorder=%.6f total=%.6f\n",r,step,tS-start,tH-tS,tE-tH,tBack-tE,end-tBack,end-start);fflush(stdout);}
 if(partial_opt==2&&d3_basis()&&!partial){if(!B.partial_eval_cache)B.partial_eval_cache=mem(n,8);memcpy(B.partial_eval_cache,eval,(size_t)n*8);B.partial_cache_valid=1;}
 if(partial)free(Hsrc);free(eval);MPI_Comm_free(&rowcomm);MPI_Comm_free(&colcomm);
}
