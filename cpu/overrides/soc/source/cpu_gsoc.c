/* GSOC v1/v2: full native spinors from distributed real/complex symmetry sectors.
 * Independent environment gate; no C3/D3 selector, geometry or compact aliasing. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "mpi.h"
#include "openmx_common.h"
#include "cpu_gsoc.h"
#include "cpu_soc_native.h"
#include "gsoc_export_audit.h"
#include "cpu_gsoc_compact.h"
#include "gsoc_index.h"
#include "gsoc_schema.h"
#include "../audit/gsoc_operator_gate.h"
extern int numroc_(int*,int*,int*,int*,int*);
extern void descinit_(int*,int*,int*,int*,int*,int*,int*,int*,int*,int*);
extern void pdpotrf_(char*,int*,double*,int*,int*,int*,int*);
extern void pdsygst_(int*,char*,int*,double*,int*,int*,int*,double*,int*,int*,int*,double*,int*);
extern void pdtran_(int*,int*,double*,double*,int*,int*,int*,double*,double*,int*,int*,int*);
extern void pdtrsm_(char*,char*,char*,char*,int*,int*,double*,double*,int*,int*,int*,double*,int*,int*,int*);
extern void pzpotrf_(char*,int*,dcomplex*,int*,int*,int*,int*);
extern void pzhegst_(int*,char*,int*,dcomplex*,int*,int*,int*,dcomplex*,int*,int*,int*,double*,int*);
extern void pztranc_(int*,int*,dcomplex*,dcomplex*,int*,int*,int*,dcomplex*,dcomplex*,int*,int*,int*);
extern void pztrsm_(char*,char*,char*,char*,int*,int*,dcomplex*,dcomplex*,int*,int*,int*,dcomplex*,int*,int*,int*);
extern void elpa_solve_evp_real_2stage_double_impl_(int*,int*,double*,int*,double*,double*,int*,int*,int*,int*,int*,int*);
extern void elpa_solve_evp_complex_2stage_double_impl_(int*,int*,dcomplex*,int*,double*,dcomplex*,int*,int*,int*,int*,int*,int*);
extern double CpuC3_PeriodicDistanceBohr(const double*,const double*);
typedef char gsoc_complex_size_check[sizeof(dcomplex)==sizeof(GsocAuditComplex)?1:-1];
typedef struct {int n,ofs,dtype,partner,child,nr,nc,nb,desc[9];void*Q;} GSec;
typedef struct {int band,row;dcomplex z;} GRecord;
typedef struct {int a,b,ga,gb,na,nb;} GBlock;
static struct {int loaded,version,improper,n,ns,nat,nsec,rank,np,*offset,*dtype,*mp,*rp,*ci;int64_t nnz,*cp;int*ri;dcomplex*cv,*rv;GSec sec[48];double*****nh,*****im;double****olp;int*is2,*ie2;double*ko;dcomplex*ev;int projected;double tolerance;int nblocks;GBlock*blocks;} G;
static int gsoc_restart_read=-1;
void CpuGSOC_RestartState(int ok){gsoc_restart_read=ok==1?1:0;}
static uint64_t ghash(uint64_t x){x+=UINT64_C(0x9e3779b97f4a7c15);x=(x^(x>>30))*UINT64_C(0xbf58476d1ce4e5b9);x=(x^(x>>27))*UINT64_C(0x94d049bb133111eb);return x^(x>>31);}
static dcomplex gz(double r,double i){dcomplex z;z.r=r;z.i=i;return z;}
static dcomplex ga(dcomplex a,dcomplex b){return gz(a.r+b.r,a.i+b.i);}
static dcomplex gm(dcomplex a,dcomplex b){return gz(a.r*b.r-a.i*b.i,a.r*b.i+a.i*b.r);}
static dcomplex gc(dcomplex a){return gz(a.r,-a.i);}
static void gdie(const char*s){int rank;MPI_Comm_rank(mpi_comm_level1,&rank);fprintf(stderr,"CPU_GSOC_ERROR rank=%d %s\n",rank,s);fflush(stderr);MPI_Abort(mpi_comm_level1,94);abort();}
static void* gmem(size_t n,size_t bytes){if(bytes&&n>SIZE_MAX/bytes)gdie("allocation integer overflow");void*p=calloc(n?n:1,bytes);if(!p)gdie("allocation failed");return p;}
static void grd(void*p,size_t s,size_t n,FILE*f){if(fread(p,s,n,f)!=n)gdie("short GSOC payload");}
int CpuGSOC_Enabled(void){const char*s=getenv("OPENMX_CPU_GSOC");if(!s||!strcmp(s,"0"))return 0;if(strcmp(s,"1"))gdie("GSOC must be exactly0/1");return 1;}
static int genv(const char*key,int fallback,int low,int high){const char*s=getenv(key);if(!s)return fallback;char*end;long v=strtol(s,&end,10);if(*end||v<low||v>high)gdie("invalid bounded integer GSOC option");return(int)v;}
static double grealenv(const char*key,double fallback,double low,double high){const char*s=getenv(key);if(!s)return fallback;char*end;double v=strtod(s,&end);if(*end||!isfinite(v)||v<low||v>high)gdie("invalid finite GSOC tolerance");return v;}
static int gforward(void*ctx,const GsocAuditComplex*xx,GsocAuditComplex*yy){const dcomplex*x=(const dcomplex*)xx;dcomplex*y=(dcomplex*)yy;memset(y,0,(size_t)G.ns*sizeof(dcomplex));for(int j=0;j<G.ns;j++)for(int64_t a=G.cp[j];a<G.cp[j+1];a++)y[G.ri[a]]=ga(y[G.ri[a]],gm(G.cv[a],x[j]));return 0;}
static int ginverse(void*ctx,const GsocAuditComplex*xx,GsocAuditComplex*yy){const dcomplex*x=(const dcomplex*)xx;dcomplex*y=(dcomplex*)yy;for(int j=0;j<G.ns;j++){dcomplex z=gz(0,0);for(int64_t a=G.cp[j];a<G.cp[j+1];a++)z=ga(z,gm(gc(G.cv[a]),x[G.ri[a]]));y[j]=z;}return 0;}
static int gapply(void*ctx,const GsocAuditComplex*xx,GsocAuditComplex*hh,GsocAuditComplex*ss){const dcomplex*x=(const dcomplex*)xx;dcomplex*h=(dcomplex*)hh,*s=(dcomplex*)ss;memset(h,0,(size_t)G.ns*sizeof(dcomplex));memset(s,0,(size_t)G.ns*sizeof(dcomplex));
 for(int a=1;a<=Matomnum;a++){int ga0=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga0]];for(int b=0;b<=FNAN[ga0];b++){int gb=natn[ga0][b],nb=Spe_Total_CNO[WhatSpecies[gb]];for(int i=0;i<na;i++)for(int j=0;j<nb;j++){int I=G.mp[ga0]+i,J=G.mp[gb]+j;double ov=G.olp[a][b][i][j];dcomplex uu=gz(G.nh[0][a][b][i][j],G.im[0][a][b][i][j]),dd=gz(G.nh[1][a][b][i][j],G.im[1][a][b][i][j]),ud=gz(G.nh[2][a][b][i][j],G.nh[3][a][b][i][j]+G.im[2][a][b][i][j]);h[I]=ga(h[I],gm(uu,x[J]));h[I+G.n]=ga(h[I+G.n],gm(dd,x[J+G.n]));h[I]=ga(h[I],gm(ud,x[J+G.n]));h[J+G.n]=ga(h[J+G.n],gm(gc(ud),x[I]));s[I]=ga(s[I],gz(ov*x[J].r,ov*x[J].i));s[I+G.n]=ga(s[I+G.n],gz(ov*x[J+G.n].r,ov*x[J+G.n].i));}}}
 MPI_Allreduce(MPI_IN_PLACE,h,2*G.ns,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);MPI_Allreduce(MPI_IN_PLACE,s,2*G.ns,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);return 0;}
/* Exact Hermitian assembly action used by gentry, independent of sector code.
 * Directed UD is already paired with its conjugate reverse in native storage. */
static int gapply_hermitian(void*ctx,const GsocAuditComplex*xx,GsocAuditComplex*hh,GsocAuditComplex*ss){
 const dcomplex*x=(const dcomplex*)xx;dcomplex*h=(dcomplex*)hh,*s=(dcomplex*)ss;
 memset(h,0,(size_t)G.ns*sizeof(*h));memset(s,0,(size_t)G.ns*sizeof(*s));
 for(int a=1;a<=Matomnum;a++){int ga0=M2G[a],na=Spe_Total_CNO[WhatSpecies[ga0]];
  for(int b=0;b<=FNAN[ga0];b++){int gb=natn[ga0][b],nb=Spe_Total_CNO[WhatSpecies[gb]];
   for(int i=0;i<na;i++)for(int j=0;j<nb;j++){
    int I=G.mp[ga0]+i,J=G.mp[gb]+j;double ov=.5*G.olp[a][b][i][j];
    dcomplex uu=gz(.5*G.nh[0][a][b][i][j],.5*G.im[0][a][b][i][j]),dd=gz(.5*G.nh[1][a][b][i][j],.5*G.im[1][a][b][i][j]);
    dcomplex ud=gz(G.nh[2][a][b][i][j],G.nh[3][a][b][i][j]+G.im[2][a][b][i][j]);
    h[I]=ga(h[I],gm(uu,x[J]));h[J]=ga(h[J],gm(gc(uu),x[I]));
    h[I+G.n]=ga(h[I+G.n],gm(dd,x[J+G.n]));h[J+G.n]=ga(h[J+G.n],gm(gc(dd),x[I+G.n]));
    h[I]=ga(h[I],gm(ud,x[J+G.n]));h[J+G.n]=ga(h[J+G.n],gm(gc(ud),x[I]));
    for(int spin=0;spin<2;spin++){int ii=I+spin*G.n,jj=J+spin*G.n;
     s[ii]=ga(s[ii],gz(ov*x[jj].r,ov*x[jj].i));s[jj]=ga(s[jj],gz(ov*x[ii].r,ov*x[ii].i));}
   }
  }
 }
 MPI_Allreduce(MPI_IN_PLACE,h,2*G.ns,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 MPI_Allreduce(MPI_IN_PLACE,s,2*G.ns,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);return 0;
}
static GsocAuditComplex *gexport_probe=NULL;static int gexport_count=0;
static int gexport_enabled(void){const char*s=getenv("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE");return s&&!strcmp(s,"1");}
static void gexport_reference(void){
 if(!gexport_enabled())return;
 int partners[48];for(int k=0;k<G.nsec;k++)partners[k]=G.sec[k].partner;
 GsocAuditOperator op={G.ns,G.nsec,G.offset,G.dtype,NULL,gapply_hermitian,gforward,ginverse};
 gexport_count=2;gexport_probe=gmem((size_t)gexport_count*3*G.ns,sizeof(*gexport_probe));
 uint64_t seed=UINT64_C(20261102);double scale=1./sqrt(2.*G.ns);
 for(int p=0;p<gexport_count;p++){
  GsocAuditComplex*x=gexport_probe+(size_t)p*3*G.ns,*h=x+G.ns,*s=h+G.ns;
  for(int i=0;i<G.ns;i++)x[i]=gsoc_audit_c((gsoc_audit_rng(&seed)&1)?scale:-scale,(gsoc_audit_rng(&seed)&1)?scale:-scale);
  if(gsoc_export_projected_apply(&op,partners,x,h,s))gdie("cannot construct solved-operator export reference");
 }
}
void CpuGSOC_ValidateExport(void){
 if(!gexport_enabled()||!gexport_probe||!G.mp||gexport_count!=2)gdie("missing completed GSOC export reference");
 GsocAuditComplex*v=gmem((size_t)G.ns*2,sizeof(*v));double absolute[2]={0},relative[2]={0};int bad=0;
 for(int p=0;p<gexport_count;p++){
  GsocAuditComplex*x=gexport_probe+(size_t)p*3*G.ns;
  gapply(NULL,x,v,v+G.ns);
  gsoc_audit_metric(v,x+G.ns,G.ns,&absolute[0],&relative[0],&bad);
  gsoc_audit_metric(v+G.ns,x+2*G.ns,G.ns,&absolute[1],&relative[1],&bad);
 }
 if(!G.rank){printf("CPU_GSOC_EXPORT_OPERATOR_AUDIT probes=%d entrywise_proof=0 H_abs=%.12e H_rel=%.12e S_abs=%.12e S_rel=%.12e tolerance=1e-10 spatial_Reynolds_vs_solved_sector_TR_real=1\n",gexport_count,absolute[0],relative[0],absolute[1],relative[1]);fflush(stdout);}
 free(v);free(gexport_probe);gexport_probe=NULL;free(G.mp);G.mp=NULL;
 if(bad||relative[0]>1e-10||relative[1]>1e-10)gdie("spatial Reynolds export differs from actual solved sector/TR/real operator; unsupported projection");
}
static void gbasis_validate(void){
 int*stamp=gmem(G.ns,sizeof(int)),*touched=gmem(G.ns,sizeof(int));dcomplex*sum=gmem(G.ns,sizeof(dcomplex));double defect=0;
 for(int j=0;j<G.ns;j++){int used=0,tag=j+1;for(int64_t a=G.cp[j];a<G.cp[j+1];a++){int row=G.ri[a];for(int b=G.rp[row];b<G.rp[row+1];b++){int k=G.ci[b];if(stamp[k]!=tag){stamp[k]=tag;sum[k]=gz(0,0);touched[used++]=k;}sum[k]=ga(sum[k],gm(gc(G.cv[a]),G.rv[b]));}}for(int q=0;q<used;q++){int k=touched[q];dcomplex z=sum[k];if(k==j)z.r-=1;double er=hypot(z.r,z.i);if(er>defect)defect=er;}}
 if(defect>1e-10)gdie("full sparse GSOC basis Gram check failed");
 double tr=0;memset(sum,0,(size_t)G.ns*sizeof(dcomplex));memset(stamp,0,(size_t)G.ns*sizeof(int));
 for(int k=0;k<G.nsec;k++)if(G.sec[k].partner>=0){GSec*t=G.sec+k,*src=G.sec+t->partner;for(int j=0;j<t->n;j++){int tag=t->ofs+j+1,used=0;for(int side=0;side<2;side++){int col=(side?t->ofs:src->ofs)+j;for(int64_t a=G.cp[col];a<G.cp[col+1];a++){int row=G.ri[a];dcomplex z=G.cv[a];if(!side){z=gc(z);if(row<G.n){row+=G.n;z.r=-z.r;z.i=-z.i;}else row-=G.n;}else{z.r=-z.r;z.i=-z.i;}if(stamp[row]!=tag){stamp[row]=tag;sum[row]=gz(0,0);touched[used++]=row;}sum[row]=ga(sum[row],z);}}for(int q=0;q<used;q++){double er=hypot(sum[touched[q]].r,sum[touched[q]].i);if(er>tr)tr=er;}}}
 if(tr>1e-10)gdie("TR partner physical basis identity failed");free(stamp);free(touched);free(sum);if(!G.rank){printf("CPU_GSOC_BASIS_GRAM max=%.12e TR_partner_map=%.12e\n",defect,tr);fflush(stdout);}}
static void gload(void){
 if(G.loaded)gdie("GSOC supports one physical solve per process");MPI_Comm_rank(mpi_comm_level1,&G.rank);MPI_Comm_size(mpi_comm_level1,&G.np);const char*path=getenv("OPENMX_CPU_GSOC_BASIS");if(!path)gdie("missing GSOC_BASIS");FILE*f=fopen(path,"rb");if(!f)gdie("cannot open GSOC basis");int32_t hd[6];grd(hd,4,6,f);if(hd[0]!=0x47534f43||!gsoc_schema_version(hd[1])||hd[2]!=atomnum||hd[3]<1||hd[3]>INT_MAX/4||hd[4]<1||hd[4]>48||hd[5]<1||hd[5]>192)gdie("invalid GSOC v1/v2 header");G.version=hd[1];G.nat=hd[2];G.n=hd[3];G.ns=2*G.n;G.nsec=hd[4];grd(&G.nnz,8,1,f);if(G.nnz<1||G.nnz>INT_MAX)gdie("GSOC sparse index capacity");G.offset=gmem(G.nsec+1,sizeof(int));G.dtype=gmem(G.nsec,sizeof(int));int sizes[48],partners[48];grd(sizes,4,G.nsec,f);grd(G.dtype,4,G.nsec,f);grd(partners,4,G.nsec,f);
 for(int k=0;k<G.nsec;k++){GSec*t=G.sec+k;t->n=sizes[k];t->ofs=G.offset[k];t->dtype=G.dtype[k];t->partner=partners[k];t->child=-1;if(t->n<1||t->n>G.ns-t->ofs||(t->dtype!=1&&t->dtype!=2)||t->partner<-1||t->partner>=k)gdie("invalid GSOC sector schema");G.offset[k+1]=t->ofs+t->n;if(t->partner>=0){GSec*s=G.sec+t->partner;if(s->partner>=0||s->n!=t->n||s->dtype!=t->dtype||s->child>=0)gdie("invalid TR independent-sector graph");s->child=k;}}
 if(G.offset[G.nsec]!=G.ns)gdie("sector sizes do not span full spinor space");double cell[9],nativecell[9];grd(cell,8,9,f);for(int i=0;i<3;i++)for(int j=0;j<3;j++){nativecell[3*i+j]=tv[i+1][j+1];if(!isfinite(cell[3*i+j])||fabs(cell[3*i+j]-nativecell[3*i+j]*BohrR)>1e-8)gdie("basis cell differs from native input");}
 int*counts=gmem(G.nat,sizeof(int)),*species=gmem(G.nat,sizeof(int));grd(counts,4,G.nat,f);grd(species,4,G.nat,f);G.mp=gmem(atomnum+1,sizeof(int));int off=0;for(int a=1;a<=atomnum;a++){G.mp[a]=off;int n=Spe_Total_CNO[WhatSpecies[a]];if(counts[a-1]!=n||species[a-1]!=WhatSpecies[a])gdie("basis species/AO ordering mismatch");off+=n;}if(off!=G.n)gdie("native scalar AO count mismatch");free(counts);free(species);
 double gt=grealenv("OPENMX_CPU_GSOC_GEOMETRY_TOL_ANG",5e-5,1e-12,0.02),maxgeom=0;int*perm=gmem(G.nat,sizeof(int)),*seen=gmem(G.nat,sizeof(int));for(int op=0;op<hd[5];op++){double rot[9],shift[3];grd(rot,8,9,f);grd(shift,8,3,f);grd(perm,4,G.nat,f);memset(seen,0,(size_t)G.nat*sizeof(int));for(int i=0;i<9;i++)if(!isfinite(rot[i]))gdie("nonfinite geometry rotation");for(int i=0;i<3;i++){if(!isfinite(shift[i]))gdie("nonfinite geometry shift");for(int j=0;j<3;j++){double v=0;for(int l=0;l<3;l++)v+=rot[3*l+i]*rot[3*l+j];if(fabs(v-(i==j))>1e-10)gdie("nonorthogonal geometry operation");}}double determinant=rot[0]*(rot[4]*rot[8]-rot[5]*rot[7])-rot[1]*(rot[3]*rot[8]-rot[5]*rot[6])+rot[2]*(rot[3]*rot[7]-rot[4]*rot[6]);if(!gsoc_schema_determinant(G.version,determinant))gdie("geometry determinant violates declared GSOC schema version");if(determinant<0)G.improper++;for(int a=0;a<G.nat;a++){int b=perm[a];if(b<0||b>=G.nat||seen[b]++||WhatSpecies[a+1]!=WhatSpecies[b+1])gdie("invalid geometry atom bijection");double dr[3];for(int i=0;i<3;i++){double x=0;for(int j=0;j<3;j++)x+=rot[3*i+j]*Gxyz[a+1][j+1]+shift[j]*tv[j+1][i+1];dr[i]=Gxyz[b+1][i+1]-x;}double er=CpuC3_PeriodicDistanceBohr(dr,nativecell)*BohrR;if(!isfinite(er)||er>gt)gdie("geometry exceeds declared GSOC tolerance");if(er>maxgeom)maxgeom=er;}}free(perm);free(seen);
 G.cp=gmem(G.ns+1,sizeof(int64_t));G.ri=gmem(G.nnz,sizeof(int));G.cv=gmem(G.nnz,sizeof(dcomplex));grd(G.cp,8,G.ns+1,f);grd(G.ri,4,G.nnz,f);grd(G.cv,sizeof(dcomplex),G.nnz,f);if(fgetc(f)!=EOF)gdie("trailing GSOC payload");fclose(f);if(G.cp[0]!=0||G.cp[G.ns]!=G.nnz)gdie("invalid CSC extents");G.rp=gmem(G.ns+1,sizeof(int));G.ci=gmem(G.nnz,sizeof(int));G.rv=gmem(G.nnz,sizeof(dcomplex));for(int j=0;j<G.ns;j++){if(G.cp[j]<0||G.cp[j]>=G.cp[j+1]||G.cp[j+1]>G.nnz)gdie("invalid CSC pointer");for(int64_t a=G.cp[j];a<G.cp[j+1];a++){int i=G.ri[a];if(i<0||i>=G.ns||!isfinite(G.cv[a].r)||!isfinite(G.cv[a].i))gdie("invalid CSC entry");G.rp[i+1]++;}}for(int i=0;i<G.ns;i++)G.rp[i+1]+=G.rp[i];int*pos=gmem(G.ns,sizeof(int));memcpy(pos,G.rp,(size_t)G.ns*sizeof(int));for(int j=0;j<G.ns;j++)for(int64_t a=G.cp[j];a<G.cp[j+1];a++){int i=G.ri[a],at=pos[i]++;G.ci[at]=j;G.rv[at]=G.cv[a];}free(pos);gbasis_validate();
 if(G.rank!=my_prow2*np_cols2+my_pcol2||G.np!=np_rows2*np_cols2)gdie("native BLACS/MPI mapping is not supported row-major layout");for(int k=0;k<G.nsec;k++){GSec*t=G.sec+k;int zero=0,info;t->nb=nblk2;while(t->nb>1&&(t->nb*np_rows2>t->n||t->nb*np_cols2>t->n))t->nb/=2;t->nr=numroc_(&t->n,&t->nb,&my_prow2,&zero,&np_rows2);t->nc=numroc_(&t->n,&t->nb,&my_pcol2,&zero,&np_cols2);if(t->nr<1||t->nc<1)gdie("too many ranks for nonempty ELPA sector");int ld=t->nr;descinit_(t->desc,&t->n,&t->n,&t->nb,&t->nb,&zero,&zero,&ictxt1_2,&ld,&info);if(info)gdie("sector descriptor creation");}
 G.loaded=1;if(!G.rank){printf("CPU_GSOC_BASIS version=%d improper_operations=%d spin_lift=%s scalar_N=%d spinor_N=%d sectors=%d nnz=%lld geometry_max_Ang=%.12e caller_dense_buffers_used=0\n",G.version,G.improper,G.version==2?"axial_detR_R":"proper_R",G.n,G.ns,G.nsec,(long long)G.nnz,maxgeom);fflush(stdout);}}
/* Contribution to one independent sector. A physical TR partner is averaged
 * into its representative before any real gauge projection. No full AO matrix. */
static void gentry_one(double*buf,int j0,int j1,GSec*t,int I,int J,dcomplex h,double s){
 int choices=t->child>=0?2:1;
 for(int choice=0;choice<choices;choice++){GSec*v=choice?G.sec+t->child:t;double weight=choices==2?0.5:1.;
  for(int b=G.rp[J];b<G.rp[J+1];b++){int col=G.ci[b]-v->ofs;if(col<j0||col>=j1)continue;
   for(int a=G.rp[I];a<G.rp[I+1];a++){int row=G.ci[a]-v->ofs;if(row<0||row>=t->n)continue;dcomplex z=gm(gc(G.rv[a]),G.rv[b]),hv=gm(z,h);if(choice){hv=gc(hv);z=gc(z);}size_t at=((size_t)(col-j0)*t->n+row)*4;buf[at]+=weight*hv.r;buf[at+1]+=weight*hv.i;buf[at+2]+=weight*z.r*s;buf[at+3]+=weight*z.i*s;}
  }
 }
}
/* Hermitian Reynolds projection is explicit and applied equally to H/S. */
static void gentry(double*buf,int j0,int j1,GSec*t,int I,int J,dcomplex h,double s){
 dcomplex half=gz(0.5*h.r,0.5*h.i);gentry_one(buf,j0,j1,t,I,J,half,0.5*s);gentry_one(buf,j0,j1,t,J,I,gc(half),0.5*s);
}
static int grank(int row,int col,GSec*t){return((row/t->nb)%np_rows2)*np_cols2+(col/t->nb)%np_cols2;}
/* Bucket native neighbor blocks by the reduced columns they can touch.
 * Both AO ends/spins are included because assembly explicitly Hermitizes. */
static void gblocks_init(void){
 int count=0;for(int a=1;a<=Matomnum;a++){int ga0=M2G[a];if(FNAN[ga0]>INT_MAX-count-1)gdie("native block count overflow");count+=FNAN[ga0]+1;}G.nblocks=count;G.blocks=gmem(count,sizeof(GBlock));int at=0;for(int a=1;a<=Matomnum;a++){int ga0=M2G[a];for(int b=0;b<=FNAN[ga0];b++){GBlock*v=G.blocks+at++;v->a=a;v->b=b;v->ga=ga0;v->gb=natn[ga0][b];v->na=Spe_Total_CNO[WhatSpecies[v->ga]];v->nb=Spe_Total_CNO[WhatSpecies[v->gb]];}}
}
static int gblock_tiles(GSec*t,GBlock*v,int tile,int*mark,int tag,int*touched){
 int used=0;for(int end=0;end<2;end++){int atom=end?v->gb:v->ga,n=end?v->nb:v->na;for(int spin=0;spin<2;spin++)for(int i=0;i<n;i++){int row=G.mp[atom]+i+spin*G.n;for(int a=G.rp[row];a<G.rp[row+1];a++){int global=G.ci[a],local=global-t->ofs;if(local<0||local>=t->n){if(t->child<0)continue;local=global-G.sec[t->child].ofs;if(local<0||local>=t->n)continue;}int id=local/tile;if(mark[id]!=tag){mark[id]=tag;touched[used++]=id;}}}}return used;
}
static void gblock_buckets(GSec*t,int tile,int**offset,int**ids){
 int ntile=(t->n+tile-1)/tile;int*mark=gmem(ntile,sizeof(int)),*touched=gmem(ntile,sizeof(int)),*ptr=gmem(ntile+1,sizeof(int));size_t total=0;for(int b=0;b<G.nblocks;b++){int used=gblock_tiles(t,G.blocks+b,tile,mark,b+1,touched);for(int q=0;q<used;q++){if(ptr[touched[q]+1]==INT_MAX)gdie("block bucket count overflow");ptr[touched[q]+1]++;}total+=used;if(total>INT_MAX)gdie("block bucket integer capacity");}for(int j=1;j<=ntile;j++)ptr[j]+=ptr[j-1];int*list=gmem(total,sizeof(int)),*pos=gmem(ntile,sizeof(int));memcpy(pos,ptr,(size_t)ntile*sizeof(int));memset(mark,0,(size_t)ntile*sizeof(int));for(int b=0;b<G.nblocks;b++){int used=gblock_tiles(t,G.blocks+b,tile,mark,b+1,touched);for(int q=0;q<used;q++)list[pos[touched[q]]++]=b;}*offset=ptr;*ids=list;unsigned long long local=total,all=0;MPI_Reduce(&local,&all,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,mpi_comm_level1);if(!G.rank){printf("CPU_GSOC_BLOCK_BUCKETS sector=%d tiles=%d block_visits=%llu full_B_rescans=0 bucket_bytes_local=%zu\n",(int)(t-G.sec),ntile,all,total*sizeof(int)+(size_t)(ntile+1)*sizeof(int));fflush(stdout);}free(mark);free(touched);free(pos);
}
static void gasm(GSec*t,void*H,void*S){
 int tile=genv("OPENMX_CPU_GSOC_TILE",32,1,128),fields=t->dtype==1?2:4,p=G.np;size_t cap=(size_t)t->n*tile*4;
 if(cap>INT_MAX)gdie("GSOC assembly tile exceeds MPI count capacity");double*buf=gmem(cap,sizeof(double)),*sum=G.rank?NULL:gmem(cap,sizeof(double)),*pack=G.rank?NULL:gmem(cap,sizeof(double));
 size_t recvcap=(size_t)t->nr*tile*fields;double*recv=gmem(recvcap,sizeof(double));int*ct=gmem(p,sizeof(int)),*dp=gmem(p,sizeof(int)),*pos=gmem(p,sizeof(int));double stats[6]={0};int*bucket_ptr=NULL,*bucket_ids=NULL;int indexed=genv("OPENMX_CPU_GSOC_BLOCK_BUCKETS",1,0,1);if(indexed)gblock_buckets(t,tile,&bucket_ptr,&bucket_ids);
 for(int j0=0;j0<t->n;j0+=tile){int j1=j0+tile;if(j1>t->n)j1=t->n;size_t nvalues=(size_t)t->n*(j1-j0)*4;memset(buf,0,nvalues*sizeof(double));
  for(int bucket=indexed?bucket_ptr[j0/tile]:0;bucket<(indexed?bucket_ptr[j0/tile+1]:G.nblocks);bucket++){GBlock*v=G.blocks+(indexed?bucket_ids[bucket]:bucket);int a=v->a,b=v->b,ga0=v->ga,gb=v->gb,na=v->na,nb=v->nb;for(int i=0;i<na;i++)for(int j=0;j<nb;j++){int I=G.mp[ga0]+i,J=G.mp[gb]+j;double ov=G.olp[a][b][i][j];gentry(buf,j0,j1,t,I,J,gz(G.nh[0][a][b][i][j],G.im[0][a][b][i][j]),ov);gentry(buf,j0,j1,t,I+G.n,J+G.n,gz(G.nh[1][a][b][i][j],G.im[1][a][b][i][j]),ov);dcomplex ud=gz(G.nh[2][a][b][i][j],G.nh[3][a][b][i][j]+G.im[2][a][b][i][j]);gentry(buf,j0,j1,t,I,J+G.n,ud,0);gentry(buf,j0,j1,t,J+G.n,I,gc(ud),0);}}
  MPI_Reduce(buf,sum,(int)nvalues,MPI_DOUBLE,MPI_SUM,0,mpi_comm_level1);memset(ct,0,(size_t)p*sizeof(int));for(int j=j0;j<j1;j++)for(int i=0;i<t->n;i++){int dest=grank(i,j,t);if(dest<0||dest>=p||ct[dest]>INT_MAX-fields)gdie("assembly rank/count overflow");ct[dest]+=fields;}dp[0]=0;for(int r=1;r<p;r++){if(dp[r-1]>INT_MAX-ct[r-1])gdie("assembly displacement overflow");dp[r]=dp[r-1]+ct[r-1];}memcpy(pos,dp,(size_t)p*sizeof(int));
  if(!G.rank)for(int j=j0;j<j1;j++)for(int i=0;i<t->n;i++){size_t from=((size_t)(j-j0)*t->n+i)*4;double hr=sum[from],hi=sum[from+1],sr=sum[from+2],si=sum[from+3];if(!isfinite(hr)||!isfinite(hi)||!isfinite(sr)||!isfinite(si))gdie("nonfinite assembled H/S");if(fabs(hi)>stats[0])stats[0]=fabs(hi);if(fabs(si)>stats[1])stats[1]=fabs(si);stats[2]+=hr*hr;stats[3]+=hi*hi;stats[4]+=sr*sr;stats[5]+=si*si;int dest=grank(i,j,t),at=pos[dest];if(at<dp[dest]||at>dp[dest]+ct[dest]-fields||(size_t)at+fields>cap)gdie("root assembly pack bounds");pack[at++]=hr;if(fields==4)pack[at++]=hi;pack[at++]=sr;if(fields==4)pack[at++]=si;pos[dest]=at;}
  if((size_t)ct[G.rank]>recvcap)gdie("assembly receive capacity mismatch");MPI_Scatterv(pack,ct,dp,MPI_DOUBLE,recv,ct[G.rank],MPI_DOUBLE,0,mpi_comm_level1);int at=0;for(int j=j0;j<j1;j++)if((j/t->nb)%np_cols2==my_pcol2)for(int i=0;i<t->n;i++)if((i/t->nb)%np_rows2==my_prow2){int il=(i/(t->nb*np_rows2))*t->nb+i%t->nb,jl=(j/(t->nb*np_cols2))*t->nb+j%t->nb;size_t to=(size_t)jl*t->nr+il;if(to>=(size_t)t->nr*t->nc||at>ct[G.rank]-fields)gdie("assembly local unpack bounds");if(fields==2){((double*)H)[to]=recv[at++];((double*)S)[to]=recv[at++];}else{((dcomplex*)H)[to]=gz(recv[at],recv[at+1]);at+=2;((dcomplex*)S)[to]=gz(recv[at],recv[at+1]);at+=2;}}if(at!=ct[G.rank])gdie("assembly local unpack count");
 }
 MPI_Bcast(stats,6,MPI_DOUBLE,0,mpi_comm_level1);if(t->dtype==1){double hr=sqrt(stats[3])/fmax(sqrt(stats[2]),DBL_MIN),sr=sqrt(stats[5])/fmax(sqrt(stats[4]),DBL_MIN);if(!G.rank){printf("CPU_GSOC_REAL_GAUGE sector=%d H_imag_max=%.12e H_imag_rel=%.12e S_imag_max=%.12e S_imag_rel=%.12e projected=%d\n",(int)(t-G.sec),stats[0],hr,stats[1],sr,G.projected);fflush(stdout);}if(!isfinite(hr)||!isfinite(sr)||hr>G.tolerance||sr>G.tolerance)gdie("real sector exceeds declared imaginary covariance tolerance");}
 free(buf);free(sum);free(pack);free(recv);free(ct);free(dp);free(pos);free(bucket_ptr);free(bucket_ids);
}
static void gsolve_sector(GSec*t,double*eval,int rows,int cols,int world){
 int one=1,info,nev=t->n,wanted=nev;size_t count=(size_t)t->nr*t->nc,bytes=t->dtype==1?sizeof(double):sizeof(dcomplex);void*H=gmem(count,bytes),*S=gmem(count,bytes),*Q=gmem(count,bytes);double a0=MPI_Wtime();gasm(t,H,S);double a1=MPI_Wtime(),scale=1;
 if(t->dtype==1){double*h=H,*s=S,*q=Q,alpha=1,beta=0;pdpotrf_("L",&t->n,s,&one,&one,t->desc,&info);if(info)gdie("real sector overlap Cholesky");pdsygst_(&one,"L",&t->n,h,&one,&one,t->desc,s,&one,&one,t->desc,&scale,&info);if(info||!isfinite(scale)||scale<=0)gdie("real generalized transform");
  for(int j=0;j<t->nc;j++){int gj=gsoc_global(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++)if(gsoc_global(i,t->nb,my_prow2,np_rows2)<gj)h[(size_t)j*t->nr+i]=0;}
  pdtran_(&t->n,&t->n,&alpha,h,&one,&one,t->desc,&beta,q,&one,&one,t->desc);for(int j=0;j<t->nc;j++){int gj=gsoc_global(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++)if(gsoc_global(i,t->nb,my_prow2,np_rows2)!=gj)h[(size_t)j*t->nr+i]+=q[(size_t)j*t->nr+i];}
  double a2=MPI_Wtime();elpa_solve_evp_real_2stage_double_impl_(&t->n,&nev,h,&t->nr,eval+t->ofs,q,&t->nr,&t->nb,&t->nc,&rows,&cols,&world);if(nev!=wanted)gdie("real ELPA returned state count changed");double a3=MPI_Wtime();pdtrsm_("L","L","T","N",&t->n,&wanted,&alpha,s,&one,&one,t->desc,q,&one,&one,t->desc);printf("CPU_GSOC_PHASE rank=%d sector=%d dtype=real64 n=%d assemble=%.9f transform=%.9f eigen=%.9f back=%.9f\n",G.rank,(int)(t-G.sec),t->n,a1-a0,a2-a1,a3-a2,MPI_Wtime()-a3);
 }else{dcomplex*h=H,*s=S,*q=Q,alpha=gz(1,0),beta=gz(0,0);pzpotrf_("L",&t->n,s,&one,&one,t->desc,&info);if(info)gdie("complex sector overlap Cholesky");pzhegst_(&one,"L",&t->n,h,&one,&one,t->desc,s,&one,&one,t->desc,&scale,&info);if(info||!isfinite(scale)||scale<=0)gdie("complex generalized transform");for(int j=0;j<t->nc;j++){int gj=gsoc_global(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++){int gi=gsoc_global(i,t->nb,my_prow2,np_rows2);if(gi<gj)h[(size_t)j*t->nr+i]=gz(0,0);else if(gi==gj)h[(size_t)j*t->nr+i].i=0;}}
  pztranc_(&t->n,&t->n,&alpha,h,&one,&one,t->desc,&beta,q,&one,&one,t->desc);for(int j=0;j<t->nc;j++){int gj=gsoc_global(j,t->nb,my_pcol2,np_cols2);for(int i=0;i<t->nr;i++)if(gsoc_global(i,t->nb,my_prow2,np_rows2)!=gj)h[(size_t)j*t->nr+i]=ga(h[(size_t)j*t->nr+i],q[(size_t)j*t->nr+i]);}
  double a2=MPI_Wtime();elpa_solve_evp_complex_2stage_double_impl_(&t->n,&nev,h,&t->nr,eval+t->ofs,q,&t->nr,&t->nb,&t->nc,&rows,&cols,&world);if(nev!=wanted)gdie("complex ELPA returned state count changed");double a3=MPI_Wtime();pztrsm_("L","L","C","N",&t->n,&wanted,&alpha,s,&one,&one,t->desc,q,&one,&one,t->desc);printf("CPU_GSOC_PHASE rank=%d sector=%d dtype=complex128 n=%d assemble=%.9f transform=%.9f eigen=%.9f back=%.9f\n",G.rank,(int)(t-G.sec),t->n,a1-a0,a2-a1,a3-a2,MPI_Wtime()-a3);
 }
 MPI_Bcast(eval+t->ofs,t->n,MPI_DOUBLE,0,mpi_comm_level1);for(int j=0;j<t->n;j++){eval[t->ofs+j]*=scale;if(!isfinite(eval[t->ofs+j]))gdie("nonfinite returned eigenvalue");}if(t->child>=0)memcpy(eval+G.sec[t->child].ofs,eval+t->ofs,(size_t)t->n*sizeof(double));t->Q=Q;free(H);free(S);fflush(stdout);
}
typedef struct{double e;int source;} GEigen;
static int gecmp(const void*aa,const void*bb){const GEigen*a=aa,*b=bb;if(a->e<b->e)return-1;if(a->e>b->e)return 1;return(a->source>b->source)-(a->source<b->source);}
static int gfinalowner(int band){for(int r=0;r<G.np;r++)if(G.is2[r]<=band+1&&band+1<=G.ie2[r])return r;gdie("missing final native band owner");return 0;}
static void greturn(double*eval,int maxn){
 GEigen*sort=gmem(G.ns,sizeof(GEigen));int*inverse=gmem(G.ns,sizeof(int));for(int j=0;j<G.ns;j++){sort[j].e=eval[j];sort[j].source=j;}qsort(sort,G.ns,sizeof(GEigen),gecmp);for(int j=0;j<G.ns;j++){G.ko[j+1]=sort[j].e;inverse[sort[j].source]=j;}int*sorted_sector=gmem(G.ns,sizeof(int));uint64_t expected_hash[48]={0};for(int k=0;k<G.nsec;k++){GSec*t=G.sec+k;for(int row=0;row<t->n;row++){sorted_sector[inverse[t->ofs+row]]=k;expected_hash[k]+=ghash((uint64_t)row);}}int local=G.ie2[G.rank]>=G.is2[G.rank]?G.ie2[G.rank]-G.is2[G.rank]+1:0;memset(G.ev,0,(size_t)local*G.ns*sizeof(dcomplex));int*received=gmem(local,sizeof(int));uint64_t*row_sum=gmem(local,sizeof(uint64_t)),*row_hash=gmem(local,sizeof(uint64_t));int p=G.np;int*ct=gmem(p,sizeof(int)),*rc=gmem(p,sizeof(int)),*dp=gmem(p,sizeof(int)),*rd=gmem(p,sizeof(int)),*pos=gmem(p,sizeof(int));
 for(int k=0;k<G.nsec;k++){GSec*t=G.sec+k;if(t->partner>=0)continue;int choicecount=t->child>=0?2:1;
  /* Small bounded column packets; no global dense eigenvector transpose. */
  for(int b0=0;b0<t->n;b0+=16){int b1=b0+16;if(b1>t->n)b1=t->n;memset(ct,0,(size_t)p*sizeof(int));
   for(int j=0;j<t->nc;j++){int b=gsoc_global(j,t->nb,my_pcol2,np_cols2);if(b<b0||b>=b1)continue;for(int choice=0;choice<choicecount;choice++){GSec*v=choice?G.sec+t->child:t;int band=inverse[v->ofs+b];if(band>=maxn)continue;int dest=gfinalowner(band);for(int i=0;i<t->nr;i++)if(gsoc_global(i,t->nb,my_prow2,np_rows2)<t->n){if(ct[dest]==INT_MAX)gdie("coefficient count overflow");ct[dest]++;}}}
   MPI_Alltoall(ct,1,MPI_INT,rc,1,MPI_INT,mpi_comm_level1);dp[0]=rd[0]=0;for(int r=1;r<p;r++){if(dp[r-1]>INT_MAX-ct[r-1]||rd[r-1]>INT_MAX-rc[r-1])gdie("coefficient displacement overflow");dp[r]=dp[r-1]+ct[r-1];rd[r]=rd[r-1]+rc[r-1];}if(dp[p-1]>INT_MAX-ct[p-1]||rd[p-1]>INT_MAX-rc[p-1])gdie("coefficient total count overflow");int nsend=dp[p-1]+ct[p-1],nrecv=rd[p-1]+rc[p-1];GRecord*send=gmem(nsend,sizeof(GRecord)),*recv=gmem(nrecv,sizeof(GRecord));memcpy(pos,dp,(size_t)p*sizeof(int));
   for(int j=0;j<t->nc;j++){int b=gsoc_global(j,t->nb,my_pcol2,np_cols2);if(b<b0||b>=b1)continue;for(int choice=0;choice<choicecount;choice++){GSec*v=choice?G.sec+t->child:t;int band=inverse[v->ofs+b];if(band>=maxn)continue;int dest=gfinalowner(band);for(int i=0;i<t->nr;i++){int row=gsoc_global(i,t->nb,my_prow2,np_rows2);if(row>=t->n)continue;int at=pos[dest]++;if(at<dp[dest]||at>=dp[dest]+ct[dest])gdie("coefficient pack bounds");GRecord*z=send+at;z->band=band;z->row=v->ofs+row;z->z=t->dtype==1?gz(((double*)t->Q)[(size_t)j*t->nr+i],0):((dcomplex*)t->Q)[(size_t)j*t->nr+i];if(choice)z->z=gc(z->z);}}}
   for(int r=0;r<p;r++){if(ct[r]>INT_MAX/(int)sizeof(GRecord)||rc[r]>INT_MAX/(int)sizeof(GRecord)||dp[r]>INT_MAX/(int)sizeof(GRecord)||rd[r]>INT_MAX/(int)sizeof(GRecord))gdie("coefficient MPI byte capacity");ct[r]*=sizeof(GRecord);rc[r]*=sizeof(GRecord);dp[r]*=sizeof(GRecord);rd[r]*=sizeof(GRecord);}
   MPI_Alltoallv(send,ct,dp,MPI_BYTE,recv,rc,rd,MPI_BYTE,mpi_comm_level1);
   for(int at=0;at<nrecv;at++){GRecord*z=recv+at;int slot=z->band+1-G.is2[G.rank];if(slot<0||slot>=local||z->row<0||z->row>=G.ns)gdie("native coefficient owner mismatch");GSec*logical=G.sec+sorted_sector[z->band];int localrow=z->row-logical->ofs;if(localrow<0||localrow>=logical->n)gdie("returned coefficient sector mismatch");received[slot]++;row_sum[slot]+=(uint64_t)localrow;row_hash[slot]+=ghash((uint64_t)localrow);size_t base=(size_t)slot*G.ns;for(int64_t a=G.cp[z->row];a<G.cp[z->row+1];a++){int row=G.ri[a];G.ev[base+row]=ga(G.ev[base+row],gm(G.cv[a],z->z));}}
   free(send);free(recv);
  }
  free(t->Q);t->Q=NULL;
 }
 for(int slot=0;slot<local;slot++){int band=G.is2[G.rank]+slot-1,k=sorted_sector[band];uint64_t n=G.sec[k].n;if(received[slot]!=G.sec[k].n||row_sum[slot]!=n*(n-1)/2||row_hash[slot]!=expected_hash[k])gdie("missing/duplicate distributed sector coefficient detected");}free(received);free(row_sum);free(row_hash);free(sorted_sector);free(ct);free(rc);free(dp);free(rd);free(pos);free(sort);free(inverse);
 if(!G.rank){printf("CPU_GSOC_NATIVE_RETURN states_total=%d states_requested=%d all_required_states_retained=1 native_band_major_EVec=1 sparse_basis_backtransform=1 global_dense_transpose=0\n",G.ns,maxn);fflush(stdout);}
}
static void gcompact_return(double*eval,int maxn){
 CpuGSOCCompactSector sec[48];memset(sec,0,sizeof(sec));
 for(int k=0;k<G.nsec;k++){GSec*t=G.sec+k;CpuGSOCCompactSector*s=sec+k;s->n=t->n;s->ofs=t->ofs;s->dtype=t->dtype;s->partner=t->partner;s->child=t->child;s->nr=t->nr;s->nc=t->nc;s->nb=t->nb;s->prow=my_prow2;s->pcol=my_pcol2;s->nprows=np_rows2;s->npcols=np_cols2;s->Q=t->Q;}
 CpuGSOC_CompactBuild(mpi_comm_level1,G.ns,G.nsec,sec,G.cp,G.ri,G.cv,eval,maxn,G.ko);
 for(int k=0;k<G.nsec;k++)G.sec[k].Q=sec[k].Q;
 if(!G.rank){printf("CPU_GSOC_NATIVE_RETURN states_total=%d states_requested=%d all_required_states_retained=1 native_band_major_EVec=0 compact_independent_Q=1 sparse_basis_backtransform=on_demand global_dense_transpose=0\n",G.ns,maxn);fflush(stdout);}
}
static int ggetvec(void*ctx,int band,GsocAuditComplex*out,double*e){if(band<1||band>G.ns)return-1;if(CpuGSOC_CompactEnabled()){*e=G.ko[band];return CpuGSOC_CompactGetPhysical(band,(dcomplex*)out);}int dest=gfinalowner(band-1);if(G.rank==dest)memcpy(out,G.ev+(size_t)(band-G.is2[G.rank])*G.ns,(size_t)G.ns*sizeof(dcomplex));MPI_Bcast(out,2*G.ns,MPI_DOUBLE,dest,mpi_comm_level1);*e=G.ko[band];return 0;}
void CpuGSOC_Release(void){for(int k=0;k<G.nsec;k++){free(G.sec[k].Q);G.sec[k].Q=NULL;}free(G.blocks);G.blocks=NULL;free(G.offset);free(G.dtype);if(!gexport_enabled())free(G.mp);free(G.cp);free(G.ri);free(G.cv);free(G.rp);free(G.ci);free(G.rv);G.offset=G.dtype=G.rp=G.ci=G.ri=NULL;if(!gexport_enabled())G.mp=NULL;G.cp=NULL;G.cv=G.rv=NULL;}
void CpuGSOC_Solve(int step,int maxn,double*ko,double*****nh,double*****im,double****overlap,int*is2,int*ie2,dcomplex*sf,dcomplex*h,dcomplex*q,dcomplex*ev){
 (void)sf;(void)h;(void)q;if(!CpuGSOC_Enabled())gdie("disabled GSOC called");if(step!=1||DFTSCF_loop!=1||Solver!=2||SpinP_switch!=3||SO_switch!=1||MD_switch!=0||Cnt_switch)gdie("GSOC requires fixed one-step Gamma cluster SOC and uncontracted basis");if(Hub_U_switch||Constraint_NCS_switch||Zeeman_NCS_switch||Zeeman_NCO_switch)gdie("GSOC v1 excludes Hubbard/constraint/Zeeman spin symmetry extensions");if(Use_of_Collinear_Restart&&(Scf_RestartFromFile!=1||SpinP_switch_RestartFiles!=0||gsoc_restart_read!=1))gdie("GSOC c2n requires actual successful restart reading and a nonmagnetic scalar checkpoint");if(!Use_of_Collinear_Restart&&!genv("OPENMX_CPU_GSOC_ALLOW_FRESH",0,0,1))gdie("GSOC requires collinear c2n checkpoint or explicit reduced fresh-smoke flag");if(!nh||!im||!overlap||!ev||!ko||!is2||!ie2)gdie("missing native SOC arguments");double start=MPI_Wtime();gload();if(!Use_of_Collinear_Restart&&G.ns>512)gdie("fresh GSOC is restricted to <=512 spinors; production requires a validated c2n checkpoint");if(maxn<1||maxn>G.ns)gdie("native MaxN out of full spinor range");G.nh=nh;G.im=im;G.olp=overlap;gblocks_init();G.ko=ko;G.ev=ev;G.is2=is2;G.ie2=ie2;G.projected=genv("OPENMX_CPU_GSOC_PROJECT",0,0,1);G.tolerance=grealenv("OPENMX_CPU_GSOC_COV_TOL",1e-10,1e-14,1e-3);if(!G.projected&&G.tolerance>1e-10)gdie("relaxed covariance requires explicit GSOC_PROJECT=1");
 GsocAuditOperator op={G.ns,G.nsec,G.offset,G.dtype,NULL,gapply,gforward,ginverse};GsocAuditGateResult gate;int probes=genv("OPENMX_CPU_GSOC_PROBES",4,1,16);if(gsoc_audit_operator_gate(&op,probes,UINT64_C(20261028),&gate)||gate.nonfinite)gdie("raw covariance operator measurement failed");if(!G.rank){printf("CPU_GSOC_RAW_GATE sampled_probes=%d entrywise_proof=0 projected=%d tolerance=%.12e H_comm_abs=%.12e H_comm_rel=%.12e S_comm_abs=%.12e S_comm_rel=%.12e H_TR_abs=%.12e H_TR_rel=%.12e S_TR_abs=%.12e S_TR_rel=%.12e H_real_abs=%.12e H_real_rel=%.12e S_real_abs=%.12e S_real_rel=%.12e\n",probes,G.projected,G.tolerance,gate.comm_abs[0],gate.comm_rel[0],gate.comm_abs[1],gate.comm_rel[1],gate.tr_abs[0],gate.tr_rel[0],gate.tr_abs[1],gate.tr_rel[1],gate.real_abs[0],gate.real_rel[0],gate.real_abs[1],gate.real_rel[1]);fflush(stdout);}for(int field=0;field<2;field++)if(gate.comm_rel[field]>G.tolerance||gate.tr_rel[field]>G.tolerance||gate.real_rel[field]>G.tolerance)gdie("raw H/S covariance exceeds declared threshold");GsocAuditHermitianResult herm;if(gsoc_audit_hermiticity(&op,probes,UINT64_C(20261029),&herm)||herm.nonfinite)gdie("raw Hermiticity operator measurement failed");if(!G.rank){printf("CPU_GSOC_RAW_HERMITICITY sampled_probes=%d H_abs=%.12e H_rel=%.12e S_abs=%.12e S_rel=%.12e\n",probes,herm.absolute[0],herm.relative[0],herm.absolute[1],herm.relative[1]);fflush(stdout);}for(int field=0;field<2;field++)if(herm.relative[field]>G.tolerance)gdie("raw H/S Hermiticity exceeds declared threshold");
 MPI_Comm rowscomm,colscomm;MPI_Comm_split(mpi_comm_level1,my_pcol2,my_prow2,&rowscomm);MPI_Comm_split(mpi_comm_level1,my_prow2,my_pcol2,&colscomm);int rows=MPI_Comm_c2f(rowscomm),cols=MPI_Comm_c2f(colscomm),world=MPI_Comm_c2f(mpi_comm_level1);double*eval=gmem(G.ns,sizeof(double));for(int k=0;k<G.nsec;k++)if(G.sec[k].partner<0)gsolve_sector(G.sec+k,eval,rows,cols,world);if(CpuGSOC_CompactEnabled())gcompact_return(eval,maxn);else greturn(eval,maxn);free(eval);MPI_Comm_free(&rowscomm);MPI_Comm_free(&colscomm);
 double electrons=-system_charge;for(int atom=1;atom<=atomnum;atom++)electrons+=Spe_Core_Charge[WhatSpecies[atom]];if(!isfinite(electrons)||electrons<0||electrons>G.ns)gdie("invalid physical electron count for eigenpair sampling");int homo=(int)ceil(electrons);int samples[8],count=0,initial[8]={1,2,(maxn+1)/2,(maxn+1)/2+1,maxn-1,maxn,homo,homo+1};for(int i=0;i<8;i++){int band=initial[i],duplicate=0;if(band<1||band>maxn)continue;for(int j=0;j<count;j++)if(samples[j]==band)duplicate=1;if(!duplicate)samples[count++]=band;}GsocAuditEigenResult result;if(gsoc_audit_eigenpairs(&op,ggetvec,samples,count,&result)||result.nonfinite)gdie("native physical eigenpair audit failed");if(!G.rank){printf("CPU_GSOC_NATIVE_AUDIT sampled_bands=%d homo_sample=%d lumo_sample=%d electrons=%.9f raw_residual_max=%.12e raw_residual_relative=%.12e raw_S_orth=%.12e full_state_drop=0 true_SOC_energy_forces_native=1 seconds=%.9f\n",count,homo,homo+1,electrons,result.raw_abs,result.raw_relative,result.Sorth,MPI_Wtime()-start);fflush(stdout);}double allowance=G.projected?fmax(1e-8,20*G.tolerance):1e-8;if(result.raw_relative>allowance||result.Sorth>allowance)gdie("native spinor residual/Sorth exceeds declared projection allowance");gexport_reference();CpuSOCNative_MarkSolved(step);CpuGSOC_Release();
}
