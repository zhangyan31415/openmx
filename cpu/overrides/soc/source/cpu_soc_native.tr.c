/* Distributed finite-range real-space C3/D3. Original local blocks plus one
 * accumulator and one send/receive group only: no replicated dense N*N and
 * no six full operator copies. Active H' = D_i H D_j^T, R'=W R+s_j-s_i. */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>
#include <mkl.h>
#include "openmx_common.h"
#include "cpu_soc_native.h"
#include "cpu_gsoc.h"
#include "cpu_c3.h"
#include "soc_spin.h"
#include "cpu_soc_tr_guard.h"
#ifndef CPU_BAND_USES_SPARSE_HS
#define CPU_BAND_USES_SPARSE_HS() 0
#endif

typedef struct { int a,b,r[3],h,na,nb,owner; } RSBlock;
typedef struct {
  int *sc,*sd,*rc,*rd,*order,*target,*recv_target,*recv_offset;
  int *bc,*bd,*brc,*brd;
} RSRoute;
#define RS_MAX_OPS 192
#define RS_REVERSE RS_MAX_OPS
static int rs_generic=0;
static int rs_ng=-1;
static int rs_mode=-1,rs_ready=0,rs_s_done=0,rs_rank,rs_np,rs_nat,rs_n,rs_ns,rs_nb,rs_nlocal,rs_maxorb;
static int *rs_type,*rs_norb,*rs_local_id,*rs_g2l,*rs_offset,*rs_hash,rs_hashsize;
static int *rs_perm[RS_MAX_OPS],*rs_img[RS_MAX_OPS],rs_w[RS_MAX_OPS][9];
static double rs_cell[9],*rs_xyz,*rs_d[RS_MAX_OPS],rs_group_error;
static size_t rs_v,rs_bytes=0;
static uint64_t rs_topology,rs_map_fnv=UINT64_C(1469598103934665603);
static RSBlock *rs_block;
static RSRoute rs_route[RS_MAX_OPS+1];
static double *rs_orig,*rs_sum,*rs_send,*rs_recv,*rs_tmp,*rs_mat;
static double ****rs_overlap;
static int rs_last_step=-1,rs_solved=0,rs_exported=0,rs_restart_ok=0;
static double rs_spin[RS_MAX_OPS][8],rs_spin_t[RS_MAX_OPS][64];
static double *****rs_hptr,*****rs_iptr;
static uint64_t rs_solve_hash;
static int rs_band_prepared=0,rs_band_export_handoff=0;
static uint64_t rs_band_hash;
static double *rs_soc_raw[8],*rs_soc_projected[8];

static void rs_die(const char *msg) {
  int rank=0; MPI_Comm_rank(mpi_comm_level1,&rank);
  fprintf(stderr,"CPU_SOC_RS_REJECT rank=%d reason=%s\n",rank,msg);fflush(stderr);
  MPI_Abort(mpi_comm_level1,96); exit(96);
}
static void *rs_alloc(size_t n,size_t width) {
  void *p=calloc(n?n:1,width);if(!p)rs_die("allocation failed");rs_bytes+=(n?n:1)*width;return p;
}
static void rs_read(FILE *f,void *ptr,size_t width,size_t count) {
  size_t i,total=width*count;unsigned char *b=ptr;
  if(fread(ptr,width,count,f)!=count)rs_die("truncated metadata");
  for(i=0;i<total;i++){rs_map_fnv^=b[i];rs_map_fnv*=UINT64_C(1099511628211);}
}
static uint64_t rs_mix(uint64_t h,uint32_t x){h^=x;return h*UINT64_C(1099511628211);}
static uint64_t rs_keyhash(int a,int b,const int *r) {
  uint64_t h=UINT64_C(1469598103934665603);int j;
  h=rs_mix(h,a);h=rs_mix(h,b);for(j=0;j<3;j++)h=rs_mix(h,(uint32_t)r[j]);return h;
}
static int rs_find(int a,int b,const int *r) {
  unsigned int k=(unsigned int)rs_keyhash(a,b,r)&(rs_hashsize-1);
  for(;;){int id=rs_hash[k];RSBlock *x;
    if(id<0)return -1;x=&rs_block[id];
    if(x->a==a&&x->b==b&&x->r[0]==r[0]&&x->r[1]==r[1]&&x->r[2]==r[2])return id;
    k=(k+1)&(rs_hashsize-1);
  }
}
static uint64_t rs_current_topology(void) {
  uint64_t h=UINT64_C(1469598103934665603);int a,j,k;
  h=rs_mix(h,Matomnum);for(a=1;a<=Matomnum;a++)h=rs_mix(h,M2G[a]);
  for(a=1;a<=atomnum;a++){
    h=rs_mix(h,WhatSpecies[a]);h=rs_mix(h,G2ID[a]);h=rs_mix(h,FNAN[a]);
    for(j=0;j<=FNAN[a];j++){h=rs_mix(h,natn[a][j]);for(k=1;k<=3;k++)h=rs_mix(h,(uint32_t)atv_ijk[ncn[a][j]][k]);}
  }return h;
}
int CpuSOCNative_Enabled(void) {
 const char *group=getenv("OPENMX_CPU_SOC_NATIVE_GROUP");
 const char *legacy=getenv("OPENMX_CPU_SOC_NATIVE_D3");
 const char *d3=getenv("OPENMX_CPU_D3");
 int ng;
 const char *generic=getenv("OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE");
 if(generic&&strcmp(generic,"0")&&strcmp(generic,"1"))rs_die("GSOC_EXPORT_SYMMETRIZE must be 0 or 1");
 if(generic&&!strcmp(generic,"1")){
  const char*project=getenv("OPENMX_CPU_GSOC_PROJECT");
  if(!CpuGSOC_Enabled()||!project||strcmp(project,"1"))rs_die("generic symmetrized export requires explicit GSOC_PROJECT=1");
  if(group||(legacy&&strcmp(legacy,"0"))||(d3&&strcmp(d3,"0")))rs_die("generic export conflicts with legacy solver selectors");
  const char*name=getenv("OPENMX_CPU_SOC_NATIVE_MAP");int32_t hd[8];
  if(!name||name[0]!='/')rs_die("generic export requires an absolute native map");
  FILE*f=fopen(name,"rb");if(!f||fread(hd,4,8,f)!=8)rs_die("generic export map header unavailable");fclose(f);
  if(hd[0]!=0x52534433||hd[1]!=3||hd[4]<1||hd[4]>RS_MAX_OPS)rs_die("generic export requires version3 bounded SU2 map");
  ng=hd[4];rs_generic=1;
 }else{
 if(legacy&&strcmp(legacy,"0")&&strcmp(legacy,"1"))rs_die("OPENMX_CPU_SOC_NATIVE_D3 must be 0 or 1");
 if(group){
  if(!strcmp(group,"C3"))ng=3;
  else if(!strcmp(group,"D3"))ng=6;
  else rs_die("OPENMX_CPU_SOC_NATIVE_GROUP must be C3 or D3");
  if(legacy&&(!strcmp(legacy,"0")||ng!=6))rs_die("native GROUP and legacy D3 selectors conflict");
 }else{
  if(!legacy||!strcmp(legacy,"0"))return 0;
  ng=6;
 }
 if(d3&&strcmp(d3,"0")&&strcmp(d3,"1"))rs_die("OPENMX_CPU_D3 must be 0 or 1 for native export");
 if((ng==6)!=!!(d3&&!strcmp(d3,"1")))rs_die("native group and solver D3 selector disagree");
 }
 if(rs_ng!=-1&&rs_ng!=ng)rs_die("native group changed after activation");
 rs_ng=ng;
 if(SpinP_switch!=3||SO_switch!=1||Cnt_switch||Hub_U_switch||Solver!=2||MD_switch!=0||SCF_Control_Temp||DFTSCF_loop!=1)rs_die("requires uncontracted fixed-geometry Cluster SOC, Hubbard off, exactly one step");
 if((Band_disp_switch&&!CPU_BAND_USES_SPARSE_HS())||MO_fileout||Dos_fileout||unfold_electronic_band||Wannier_Func_Calc||Wannier90_fileout||CWF_fileout_flag||pop_anal_aow_flag)rs_die("auxiliary raw-matrix postprocessing must be disabled; use native scfout H/S");
 if(Constraint_NCS_switch||Zeeman_NCS_switch||Zeeman_NCO_switch)rs_die("magnetic constraints/fields not supported");
 if((rs_generic?CpuC3_Mode()!=0:CpuC3_Mode()!=3)||(!HS_fileout&&!CPU_BAND_USES_SPARSE_HS()))rs_die("native-export requires a projected solver and H/S output or targeted bands");
 rs_mode=rs_generic?2:1;return 1;
}
static void rs_geometry_guard(void) {
  int a,j;double error=0;
  if(atomnum!=rs_nat)rs_die("atom count changed; regenerate metadata");
  for(a=0;a<3;a++)for(j=0;j<3;j++){double v=tv[a+1][j+1]*BohrR;if(!isfinite(v)||!isfinite(rs_cell[3*a+j]))rs_die("nonfinite cell geometry");error=fmax(error,fabs(v-rs_cell[3*a+j]));}
  for(a=0;a<rs_nat;a++)for(j=0;j<3;j++){double v=Gxyz[a+1][j+1]*BohrR;if(!isfinite(v)||!isfinite(rs_xyz[3*a+j]))rs_die("nonfinite atom geometry");error=fmax(error,fabs(v-rs_xyz[3*a+j]));}
  if(error>2e-7)rs_die("geometry differs from metadata (2e-7 Ang bound); regenerate map");
  if(rs_current_topology()!=rs_topology)rs_die("basis/distribution/neighbor support changed; fixed-S cache invalid");
}
static void rs_phase(int step,const char *label,double elapsed) {
  double lo,hi;unsigned long long local=(unsigned long long)rs_bytes,maxbytes,sumbytes;
  MPI_Reduce(&elapsed,&lo,1,MPI_DOUBLE,MPI_MIN,Host_ID,mpi_comm_level1);
  MPI_Reduce(&elapsed,&hi,1,MPI_DOUBLE,MPI_MAX,Host_ID,mpi_comm_level1);
  MPI_Reduce(&local,&maxbytes,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,Host_ID,mpi_comm_level1);
  MPI_Reduce(&local,&sumbytes,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,Host_ID,mpi_comm_level1);
  if(rs_rank==Host_ID){printf("CPU_SOC_RS_PHASE step=%d stage=%s wall_min=%.9f wall_max=%.9f spread=%.9f extra_bytes_max=%llu extra_bytes_sum=%llu\n",step,label,lo,hi,hi-lo,maxbytes,sumbytes);fflush(stdout);}
}
static double *rs_rotation(int g,int species) {
  int s;size_t o=0;for(s=0;s<species;s++)o+=(size_t)rs_norb[s]*rs_norb[s];return rs_d[g]+o;
}
static void rs_make_route(int g) {
  RSRoute *t=&rs_route[g];int q,id,i,j,d,*cursor,*sent,*seen;
  t->sc=rs_alloc(rs_np,sizeof(int));t->sd=rs_alloc(rs_np,sizeof(int));t->rc=rs_alloc(rs_np,sizeof(int));t->rd=rs_alloc(rs_np,sizeof(int));
  t->bc=rs_alloc(rs_np,sizeof(int));t->bd=rs_alloc(rs_np,sizeof(int));t->brc=rs_alloc(rs_np,sizeof(int));t->brd=rs_alloc(rs_np,sizeof(int));
  t->order=rs_alloc(rs_nlocal,sizeof(int));t->target=rs_alloc(rs_nlocal,sizeof(int));t->recv_target=rs_alloc(rs_nlocal,sizeof(int));t->recv_offset=rs_alloc(rs_nlocal,sizeof(int));
  for(q=0;q<rs_nlocal;q++){
    RSBlock *b=&rs_block[rs_local_id[q]];int a,bb,r[3];
    if(g==RS_REVERSE){a=b->b;bb=b->a;for(i=0;i<3;i++)r[i]=-b->r[i];}
    else{a=rs_perm[g][b->a-1]+1;bb=rs_perm[g][b->b-1]+1;
      for(i=0;i<3;i++){r[i]=rs_img[g][3*(b->b-1)+i]-rs_img[g][3*(b->a-1)+i];for(j=0;j<3;j++)r[i]+=rs_w[g][3*i+j]*b->r[j];}}
    id=rs_find(a,bb,r);
    if(id<0){char msg[300];snprintf(msg,sizeof(msg),"missing %s closure g=%d source=(%d,%d,%d,%d,%d) target=(%d,%d,%d,%d,%d)",g==RS_REVERSE?"Hermitian reverse":"group neighbor",g,b->a,b->b,b->r[0],b->r[1],b->r[2],a,bb,r[0],r[1],r[2]);rs_die(msg);}
    if((g<rs_ng&&(rs_block[id].na!=b->na||rs_block[id].nb!=b->nb))||(g==RS_REVERSE&&(rs_block[id].na!=b->nb||rs_block[id].nb!=b->na)))rs_die("orbital shape not preserved by map");
    t->target[q]=id;d=rs_block[id].owner;t->bc[d]++;t->sc[d]+=b->na*b->nb;
  }
  for(i=1;i<rs_np;i++){t->sd[i]=t->sd[i-1]+t->sc[i-1];t->bd[i]=t->bd[i-1]+t->bc[i-1];}
  cursor=calloc(rs_np,sizeof(int));sent=calloc(rs_nlocal?rs_nlocal:1,sizeof(int));if(!cursor||!sent)rs_die("route scratch allocation");
  for(q=0;q<rs_nlocal;q++){d=rs_block[t->target[q]].owner;i=t->bd[d]+cursor[d]++;t->order[i]=q;sent[i]=t->target[q];}
  MPI_Alltoall(t->sc,1,MPI_INT,t->rc,1,MPI_INT,mpi_comm_level1);MPI_Alltoall(t->bc,1,MPI_INT,t->brc,1,MPI_INT,mpi_comm_level1);
  for(i=1;i<rs_np;i++){t->rd[i]=t->rd[i-1]+t->rc[i-1];t->brd[i]=t->brd[i-1]+t->brc[i-1];}
  if(t->rd[rs_np-1]+t->rc[rs_np-1]!=(int)rs_v||t->brd[rs_np-1]+t->brc[rs_np-1]!=rs_nlocal)rs_die("map is not bijective on local finite-range support");
  MPI_Alltoallv(sent,t->bc,t->bd,MPI_INT,t->recv_target,t->brc,t->brd,MPI_INT,mpi_comm_level1);
  seen=calloc(rs_nlocal?rs_nlocal:1,sizeof(int));j=0;
  for(i=0;i<rs_nlocal;i++){id=t->recv_target[i];if(id<0||id>=rs_nb||rs_g2l[id]<0)rs_die("received route targets another owner");q=rs_g2l[id];if(seen[q]++)rs_die("duplicate receive target");t->recv_offset[i]=j;j+=rs_block[id].na*rs_block[id].nb;}
  if(j!=(int)rs_v)rs_die("receive block/element order mismatch");free(cursor);free(sent);free(seen);
}
static void rs_validate_metadata(void) {
  int g,h,k,a,i,j,l,s,*seen=calloc(rs_nat,sizeof(int));
  if(!seen)rs_die("metadata validation scratch");
  for(i=0;i<9;i++)if(rs_w[0][i]!=(i/3==i%3))rs_die("operation0 W is not identity");
  for(a=0;a<rs_nat;a++){if(rs_perm[0][a]!=a)rs_die("operation0 permutation not identity");for(i=0;i<3;i++)if(rs_img[0][3*a+i])rs_die("operation0 image not zero");}
  for(s=0;s<rs_ns;s++){int n=rs_norb[s];double *d=rs_rotation(0,s);for(i=0;i<n;i++)for(j=0;j<n;j++)if(fabs(d[i*n+j]-(i==j))>1e-12)rs_die("operation0 orbital matrix not identity");}
  for(g=0;g<rs_ng;g++){memset(seen,0,rs_nat*sizeof(int));for(a=0;a<rs_nat;a++)if(seen[rs_perm[g][a]]++)rs_die("metadata permutation not bijective");}
  for(g=0;g<rs_ng;g++)for(h=0;h<rs_ng;h++){
    int target=-1;
    for(k=0;k<rs_ng;k++){
      int match=1,shift[3]={0,0,0};
      for(i=0;i<3;i++)for(j=0;j<3;j++){int w=0;for(l=0;l<3;l++)w+=rs_w[g][3*i+l]*rs_w[h][3*l+j];if(w!=rs_w[k][3*i+j])match=0;}
      for(a=0;a<rs_nat&&match;a++){
        if(rs_perm[g][rs_perm[h][a]]!=rs_perm[k][a])match=0;
        for(i=0;i<3;i++){int im=rs_img[g][3*rs_perm[h][a]+i]-rs_img[k][3*a+i];for(j=0;j<3;j++)im+=rs_w[g][3*i+j]*rs_img[h][3*a+j];if(a==0)shift[i]=im;else if(im!=shift[i])match=0;}
      }
      if(match){if(target>=0)rs_die("duplicate metadata group operations");target=k;}
    }
    if(target<0)rs_die("integer/permutation/relative-image group composition fails");
    for(i=0;i<8;i++)for(j=0;j<8;j++){double z=0;for(l=0;l<8;l++)z+=rs_spin_t[g][8*i+l]*rs_spin_t[h][8*l+j];if(fabs(z-rs_spin_t[target][8*i+j])>1e-11)rs_die("SU2 adjoint group composition fails");}

    for(s=0;s<rs_ns;s++){
      int n=rs_norb[s];double *x=rs_rotation(g,s),*y=rs_rotation(h,s),*z=rs_rotation(target,s);
      for(i=0;i<n;i++)for(j=0;j<n;j++){double v=0;for(l=0;l<n;l++)v+=x[i*n+l]*y[l*n+j];rs_group_error=fmax(rs_group_error,fabs(v-z[i*n+j]));}
    }
  }
  free(seen);if(rs_group_error>1e-11)rs_die("orbital group composition is not exact selected group");
}
static void rs_setup(void) {
  const char *name;FILE *f;int32_t h[8];int g,a,j,k,b,mc,q,id;size_t nd=0;double start=MPI_Wtime(),orth=0;
  if(rs_ready){rs_geometry_guard();return;}
  MPI_Comm_rank(mpi_comm_level1,&rs_rank);MPI_Comm_size(mpi_comm_level1,&rs_np);
  name=getenv("OPENMX_CPU_SOC_NATIVE_MAP");if(!name||name[0]!='/')rs_die("absolute OPENMX_CPU_SOC_NATIVE_MAP required");
  f=fopen(name,"rb");if(!f)rs_die("cannot open metadata");rs_read(f,h,4,8);
  if(h[0]!=0x52534433||(h[1]!=2&&h[1]!=3)||h[2]!=atomnum||h[4]!=rs_ng||h[5]<1||h[5]>SpeciesNum||h[6]!=0||h[7]!=0)rs_die("metadata header invalid");
  rs_nat=h[2];rs_n=h[3];rs_ns=h[5];rs_read(f,rs_cell,8,9);
  rs_xyz=rs_alloc(3*rs_nat,sizeof(double));rs_type=rs_alloc(rs_nat,sizeof(int));rs_norb=rs_alloc(rs_ns,sizeof(int));
  rs_read(f,rs_xyz,8,3*rs_nat);rs_read(f,rs_type,4,rs_nat);rs_read(f,rs_norb,4,rs_ns);
  for(a=0;a<rs_nat;a++)if(rs_type[a]<0||rs_type[a]>=rs_ns)rs_die("invalid metadata species index");
  for(a=0;a<rs_ns;a++){if(rs_norb[a]<1||rs_norb[a]>128)rs_die("unsupported orbital block size");nd+=(size_t)rs_norb[a]*rs_norb[a];if(rs_norb[a]>rs_maxorb)rs_maxorb=rs_norb[a];}
  for(g=0;g<rs_ng;g++){
    rs_perm[g]=rs_alloc(rs_nat,sizeof(int));rs_img[g]=rs_alloc(3*rs_nat,sizeof(int));rs_d[g]=rs_alloc(nd,sizeof(double));
    rs_read(f,rs_w[g],4,9);rs_read(f,rs_perm[g],4,rs_nat);rs_read(f,rs_img[g],4,3*rs_nat);rs_read(f,rs_d[g],8,nd);
    {size_t z;for(z=0;z<nd;z++)if(!isfinite(rs_d[g][z]))rs_die("nonfinite orbital metadata");}
    for(a=0;a<rs_nat;a++)if(rs_perm[g][a]<0||rs_perm[g][a]>=rs_nat||rs_type[a]!=rs_type[rs_perm[g][a]])rs_die("invalid species permutation");
    for(a=0;a<rs_ns;a++){double *d=rs_rotation(g,a);int n=rs_norb[a];for(j=0;j<n;j++)for(k=0;k<n;k++){double v=0;for(b=0;b<n;b++)v+=d[j*n+b]*d[k*n+b];orth=fmax(orth,fabs(v-(j==k)));}}
  }
  rs_read(f,rs_spin,8,8*rs_ng);
  for(g=0;g<rs_ng;g++)if(SocSpin_Transform(rs_spin[g],rs_spin_t[g]))rs_die("nonunitary SU2 metadata");
  for(a=0;a<8;a++)for(j=0;j<8;j++)if(fabs(rs_spin_t[0][8*a+j]-(a==j))>1e-12)rs_die("spin operation0 not identity");
  if(fgetc(f)!=EOF)rs_die("metadata has trailing bytes");fclose(f);
  if(orth>1e-11)rs_die("orbital metadata is not exact orthogonal group; regenerate nominal-generator map");
  rs_validate_metadata();
  b=0;for(a=1;a<=rs_nat;a++){if(rs_type[a-1]<0||rs_type[a-1]>=rs_ns||Spe_Total_CNO[WhatSpecies[a]]!=rs_norb[rs_type[a-1]])rs_die("native basis differs from metadata");b+=Spe_Total_CNO[WhatSpecies[a]];rs_nb+=FNAN[a]+1;}
  if(b!=rs_n)rs_die("native basis dimension mismatch");
  rs_block=rs_alloc(rs_nb,sizeof(RSBlock));rs_g2l=rs_alloc(rs_nb,sizeof(int));for(a=0;a<rs_nb;a++)rs_g2l[a]=-1;
  id=0;for(a=1;a<=rs_nat;a++)for(j=0;j<=FNAN[a];j++){
    RSBlock *x=&rs_block[id++];x->a=a;x->b=natn[a][j];x->h=j;x->na=Spe_Total_CNO[WhatSpecies[a]];x->nb=Spe_Total_CNO[WhatSpecies[x->b]];x->owner=G2ID[a];for(k=0;k<3;k++)x->r[k]=atv_ijk[ncn[a][j]][k+1];
  }
  rs_hashsize=1;while(rs_hashsize<2*rs_nb)rs_hashsize*=2;rs_hash=rs_alloc(rs_hashsize,sizeof(int));for(a=0;a<rs_hashsize;a++)rs_hash[a]=-1;
  for(id=0;id<rs_nb;id++){RSBlock *x=&rs_block[id];unsigned int pos=(unsigned int)rs_keyhash(x->a,x->b,x->r)&(rs_hashsize-1);if(rs_find(x->a,x->b,x->r)>=0)rs_die("duplicate native neighbor key");while(rs_hash[pos]>=0)pos=(pos+1)&(rs_hashsize-1);rs_hash[pos]=id;}
  for(mc=1;mc<=Matomnum;mc++)rs_nlocal+=FNAN[M2G[mc]]+1;
  rs_local_id=rs_alloc(rs_nlocal,sizeof(int));rs_offset=rs_alloc(rs_nlocal+1,sizeof(int));q=0;
  for(mc=1;mc<=Matomnum;mc++){a=M2G[mc];for(j=0;j<=FNAN[a];j++){int rr[3];for(k=0;k<3;k++)rr[k]=atv_ijk[ncn[a][j]][k+1];id=rs_find(a,natn[a][j],rr);if(id<0)rs_die("owned neighbor not found in global hash");rs_local_id[q]=id;rs_g2l[id]=q;rs_offset[q]=(int)rs_v;rs_v+=(size_t)rs_block[id].na*rs_block[id].nb;if(rs_v>INT_MAX)rs_die("local MPI count overflow");q++;}}
  rs_offset[rs_nlocal]=(int)rs_v;
  for(g=0;g<rs_ng;g++)rs_make_route(g);rs_make_route(RS_REVERSE);
  rs_orig=rs_alloc(rs_v,sizeof(double));rs_sum=rs_alloc(rs_v,sizeof(double));rs_send=rs_alloc(rs_v,sizeof(double));rs_recv=rs_alloc(rs_v,sizeof(double));
  rs_tmp=rs_alloc((size_t)rs_maxorb*rs_maxorb,sizeof(double));rs_mat=rs_alloc((size_t)rs_maxorb*rs_maxorb,sizeof(double));
  rs_topology=rs_current_topology();rs_ready=1;rs_geometry_guard();
  if(rs_rank==Host_ID){printf("CPU_SOC_RS_SETUP mode=%s group=%s operations=%d natoms=%d N=%d global_blocks=%d local_blocks_host=%d map=%s map_fnv64=%016llx orbital_orth_error=%.3e orbital_group_error=%.3e storage=4local_fields_plus_maps ghost_policy=owned_primary_H_S HubU=off contraction=off\n",rs_generic?"generic_export":"export",rs_generic?"GSOC":(rs_ng==6?"D3":"C3"),rs_ng,rs_nat,rs_n,rs_nb,rs_nlocal,name,(unsigned long long)rs_map_fnv,orth,rs_group_error);fflush(stdout);}
  rs_phase(0,"setup",MPI_Wtime()-start);
}
static void rs_pack(double ****x,double *packed,int write) {
  int q,i,j;for(q=0;q<rs_nlocal;q++){RSBlock *b=&rs_block[rs_local_id[q]];int mc=F_G2M[b->a],base=rs_offset[q];for(i=0;i<b->na;i++)for(j=0;j<b->nb;j++){if(write)x[mc][b->h][i][j]=packed[base+i*b->nb+j];else packed[base+i*b->nb+j]=x[mc][b->h][i][j];}}
}
static void rs_exchange(const double *source,int g,double *destination,double scale) {
  RSRoute *t=&rs_route[g];int order,q,id,i,j,pos=0;
  for(order=0;order<rs_nlocal;order++){
    q=t->order[order];RSBlock *b=&rs_block[rs_local_id[q]];const double *x=source+rs_offset[q];
    if(g==RS_REVERSE){for(i=0;i<b->nb;i++)for(j=0;j<b->na;j++)rs_send[pos++]=x[j*b->nb+i];}
    else{
      double *left=rs_rotation(g,rs_type[b->a-1]),*right=rs_rotation(g,rs_type[b->b-1]);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,b->na,b->nb,b->na,1.,left,b->na,x,b->nb,0.,rs_tmp,b->nb);
      cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,b->na,b->nb,b->nb,1.,rs_tmp,b->nb,right,b->nb,0.,rs_mat,b->nb);
      memcpy(rs_send+pos,rs_mat,sizeof(double)*(size_t)b->na*b->nb);pos+=b->na*b->nb;
    }
  }
  MPI_Alltoallv(rs_send,t->sc,t->sd,MPI_DOUBLE,rs_recv,t->rc,t->rd,MPI_DOUBLE,mpi_comm_level1);
  for(order=0;order<rs_nlocal;order++){id=t->recv_target[order];q=rs_g2l[id];RSBlock *b=&rs_block[id];pos=t->recv_offset[order];for(i=0;i<b->na*b->nb;i++)destination[rs_offset[q]+i]+=scale*rs_recv[pos+i];}
}
static void rs_project(double ****x,int step,const char *channel,int write) {
  int g;size_t i;double start=MPI_Wtime();rs_pack(x,rs_orig,0);
  for(i=0;i<rs_v;i++)rs_sum[i]=rs_orig[i]/(double)rs_ng; /* metadata identity is checked by generator */
  for(g=1;g<rs_ng;g++)rs_exchange(rs_orig,g,rs_sum,1./(double)rs_ng);
  for(i=0;i<rs_v;i++)rs_sum[i]*=.5;
  /* Reverse Hermitization commutes with the exact selected group average. */
  rs_exchange(rs_sum,RS_REVERSE,rs_sum,1.); /* reads the half-average, accumulates its transpose */
  if(write)rs_pack(x,rs_sum,1);
  rs_phase(step,channel,MPI_Wtime()-start);
}

/* Check exactly the fields consumed by the solver/export, on every owner.
 * No field is changed until after Force/Total_Energy and auxiliary paths are
 * rejected. Raw array identity alone would miss an intervening mutation. */
static uint64_t rs_field_hash(double *****h,double *****im,double ****olp){
 uint64_t hash=UINT64_C(1469598103934665603);int q,i,j,k;
 for(q=0;q<rs_nlocal;q++){RSBlock *b=&rs_block[rs_local_id[q]];int mc=F_G2M[b->a];
  for(i=0;i<b->na;i++)for(j=0;j<b->nb;j++){
   double vals[8]={h[0][mc][b->h][i][j],h[1][mc][b->h][i][j],h[2][mc][b->h][i][j],h[3][mc][b->h][i][j],im[0][mc][b->h][i][j],im[1][mc][b->h][i][j],im[2][mc][b->h][i][j],olp[mc][b->h][i][j]};
   for(k=0;k<8;k++){unsigned char *bytes=(unsigned char*)(vals+k);int z;if(!isfinite(vals[k]))rs_die("nonfinite raw operator field");for(z=0;z<8;z++){hash^=bytes[z];hash*=UINT64_C(1099511628211);}}
  }
 }return hash;
}
/* Actual raw finite-range H/S check for T=L K, L=[[0,I],[-I,0]].
 * Uses the existing owner buffers and exact Hermitian reverse route RS_REVERSE. No
 * spatial C2 map and no dense matrix are needed. Arrays remain unchanged. */
static int rs_tr_passed=0;
int CpuSOCNative_TRRequested(void){
 const char*s=getenv("OPENMX_CPU_SOC_TR");
 if(!s||!strcmp(s,"0"))return 0;
 return !strcmp(s,"1")?1:-1;
}
int CpuSOCNative_TRPassed(void){return rs_tr_passed;}
void CpuSOCNative_TRValidateConfig(void){
 int requested=CpuSOCNative_TRRequested(),lo,hi,valid,all_valid;
 MPI_Allreduce(&requested,&lo,1,MPI_INT,MPI_MIN,mpi_comm_level1);
 MPI_Allreduce(&requested,&hi,1,MPI_INT,MPI_MAX,mpi_comm_level1);
 if(lo<0||lo!=hi)rs_die("OPENMX_CPU_SOC_TR must be collective-consistent0/1");
 if(!hi)return;
 valid=Solver==2&&SpinP_switch==3&&SO_switch==1&&Cnt_switch==0&&
   MD_switch==0&&DFTSCF_loop==1&&SCF_Control_Temp==0&&CpuC3_Mode()==3&&
   Use_of_Collinear_Restart==1&&Scf_RestartFromFile==1&&SpinP_switch_RestartFiles==0&&
   Hub_U_switch==0&&Constraint_NCS_switch==0&&Zeeman_NCS_switch==0&&Zeeman_NCO_switch==0&&
   core_hole_state_flag==0&&xmcd_calc==0&&level_fileout<2&&
   !CDDF_on&&!LNAO_calc_flag&&!LNBO_calc_flag&&!NBO_switch&&
   !Band_disp_switch&&!MO_fileout&&!Dos_fileout&&!DosGauss_fileout&&
   !Wannier_Func_Calc&&!Wannier90_fileout&&!CWF_fileout_flag&&!CWF_Calc&&
   !pop_anal_aow_flag&&!unfold_electronic_band&&!Polarization_ABC_flag&&
   !xanes_gs_fileout&&!xanes_calc&&!cal_partial_charge;
 MPI_Allreduce(&valid,&all_valid,1,MPI_INT,MPI_MIN,mpi_comm_level1);
 if(!all_valid)rs_die("SOC_TR requires fixed one-step nonmagnetic c2n Cluster SOC; no contraction, magnetic fields, constraints, Hubbard or auxiliary states");
 valid=CpuSOCNative_Enabled();
 MPI_Allreduce(&valid,&all_valid,1,MPI_INT,MPI_MIN,mpi_comm_level1);
 if(!all_valid)rs_die("SOC_TR requires active native C3/D3 projected H/S validation");
}
/* Each pair is stored in rs_orig/rs_sum. A common absolute floor avoids
 * rejecting machine-level noise when an entire spin channel is near zero. */
static void rs_tr_metric(int which,double sign,double sum[14],double maxima[14],int*bad){
 for(size_t i=0;i<rs_v;i++){
  double a=rs_orig[i],b=rs_sum[i],res,scale2,err2,scale;
  if(!isfinite(a)||!isfinite(b)){*bad=1;continue;}
  res=a+sign*b;err2=res*res;scale2=.5*(a*a+b*b);scale=fmax(fabs(a),fabs(b));
  if(!isfinite(res)||!isfinite(err2)||!isfinite(scale2)){*bad=1;continue;}
  sum[2*which]+=err2;sum[2*which+1]+=scale2;
  maxima[2*which]=fmax(maxima[2*which],fabs(res));
  maxima[2*which+1]=fmax(maxima[2*which+1],scale);
 }
}
static int rs_tr_measure_raw(int step,double*****h,double*****im,double****olp,int emit){
 const double absolute_tolerance=1.e-10,relative_tolerance=1.e-9;
 const char*names[7]={"Huu_minus_conjHdd_real","Huu_minus_conjHdd_imag",
                     "Hud_plus_reverseT_real","Hud_plus_reverseT_imag","S_minus_reverseT",
                     "Huu_minus_reverseT_real","Huu_plus_reverseT_imag"};
 double local_sum[14]={0},global_sum[14],local_max[14]={0},global_max[14];
 unsigned long long local_count=(unsigned long long)rs_v,count=0;
 int bad=0,global_bad=0,pass=1;
 rs_pack(h[0],rs_orig,0);rs_pack(h[1],rs_sum,0);rs_tr_metric(0,-1.,local_sum,local_max,&bad);
 rs_pack(im[0],rs_orig,0);rs_pack(im[1],rs_sum,0);rs_tr_metric(1,1.,local_sum,local_max,&bad);
 rs_pack(h[2],rs_orig,0);memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_orig,RS_REVERSE,rs_sum,1.);rs_tr_metric(2,1.,local_sum,local_max,&bad);
 rs_pack(h[3],rs_orig,0);rs_pack(im[2],rs_sum,0);for(size_t i=0;i<rs_v;i++)rs_orig[i]+=rs_sum[i];
 memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_orig,RS_REVERSE,rs_sum,1.);rs_tr_metric(3,1.,local_sum,local_max,&bad);
 rs_pack(olp,rs_orig,0);memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_orig,RS_REVERSE,rs_sum,1.);rs_tr_metric(4,-1.,local_sum,local_max,&bad);
 /* TR alone does not imply Hermiticity of a same-spin diagonal block.
  * Check Huu explicitly; TR consistency transfers that condition to Hdd.
  * The off-diagonal spin partner is the native conjugate reverse by design. */
 rs_pack(h[0],rs_orig,0);memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_orig,RS_REVERSE,rs_sum,1.);rs_tr_metric(5,-1.,local_sum,local_max,&bad);
 rs_pack(im[0],rs_orig,0);memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_orig,RS_REVERSE,rs_sum,1.);rs_tr_metric(6,1.,local_sum,local_max,&bad);
 MPI_Allreduce(&bad,&global_bad,1,MPI_INT,MPI_MAX,mpi_comm_level1);
 MPI_Allreduce(&local_count,&count,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,mpi_comm_level1);
 MPI_Allreduce(local_sum,global_sum,14,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 MPI_Allreduce(local_max,global_max,14,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);
 if(global_bad||!count)pass=0;
 for(int k=0;k<7;k++){
  double rms=count?sqrt(global_sum[2*k]/(double)count):0.;
  double scale=count?sqrt(global_sum[2*k+1]/(double)count):0.;
  double relative=scale>0?rms/scale:(rms==0?0.:INFINITY);
  int finite=isfinite(rms)&&isfinite(scale)&&isfinite(global_max[2*k])&&isfinite(global_max[2*k+1]);
  int accepted=!global_bad&&count&&finite&&
    rms<=absolute_tolerance+relative_tolerance*scale&&
    global_max[2*k]<=absolute_tolerance+relative_tolerance*global_max[2*k+1];
  if(!accepted)pass=0;
  if(emit&&rs_rank==Host_ID)printf("CPU_SOC_TR_RAW step=%d check=%s scalar_elements=%llu rms=%.17g relative_rms=%.17g max_abs=%.17g rms_scale=%.17g max_scale=%.17g atol=%.3e rtol=%.3e nonfinite=%d accepted=%d\n",step,names[k],count,rms,relative,global_max[2*k],scale,global_max[2*k+1],absolute_tolerance,relative_tolerance,global_bad,accepted);
 }
 if(emit&&rs_rank==Host_ID){printf("CPU_SOC_TR_RAW_RESULT step=%d passed=%d storage=owner_neighbor_blocks reverse_route=dedicated overlap_imaginary=absent_real_storage extra_payload_arrays=0\n",step,pass);fflush(stdout);}
 return pass;
}
static void rs_tr_guard_raw(int step,double*****h,double*****im,double****olp){
 rs_tr_passed=0;
 if(!rs_tr_measure_raw(step,h,im,olp,1))rs_die("actual raw H/S time-reversal guard failed");
 rs_tr_passed=1;
}

void CpuSOCNative_RestartState(int ok){rs_restart_ok=ok;}
void CpuSOCNative_MarkSolve(int step,double *****h,double *****im,double ****olp){
 CpuSOCNative_TRValidateConfig();
 if(!CpuSOCNative_Enabled())return;
 if(!Use_of_Collinear_Restart||Scf_RestartFromFile!=1||SpinP_switch_RestartFiles!=0||rs_restart_ok!=1)rs_die("requires successfully read nonmagnetic NSOC charge via c2n");
 if(step!=1||rs_last_step!=-1||h!=H||im!=iHNL||olp!=OLP[0])rs_die("requires unique step1 solve on global uncontracted H/iHNL/OLP0");
 rs_setup();rs_geometry_guard();rs_hptr=h;rs_iptr=im;rs_overlap=olp;rs_solve_hash=rs_field_hash(h,im,olp);rs_last_step=step;
 rs_tr_passed=0;if(CpuSOCNative_TRRequested()==1)rs_tr_guard_raw(step,h,im,olp);
 if(rs_rank==Host_ID){printf("CPU_SOC_RS_SOLVE step=1 operator=PROJECTED_%s_GAMMA operations=%d raw_fields_preserved=1\n",rs_generic?"GSOC":(rs_ng==6?"D3":"C3"),rs_ng);fflush(stdout);}
}
void CpuSOCNative_MarkSolved(int step){
 if(!CpuSOCNative_Enabled())return;
 if(step!=1||rs_last_step!=1||rs_solved)rs_die("invalid solve completion order");rs_solved=1;
}
static void rs_spin_project(void){
 int f,g,out;size_t i;double start=MPI_Wtime();
 for(f=0;f<8;f++){rs_soc_raw[f]=rs_alloc(rs_v,sizeof(double));rs_soc_projected[f]=rs_alloc(rs_v,sizeof(double));}
 rs_pack(H[0],rs_soc_raw[0],0);rs_pack(iHNL[0],rs_soc_raw[1],0);rs_pack(H[1],rs_soc_raw[2],0);rs_pack(iHNL[1],rs_soc_raw[3],0);rs_pack(H[2],rs_soc_raw[4],0);
 rs_pack(H[3],rs_soc_raw[5],0);rs_pack(iHNL[2],rs_orig,0);for(i=0;i<rs_v;i++)rs_soc_raw[5][i]+=rs_orig[i];
 /* A directed (a,b,R) block is general complex. Its DU comes from the
  * other owned block (b,a,-R), never from UD(a,b,R)^dagger. */
 rs_exchange(rs_soc_raw[4],RS_REVERSE,rs_soc_raw[6],1.);rs_exchange(rs_soc_raw[5],RS_REVERSE,rs_soc_raw[7],-1.);
 for(f=0;f<8;f++)for(i=0;i<rs_v;i++)rs_soc_projected[f][i]=rs_soc_raw[f][i]/(double)rs_ng;
 for(g=1;g<rs_ng;g++)for(f=0;f<8;f++){
  memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_soc_raw[f],g,rs_sum,1.);
  for(out=0;out<8;out++){double weight=rs_spin_t[g][8*out+f]/(double)rs_ng;if(fabs(weight)<1e-18)continue;for(i=0;i<rs_v;i++)rs_soc_projected[out][i]+=weight*rs_sum[i];}
 }
 /* P and full spinor reverse Hermitization commute. Handle each offdiagonal
  * pair jointly so no updated channel is averaged against stale data. */
 for(f=0;f<6;f++){
  int from;double sign;SocSpin_Reverse(f,&from,&sign);memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_soc_projected[from],RS_REVERSE,rs_sum,sign);
  for(i=0;i<rs_v;i++)rs_soc_projected[f][i]=.5*(rs_soc_projected[f][i]+rs_sum[i]);
  if(f>=4){memset(rs_sum,0,rs_v*sizeof(double));rs_exchange(rs_soc_projected[f],RS_REVERSE,rs_sum,sign);memcpy(rs_soc_projected[from],rs_sum,rs_v*sizeof(double));}
 }
 rs_pack(H[0],rs_soc_projected[0],1);rs_pack(iHNL[0],rs_soc_projected[1],1);rs_pack(H[1],rs_soc_projected[2],1);rs_pack(iHNL[1],rs_soc_projected[3],1);rs_pack(H[2],rs_soc_projected[4],1);rs_pack(H[3],rs_soc_projected[5],1);
 memset(rs_sum,0,rs_v*sizeof(double));rs_pack(iHNL[2],rs_sum,1);
 rs_phase(1,"spinor_H_export",MPI_Wtime()-start);
 for(f=0;f<8;f++){free(rs_soc_raw[f]);free(rs_soc_projected[f]);rs_soc_raw[f]=NULL;rs_soc_projected[f]=NULL;}
 rs_bytes-=16*(rs_v?rs_v:1)*sizeof(double);
}
void CpuSOCNative_BeforeExport(void){
 if(!CpuSOCNative_Enabled())return;
 if(rs_band_prepared){
  rs_geometry_guard();
  if(!rs_exported||rs_band_export_handoff||H!=rs_hptr||iHNL!=rs_iptr||OLP[0]!=rs_overlap||rs_field_hash(H,iHNL,OLP[0])!=rs_band_hash)rs_die("projected H/S changed after targeted bands or export was repeated");
  rs_band_export_handoff=1;
  if(rs_rank==Host_ID){printf("CPU_SOC_BAND_EXPORT_HANDOFF projected_operator_hash_verified=1 unique_native_export=1\n");fflush(stdout);}return;
 }
 if(rs_last_step!=1||!rs_solved||rs_exported||H!=rs_hptr||iHNL!=rs_iptr||OLP[0]!=rs_overlap)rs_die("export is not the unique completed projected solve");
 rs_geometry_guard();if(rs_field_hash(H,iHNL,OLP[0])!=rs_solve_hash)rs_die("raw H/iHNL/S changed between projected solve and native export");
 rs_project(OLP[0],1,"scalar_S_export",1);rs_spin_project();if(rs_generic)CpuGSOC_ValidateExport();rs_exported=1;
 if(rs_rank==Host_ID){printf("CPU_SOC_RS_EXPORT step=1 representation=PROJECTED_REALSPACE_SPINOR_H_S spin_representation=SU2 operations=%d solver=%s_projected_blocks raw_checksum_unchanged=1 postprocess_external=0 canonical_UD_imag=H3 imNL2=0 derivatives_untouched_until_after_force=1\n",rs_ng,rs_generic?"GSOC":(rs_ng==6?"D3":"C3"));fflush(stdout);}
}
/* Called only by the deferred native band stage, after OutData and therefore
 * after all force/total-energy consumers of raw tensors. This is the same
 * one-shot projection/check that SCF2File would perform at this point. */
void CpuSOCNative_PrepareBands(void){
 if(!CpuSOCNative_Enabled())return;
 if(!CPU_BAND_USES_SPARSE_HS()||rs_band_prepared||rs_exported)rs_die("invalid targeted-band projection order");
 CpuSOCNative_BeforeExport();rs_band_hash=rs_field_hash(H,iHNL,OLP[0]);rs_band_prepared=1;
 if(rs_rank==Host_ID){printf("CPU_SOC_BAND_OPERATOR representation=PROJECTED_REALSPACE_SPINOR_H_S after_force=1 matching_native_export=1\n");fflush(stdout);}
}
