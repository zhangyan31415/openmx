/* Distributed finite-range real-space D3. Original local blocks plus one
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
#include "cpu_periodic_image_gauge.h"
#include "cpu_realspace_d3.h"

typedef struct { int a,b,r[3],h,na,nb,owner; } RSBlock;
typedef struct {
  int *sc,*sd,*rc,*rd,*order,*target,*recv_target,*recv_offset;
  int *bc,*bd,*brc,*brd;
} RSRoute;
/* A finite 3-D crystallographic point group has at most 48 operations.
 * The metadata header carries the actual order, so D3 remains compatible
 * while C4/D4/Oh and other validated maps use the same runtime projector. */
#define RS_MAX_OPS 48
static int rs_mode=-1,rs_ready=0,rs_s_done=0,rs_rank,rs_np,rs_nat,rs_n,rs_ns,rs_nb,rs_nlocal,rs_maxorb,rs_check,rs_ng=6;
static int *rs_type,*rs_norb,*rs_local_id,*rs_g2l,*rs_offset,*rs_hash,rs_hashsize;
static int *rs_perm[RS_MAX_OPS],*rs_img[RS_MAX_OPS],rs_w[RS_MAX_OPS][9];
static double rs_cell[9],*rs_xyz,*rs_d[RS_MAX_OPS],*rs_sraw,rs_group_error;
static size_t rs_v,rs_bytes=0;
static uint64_t rs_topology,rs_map_fnv=UINT64_C(1469598103934665603);
static RSBlock *rs_block;
static RSRoute rs_route[RS_MAX_OPS+1];
static double *rs_orig,*rs_sum,*rs_send,*rs_recv,*rs_tmp,*rs_mat;
static double ****rs_overlap;
static int rs_last_step=-1;
static int rs_reference_valid=0;static double rs_reference_ne,rs_reference_band;
void CpuRSym_SetOccupationReference(double ne,double band){rs_reference_ne=ne;rs_reference_band=band;rs_reference_valid=isfinite(ne)&&isfinite(band);}

static void rs_die(const char *msg) {
  int rank=0; MPI_Comm_rank(mpi_comm_level1,&rank);
  fprintf(stderr,"CPU_RSD3_REJECT rank=%d reason=%s\n",rank,msg);fflush(stderr);
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
int CpuRSym_Enabled(void) {
  const char *v;if(rs_mode<0){v=getenv("OPENMX_CPU_REALSPACE_D3");
    if(!v||!strcmp(v,"0")||!strcmp(v,"off"))rs_mode=0;
    else if(!strcmp(v,"export"))rs_mode=1;
    else if(!strcmp(v,"scf")||!strcmp(v,"1"))rs_mode=2;
    else rs_die("OPENMX_CPU_REALSPACE_D3 must be off,export,scf");
  }
  if(rs_mode){
    if(SpinP_switch!=0||SO_switch!=0)rs_die("only unpolarized scalar NSOC supported");
    if(Cnt_switch!=0)rs_die("contracted/optimized PAO basis requires separately validated contracted metadata");
    if(Hub_U_switch!=0)rs_die("Hubbard ghost/derivative paths not validated");
    if(rs_mode==2&&Solver!=2)rs_die("per-SCF mode requires Cluster solver");
  }return rs_mode;
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
  double lo,hi;unsigned long long local=(unsigned long long)rs_bytes,maxbytes;
  MPI_Reduce(&elapsed,&lo,1,MPI_DOUBLE,MPI_MIN,Host_ID,mpi_comm_level1);
  MPI_Reduce(&elapsed,&hi,1,MPI_DOUBLE,MPI_MAX,Host_ID,mpi_comm_level1);
  MPI_Reduce(&local,&maxbytes,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,Host_ID,mpi_comm_level1);
  if(rs_rank==Host_ID){printf("CPU_RSD3_PHASE step=%d stage=%s wall_min=%.9f wall_max=%.9f spread=%.9f extra_bytes_max=%llu\n",step,label,lo,hi,hi-lo,maxbytes);fflush(stdout);}
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
    if(g==rs_ng){a=b->b;bb=b->a;for(i=0;i<3;i++)r[i]=-b->r[i];}
    else{a=rs_perm[g][b->a-1]+1;bb=rs_perm[g][b->b-1]+1;
      for(i=0;i<3;i++){r[i]=rs_img[g][3*(b->b-1)+i]-rs_img[g][3*(b->a-1)+i];for(j=0;j<3;j++)r[i]+=rs_w[g][3*i+j]*b->r[j];}}
    id=rs_find(a,bb,r);
    if(id<0){char msg[300];snprintf(msg,sizeof(msg),"missing %s closure g=%d source=(%d,%d,%d,%d,%d) target=(%d,%d,%d,%d,%d)",g==rs_ng?"Hermitian reverse":"point-group neighbor",g,b->a,b->b,b->r[0],b->r[1],b->r[2],a,bb,r[0],r[1],r[2]);rs_die(msg);}
    if((g<rs_ng&&(rs_block[id].na!=b->na||rs_block[id].nb!=b->nb))||(g==rs_ng&&(rs_block[id].na!=b->nb||rs_block[id].nb!=b->na)))rs_die("orbital shape not preserved by map");
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
    for(s=0;s<rs_ns;s++){
      int n=rs_norb[s];double *x=rs_rotation(g,s),*y=rs_rotation(h,s),*z=rs_rotation(target,s);
      for(i=0;i<n;i++)for(j=0;j<n;j++){double v=0;for(l=0;l<n;l++)v+=x[i*n+l]*y[l*n+j];rs_group_error=fmax(rs_group_error,fabs(v-z[i*n+j]));}
    }
  }
  free(seen);if(rs_group_error>1e-11)rs_die("orbital group composition is not exact for the supplied point-group metadata");
}
static void rs_setup(void) {
  const char *name;FILE *f;int32_t h[8];int g,a,j,k,b,mc,q,id;size_t nd=0;double start=MPI_Wtime(),orth=0;
  if(rs_ready){rs_geometry_guard();return;}
  MPI_Comm_rank(mpi_comm_level1,&rs_rank);MPI_Comm_size(mpi_comm_level1,&rs_np);
  name=getenv("OPENMX_CPU_REALSPACE_MAP");if(!name||name[0]!='/')rs_die("absolute OPENMX_CPU_REALSPACE_MAP required");
  f=fopen(name,"rb");if(!f)rs_die("cannot open metadata");rs_read(f,h,4,8);
  if(h[0]!=0x52534433||h[1]!=1||h[2]!=atomnum||h[4]<1||h[4]>RS_MAX_OPS||h[5]<1||h[5]>SpeciesNum||h[6]!=0||h[7]!=0)rs_die("metadata header invalid");
  rs_nat=h[2];rs_n=h[3];rs_ng=h[4];rs_ns=h[5];rs_read(f,rs_cell,8,9);
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
  if(fgetc(f)!=EOF)rs_die("metadata has trailing bytes");fclose(f);
  if(orth>1e-11)rs_die("orbital metadata is not exact orthogonal; regenerate the point-group map");
  {
    int moved=0;char error[256];
    for(a=0;a<3;a++)for(j=0;j<3;j++)
      if(!isfinite(tv[a+1][j+1]*BohrR)||fabs(tv[a+1][j+1]*BohrR-rs_cell[3*a+j])>2e-7)
        rs_die("cell differs from metadata; regenerate map");
    if(cpu_periodic_image_rebase(rs_nat,rs_ng,rs_cell,rs_xyz,Gxyz,BohrR,
        rs_perm,rs_w,rs_img,2e-7,&moved,error,sizeof(error)))rs_die(error);
    if(moved&&rs_rank==Host_ID){printf("CPU_SYMMETRY_PERIODIC_GAUGE rebased_atoms=%d geometry_unchanged=1\n",moved);fflush(stdout);}
  }
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
  for(g=0;g<=rs_ng;g++)rs_make_route(g);
  rs_orig=rs_alloc(rs_v,sizeof(double));rs_sum=rs_alloc(rs_v,sizeof(double));rs_send=rs_alloc(rs_v,sizeof(double));rs_recv=rs_alloc(rs_v,sizeof(double));
  rs_tmp=rs_alloc((size_t)rs_maxorb*rs_maxorb,sizeof(double));rs_mat=rs_alloc((size_t)rs_maxorb*rs_maxorb,sizeof(double));
  rs_check=getenv("OPENMX_CPU_REALSPACE_CHECK")&&atoi(getenv("OPENMX_CPU_REALSPACE_CHECK"));if(rs_check)rs_sraw=rs_alloc(rs_v,sizeof(double));
  rs_topology=rs_current_topology();rs_ready=1;rs_geometry_guard();
  if(rs_rank==Host_ID){printf("CPU_RSD3_SETUP mode=%s group_ops=%d natoms=%d N=%d global_blocks=%d local_blocks_host=%d map=%s map_fnv64=%016llx orbital_orth_error=%.3e orbital_group_error=%.3e storage=4local_fields_plus_maps ghost_policy=owned_primary_H_S HubU=off contraction=off\n",rs_mode==1?"export":"scf",rs_ng,rs_nat,rs_n,rs_nb,rs_nlocal,name,(unsigned long long)rs_map_fnv,orth,rs_group_error);fflush(stdout);}
  rs_phase(0,"setup",MPI_Wtime()-start);
}
static void rs_pack(double ****x,double *packed,int write) {
  int q,i,j;for(q=0;q<rs_nlocal;q++){RSBlock *b=&rs_block[rs_local_id[q]];int mc=F_G2M[b->a],base=rs_offset[q];for(i=0;i<b->na;i++)for(j=0;j<b->nb;j++){if(write)x[mc][b->h][i][j]=packed[base+i*b->nb+j];else packed[base+i*b->nb+j]=x[mc][b->h][i][j];}}
}
static void rs_dump(int step,const char *channel,const char *stage,const double *values) {
  const char *prefix=getenv("OPENMX_CPU_REALSPACE_DUMP_PREFIX"),*wanted=getenv("OPENMX_CPU_REALSPACE_DUMP_STEP");int select=wanted?atoi(wanted):1,q;char name[4096];FILE *f;int32_t header[8];
  if(!prefix)return;if(prefix[0]!='/')rs_die("dump prefix must be absolute");if(step!=0&&select!=-1&&step!=select)return;
  snprintf(name,sizeof(name),"%s.step%04d.rank%04d.%s.%s.rsdump",prefix,step,rs_rank,channel,stage);f=fopen(name,"wb");if(!f)rs_die("cannot create diagnostic dump");
  header[0]=0x52534450;header[1]=1;header[2]=rs_nat;header[3]=rs_n;header[4]=rs_rank;header[5]=rs_nlocal;header[6]=step;header[7]=channel[0];fwrite(header,4,8,f);
  for(q=0;q<rs_nlocal;q++){RSBlock *b=&rs_block[rs_local_id[q]];int32_t key[7]={b->a-1,b->b-1,b->r[0],b->r[1],b->r[2],b->na,b->nb};fwrite(key,4,7,f);fwrite(values+rs_offset[q],8,(size_t)b->na*b->nb,f);}
  fclose(f);
}
static void rs_exchange(const double *source,int g,double *destination,double scale) {
  RSRoute *t=&rs_route[g];int order,q,id,i,j,pos=0;
  for(order=0;order<rs_nlocal;order++){
    q=t->order[order];RSBlock *b=&rs_block[rs_local_id[q]];const double *x=source+rs_offset[q];
    if(g==rs_ng){for(i=0;i<b->nb;i++)for(j=0;j<b->na;j++)rs_send[pos++]=x[j*b->nb+i];}
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
  int g;size_t i;double start=MPI_Wtime();rs_pack(x,rs_orig,0);rs_dump(step,channel,"before",rs_orig);
  for(i=0;i<rs_v;i++)rs_sum[i]=rs_orig[i]/(double)rs_ng; /* metadata identity is checked by generator */
  for(g=1;g<rs_ng;g++)rs_exchange(rs_orig,g,rs_sum,1./(double)rs_ng);
  for(i=0;i<rs_v;i++)rs_sum[i]*=.5;
  /* Reverse Hermitization commutes with the exact point-group average. */
  rs_exchange(rs_sum,rs_ng,rs_sum,1.); /* reads the half-average, accumulates its transpose */
  if(rs_check&&(channel[0]=='H'||channel[0]=='S')){
    double local[2]={0,0},global[2];
    for(i=0;i<rs_v;i++){double d=rs_orig[i]-rs_sum[i];local[0]+=rs_orig[i]*rs_orig[i];local[1]+=d*d;}
    MPI_Reduce(local,global,2,MPI_DOUBLE,MPI_SUM,Host_ID,mpi_comm_level1);
    if(rs_rank==Host_ID){
      double relative=sqrt(global[1]/fmax(global[0],1e-300));
      const char *strict=getenv("OPENMX_CPU_REALSPACE_COV_STRICT"),*tol_text=getenv("OPENMX_CPU_REALSPACE_COV_TOL");
      double tol=tol_text?strtod(tol_text,NULL):1e-8;
      printf("CPU_RSD3_COVARIANCE step=%d field=%s relative=%.12e strict=%d tolerance=%.12e\n",step,channel,relative,(strict&&atoi(strict))!=0,tol);fflush(stdout);
      if(!isfinite(relative)||(strict&&atoi(strict)&&(!isfinite(tol)||tol<=0||relative>tol)))rs_die("H/S covariance residual exceeds OPENMX_CPU_REALSPACE_COV_TOL");
    }
  }
  rs_dump(step,channel,"after",rs_sum);if(write)rs_pack(x,rs_sum,1);
  rs_phase(step,channel,MPI_Wtime()-start);
}
#include "native_panel_metadata.h"
static void cpu_panel_record_projected(int step,double ****h,double ****s);
void CpuRSym_BeforeSolve(int step,double *****hamiltonian,double ****overlap) {
  if(CpuRSym_Enabled()!=2)return;rs_setup();rs_geometry_guard();rs_overlap=overlap;
  if(!rs_s_done){if(rs_check)rs_pack(overlap,rs_sraw,0);rs_project(overlap,step,"S",1);rs_s_done=1;}
  rs_project(hamiltonian[0],step,"H",1);rs_last_step=step;cpu_panel_record_projected(step,hamiltonian[0],overlap);
}
void CpuRSym_AfterDensity(int step,double *****density,double *****energy_density) {
  double sums[4]={0,0,0,0},total[4];int q,i,j,which;
  if(CpuRSym_Enabled()!=2||!rs_check||!rs_ready)return;
  for(q=0;q<rs_nlocal;q++){RSBlock *b=&rs_block[rs_local_id[q]];int mc=F_G2M[b->a];for(i=0;i<b->na;i++)for(j=0;j<b->nb;j++){int off=rs_offset[q]+i*b->nb+j;double d=2*density[0][mc][b->h][i][j],e=2*energy_density[0][mc][b->h][i][j],s=rs_overlap[mc][b->h][i][j];sums[0]+=d*rs_sraw[off];sums[1]+=d*s;sums[2]+=e*rs_sraw[off];sums[3]+=e*s;}}
  MPI_Reduce(sums,total,4,MPI_DOUBLE,MPI_SUM,Host_ID,mpi_comm_level1);
  if(rs_rank==Host_ID){printf("CPU_RSD3_TRACE step=%d TrD_Sraw=%.15g TrD_PS=%.15g TrEDM_Sraw=%.15g TrEDM_PS=%.15g spin_factor=2\n",step,total[0],total[1],total[2],total[3]);fflush(stdout);}
  if(rs_rank==Host_ID){const char *flag=getenv("OPENMX_CPU_REALSPACE_TRACE_STRICT");if(flag&&atoi(flag)){
   double core=-system_charge;for(int a=1;a<=atomnum;a++)core+=Spe_Core_Charge[WhatSpecies[a]];
   double ntol=1e-6+1e-10*fabs(rs_reference_ne),etol=1e-7+1e-10*fabs(rs_reference_band);
   if(!rs_reference_valid)rs_die("missing finite occupation reference");
   for(int k=0;k<4;k++)if(!isfinite(total[k]))rs_die("nonfinite density trace");
   printf("CPU_RSD3_TRACE_VERIFY step=%d Ne_occ=%.15g Ne_core=%.15g band_occ=%.15g Ne_raw_error=%.12e Ne_projected_error=%.12e EDM_raw_error=%.12e EDM_projected_error=%.12e Ne_atol=1e-6 EDM_atol=1e-7 rtol=1e-10\n",step,rs_reference_ne,core,rs_reference_band,fabs(total[0]-rs_reference_ne),fabs(total[1]-rs_reference_ne),fabs(total[2]-rs_reference_band),fabs(total[3]-rs_reference_band));fflush(stdout);
   if(fabs(rs_reference_ne-core)>ntol||fabs(total[0]-rs_reference_ne)>ntol||fabs(total[1]-rs_reference_ne)>ntol||fabs(total[2]-rs_reference_band)>etol||fabs(total[3]-rs_reference_band)>etol)rs_die("original/projected S density or signed EDM trace failed occupation reference");
  }}
  for(which=0;which<2;which++){
    double norm=0,diff=0,glob[2],loc[2];size_t k;rs_project(which?energy_density[0]:density[0],step,which?"EDMcheck":"DMcheck",0);
    for(k=0;k<rs_v;k++){norm+=rs_orig[k]*rs_orig[k];diff+=(rs_orig[k]-rs_sum[k])*(rs_orig[k]-rs_sum[k]);}loc[0]=norm;loc[1]=diff;MPI_Reduce(loc,glob,2,MPI_DOUBLE,MPI_SUM,Host_ID,mpi_comm_level1);
    if(rs_rank==Host_ID){printf("CPU_RSD3_INVARIANCE step=%d field=%s relative=%.12e\n",step,which?"EDM":"DM",sqrt(glob[1]/fmax(glob[0],1e-300)));fflush(stdout);}
  }
}
void CpuRSym_BeforeExport(void) {
  int mode=CpuRSym_Enabled();if(!mode)return;rs_setup();rs_geometry_guard();
  if(mode==1){double *****h=Cnt_switch?CntH:H;double ****s=Cnt_switch?CntOLP[0]:OLP[0];rs_project(s,0,"S",1);rs_project(h[0],0,"H",1);}
  if(rs_rank==Host_ID){printf("CPU_RSD3_EXPORT mode=%s last_scf_step=%d representation=PROJECTED_REALSPACE_H_S density_matrices_not_export_projected=1 derivatives_unchanged=1\n",mode==1?"export":"scf",rs_last_step);fflush(stdout);}
}

/* Read-only native history admission probe. No stored field is overwritten. */
void CpuRSym_HistoryProbe(int step,double ******history,double ******residual)
{
 const char *v=getenv("OPENMX_CPU_NSOC_HISTORY_PROBE");
 if(!v||!atoi(v)||CpuRSym_Enabled()!=2||!rs_ready)return;
 if(!(step<=3||step==7||step==8||step==14||step==15))return;
 rs_geometry_guard();
 for(int which=0;which<2;which++)for(int slot=0;slot<=2;slot++){
  double ****field=which?residual[slot][0]:history[slot][0];
  double local[3]={0,0,0},global[3];
  rs_project(field,step,which?"HistoryRcheck":"HistoryDcheck",0);
  for(size_t i=0;i<rs_v;i++){double d=rs_orig[i]-rs_sum[i];local[0]+=rs_orig[i]*rs_orig[i];local[1]+=d*d;local[2]=fmax(local[2],fabs(d));}
  MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
  MPI_Allreduce(local+2,global+2,1,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);
  if(rs_rank==Host_ID){printf("CPU_NSOC_HISTORY_INVARIANCE step=%d field=%s slot=%d norm=%.12e defect=%.12e relative=%.12e max_abs=%.12e\n",step,which?"ResidualDM":"DM",slot,sqrt(global[0]),sqrt(global[1]),sqrt(global[1]/fmax(global[0],1e-300)),global[2]);fflush(stdout);}
 }
}

/* Group-invariant directed-block history representation. Periodic keys and
 * nontrivial AO rotations are the existing validated real-space projector's. */
static int rh_ready=0,rh_nr=0,*rh_id,*rh_size,*rh_offset;
static size_t rh_n=0;
static double *rh_weights;
static void rh_setup(void)
{
 if(rh_ready)return;
 if(CpuRSym_Enabled()!=2)rs_die("symmetry history requires per-SCF realspace group projection");
 rs_setup();rs_geometry_guard();
 int *rep=calloc(rs_nb,sizeof(int)),*mult=calloc(rs_nb,sizeof(int));
 if(!rep||!mult)rs_die("history orbit metadata allocation");
 for(int id=0;id<rs_nb;id++){
  RSBlock *b=rs_block+id;int unique[RS_MAX_OPS],count=0,minimum=id;
  for(int g=0;g<rs_ng;g++){
   int a=rs_perm[g][b->a-1]+1,bb=rs_perm[g][b->b-1]+1,r[3];
   for(int i=0;i<3;i++){r[i]=rs_img[g][3*(b->b-1)+i]-rs_img[g][3*(b->a-1)+i];for(int j=0;j<3;j++)r[i]+=rs_w[g][3*i+j]*b->r[j];}
   int target=rs_find(a,bb,r);if(target<0)rs_die("history orbit missing finite-range block");
   if(target<minimum)minimum=target;
   int found=0;for(int k=0;k<count;k++)if(unique[k]==target)found=1;
   if(!found)unique[count++]=target;
  }
  rep[id]=minimum;mult[id]=count;
 }
 for(int q=0;q<rs_nlocal;q++)if(rep[rs_local_id[q]]==rs_local_id[q])rh_nr++;
 rh_id=rs_alloc(rh_nr,sizeof(int));rh_size=rs_alloc(rh_nr,sizeof(int));rh_offset=rs_alloc(rh_nr+1,sizeof(int));
 int nr=0;
 for(int q=0;q<rs_nlocal;q++){
  int id=rs_local_id[q];if(rep[id]!=id)continue;
  rh_id[nr]=id;rh_size[nr]=mult[id];rh_offset[nr]=(int)rh_n;
  rh_n+=(size_t)rs_block[id].na*rs_block[id].nb;if(rh_n>INT_MAX)rs_die("history representative local size overflow");nr++;
 }
 rh_offset[rh_nr]=(int)rh_n;rh_weights=rs_alloc(rh_n,sizeof(double));
 for(int q=0;q<rh_nr;q++)for(int i=rh_offset[q];i<rh_offset[q+1];i++)rh_weights[i]=(double)rh_size[q];
 unsigned long long local[3]={rs_v,rh_n,(unsigned long long)rh_nr},total[3];
 MPI_Reduce(local,total,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,Host_ID,mpi_comm_level1);
 if(rs_rank==Host_ID){printf("CPU_NSOC_SYMHISTORY_SETUP group_ops=%d directed_elements=%llu representative_elements=%llu representative_blocks=%llu scalar_storage_ratio=%.9f metric=orbit_cardinality_times_AO_Frobenius history=30 full_native_slots=3\n",rs_ng,total[0],total[1],total[2],(double)total[0]/fmax((double)total[1],1.));fflush(stdout);}
 free(rep);free(mult);rh_ready=1;
}
size_t CpuRSym_HistorySize(void){rh_setup();return rh_n;}
const double *CpuRSym_HistoryWeights(void){rh_setup();return rh_weights;}
void CpuRSym_HistoryPack(double ****field,double *representative)
{
 rh_setup();rs_geometry_guard();
 for(int q=0;q<rh_nr;q++){
  RSBlock *b=rs_block+rh_id[q];int mc=F_G2M[b->a];
  for(int i=0;i<b->na;i++)for(int j=0;j<b->nb;j++){
   double v=field[mc][b->h][i][j];if(!isfinite(v))rs_die("nonfinite history representative");
   representative[rh_offset[q]+i*b->nb+j]=v;
  }
 }
}
static void rh_expand_scratch(const double *representative)
{
 rh_setup();rs_geometry_guard();memset(rs_orig,0,rs_v*sizeof(double));
 for(int q=0;q<rh_nr;q++){
  RSBlock *b=rs_block+rh_id[q];int local=rs_g2l[rh_id[q]],base=rs_offset[local];
  for(int i=0;i<b->na*b->nb;i++)rs_orig[base+i]=rh_size[q]*representative[rh_offset[q]+i];
 }
 for(size_t i=0;i<rs_v;i++)rs_sum[i]=rs_orig[i]/(double)rs_ng;
 for(int g=1;g<rs_ng;g++)rs_exchange(rs_orig,g,rs_sum,1./(double)rs_ng);
 for(size_t i=0;i<rs_v;i++)rs_sum[i]*=.5;
 rs_exchange(rs_sum,rs_ng,rs_sum,1.);
}
void CpuRSym_HistoryExpand(const double *representative,double ****field)
{
 rh_expand_scratch(representative);rs_pack(field,rs_sum,1);
}
void CpuRSym_HistoryPackChecked(double ****field,double *representative,int step,const char *label)
{
 CpuRSym_HistoryPack(field,representative);rh_expand_scratch(representative);
 double local[3]={0,0,0},global[3];
 for(int q=0;q<rs_nlocal;q++){
  RSBlock *b=rs_block+rs_local_id[q];int mc=F_G2M[b->a],base=rs_offset[q];
  for(int i=0;i<b->na;i++)for(int j=0;j<b->nb;j++){
   double v=field[mc][b->h][i][j],d=v-rs_sum[base+i*b->nb+j];
   local[0]+=v*v;local[1]+=d*d;local[2]=fmax(local[2],fabs(d));
  }
 }
 MPI_Allreduce(local,global,2,MPI_DOUBLE,MPI_SUM,mpi_comm_level1);
 MPI_Allreduce(local+2,global+2,1,MPI_DOUBLE,MPI_MAX,mpi_comm_level1);
 double norm=sqrt(global[0]),defect=sqrt(global[1]),tolerance=1e-12+1e-11*norm;
 if(rs_rank==Host_ID){printf("CPU_NSOC_SYMHISTORY_ADMISSION step=%d field=%s norm=%.12e defect=%.12e relative=%.12e max_abs=%.12e tolerance=%.12e\n",step,label,norm,defect,defect/fmax(norm,1e-300),global[2],tolerance);fflush(stdout);}
 if(!isfinite(norm)||!isfinite(defect)||defect>tolerance)rs_die("raw history field has finite non-invariant complement; refusing representative-only storage");
}

/* Sparse representative-density routes. Unlike HistoryExpand, these never
 * inject a sparse field into a full all-block permutation: only representative
 * payloads are sent. Stabilizer images are averaged, then reverse Hermiticity
 * is enforced through the existing original-owner route. */
#include "cpu_nsoc_representatives.h"
static int rp_ready=0,rp_nr=0,rp_global_nr=0,*rp_id,*rp_mult,*rp_offsets;
static size_t rp_entries=0;
static RSRoute rp_routes[RS_MAX_OPS];
static int rp_target(int id,int g){
 RSBlock *b=rs_block+id;int a=rs_perm[g][b->a-1]+1,bb=rs_perm[g][b->b-1]+1,r[3];
 for(int i=0;i<3;i++){r[i]=rs_img[g][3*(b->b-1)+i]-rs_img[g][3*(b->a-1)+i];for(int j=0;j<3;j++)r[i]+=rs_w[g][3*i+j]*b->r[j];}
 int target=rs_find(a,bb,r);if(target<0)rs_die("representative operation leaves finite-range support");return target;
}
static void rp_setup(void){
 if(rp_ready){rs_geometry_guard();return;}
 if(CpuRSym_Enabled()!=2)rs_die("representative density requires per-SCF realspace projection");
 rs_setup();rs_geometry_guard();
 int *reps=rs_alloc(rs_nb,sizeof(int)),*mult=rs_alloc(rs_nb,sizeof(int));
 for(int id=0;id<rs_nb;id++){
  int unique[RS_MAX_OPS],nr=0,minimum=id;
  for(int g=0;g<rs_ng;g++){int t=rp_target(id,g),found=0;if(t<minimum)minimum=t;for(int k=0;k<nr;k++)if(unique[k]==t)found=1;if(!found)unique[nr++]=t;}
  reps[id]=minimum;mult[id]=nr;if(minimum==id){rp_global_nr++;if(rs_block[id].owner==rs_rank)rp_nr++;}
 }
 rp_id=rs_alloc(rp_nr,sizeof(int));rp_mult=rs_alloc(rp_nr,sizeof(int));rp_offsets=rs_alloc(rp_nr+1,sizeof(int));
 int q=0;for(int id=0;id<rs_nb;id++)if(reps[id]==id&&rs_block[id].owner==rs_rank){rp_id[q]=id;rp_mult[q]=mult[id];rp_offsets[q]=(int)rp_entries;rp_entries+=(size_t)rs_block[id].na*rs_block[id].nb;if(rp_entries>INT_MAX)rs_die("representative MPI count overflow");q++;}rp_offsets[rp_nr]=(int)rp_entries;
 for(int g=0;g<rs_ng;g++){
  RSRoute*t=rp_routes+g;
  t->sc=rs_alloc(rs_np,sizeof(int));t->sd=rs_alloc(rs_np,sizeof(int));t->rc=rs_alloc(rs_np,sizeof(int));t->rd=rs_alloc(rs_np,sizeof(int));
  t->bc=rs_alloc(rs_np,sizeof(int));t->bd=rs_alloc(rs_np,sizeof(int));t->brc=rs_alloc(rs_np,sizeof(int));t->brd=rs_alloc(rs_np,sizeof(int));
  t->order=rs_alloc(rp_nr,sizeof(int));t->target=rs_alloc(rp_nr,sizeof(int));
  for(int q=0;q<rp_nr;q++){int id=rp_target(rp_id[q],g),owner=rs_block[id].owner;t->target[q]=id;t->bc[owner]++;int z=rs_block[id].na*rs_block[id].nb;if(t->sc[owner]>INT_MAX-z)rs_die("representative send count overflow");t->sc[owner]+=z;}
  for(int i=1;i<rs_np;i++){t->sd[i]=t->sd[i-1]+t->sc[i-1];t->bd[i]=t->bd[i-1]+t->bc[i-1];}
  int *cursor=rs_alloc(rs_np,sizeof(int)),*sent=rs_alloc(rp_nr,sizeof(int));
  for(int q=0;q<rp_nr;q++){int owner=rs_block[t->target[q]].owner,pos=t->bd[owner]+cursor[owner]++;t->order[pos]=q;sent[pos]=t->target[q];}
  MPI_Alltoall(t->sc,1,MPI_INT,t->rc,1,MPI_INT,mpi_comm_level1);MPI_Alltoall(t->bc,1,MPI_INT,t->brc,1,MPI_INT,mpi_comm_level1);
  for(int i=1;i<rs_np;i++){t->rd[i]=t->rd[i-1]+t->rc[i-1];t->brd[i]=t->brd[i-1]+t->brc[i-1];}
  int blocks=t->brd[rs_np-1]+t->brc[rs_np-1],elements=t->rd[rs_np-1]+t->rc[rs_np-1];
  if(elements>(int)rs_v||blocks>rs_nlocal)rs_die("representative filtered receive exceeds owned support");
  t->recv_target=rs_alloc(blocks,sizeof(int));t->recv_offset=rs_alloc(blocks,sizeof(int));
  MPI_Alltoallv(sent,t->bc,t->bd,MPI_INT,t->recv_target,t->brc,t->brd,MPI_INT,mpi_comm_level1);
  int pos=0;for(int z=0;z<blocks;z++){int id=t->recv_target[z];if(id<0||id>=rs_nb||rs_g2l[id]<0)rs_die("representative received another owner's block");t->recv_offset[z]=pos;pos+=rs_block[id].na*rs_block[id].nb;}if(pos!=elements)rs_die("representative filtered route order mismatch");free(cursor);free(sent);
 }
 unsigned long long in[3]={rp_entries,(unsigned long long)rp_nr,rs_v},out[3];MPI_Allreduce(in,out,3,MPI_UNSIGNED_LONG_LONG,MPI_SUM,mpi_comm_level1);
 if(rs_rank==Host_ID){printf("CPU_NSOC_REP_SETUP group_ops=%d global_blocks=%d evaluated_blocks=%d global_entries=%llu evaluated_entries=%llu group_route_entries=%llu zero_nonrepresentative_entries_sent=0 dense_global_extra=0 stabilizers=averaged reverse=original_owner\n",rs_ng,rs_nb,rp_global_nr,out[2],out[0],out[0]*(unsigned long long)rs_ng);fflush(stdout);}free(reps);free(mult);rp_ready=1;
}
CpuNOBlock *CpuRSym_RepresentativeBlocks(int*count){
 rp_setup();int *offset=rs_alloc(rs_nat+1,sizeof(int));for(int a=1;a<=rs_nat;a++)offset[a]=offset[a-1]+rs_norb[rs_type[a-1]];
 CpuNOBlock*b=rs_alloc(rp_global_nr,sizeof(*b));int at=0;
 for(int id=0;id<rs_nb;id++){int minimum=id;for(int g=0;g<rs_ng;g++){int t=rp_target(id,g);if(t<minimum)minimum=t;}if(minimum!=id)continue;RSBlock*x=rs_block+id;b[at++]=(CpuNOBlock){offset[x->a-1],x->na,offset[x->b-1],x->nb,x->owner};}
 if(at!=rp_global_nr)rs_die("representative descriptor count mismatch");free(offset);*count=at;return b;
}
void CpuRSym_RepresentativeExpand(const double*values,size_t entries,double****field){
 rp_setup();if(entries!=rp_entries||(!values&&entries))rs_die("representative input shape mismatch");for(size_t i=0;i<entries;i++)if(!isfinite(values[i]))rs_die("nonfinite representative density input");memset(rs_sum,0,rs_v*sizeof(double));
 for(int g=0;g<rs_ng;g++){
  RSRoute*t=rp_routes+g;int pos=0;
  for(int z=0;z<rp_nr;z++){
   int q=t->order[z];RSBlock*b=rs_block+rp_id[q];const double*x=values+rp_offsets[q];double *left=rs_rotation(g,rs_type[b->a-1]),*right=rs_rotation(g,rs_type[b->b-1]);
   cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,b->na,b->nb,b->na,1.,left,b->na,x,b->nb,0.,rs_tmp,b->nb);
   cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,b->na,b->nb,b->nb,(double)rp_mult[q]/rs_ng,rs_tmp,b->nb,right,b->nb,0.,rs_mat,b->nb);
   memcpy(rs_send+pos,rs_mat,(size_t)b->na*b->nb*sizeof(double));pos+=b->na*b->nb;
  }
  MPI_Alltoallv(rs_send,t->sc,t->sd,MPI_DOUBLE,rs_recv,t->rc,t->rd,MPI_DOUBLE,mpi_comm_level1);
  int blocks=t->brd[rs_np-1]+t->brc[rs_np-1];for(int z=0;z<blocks;z++){int id=t->recv_target[z],q=rs_g2l[id],base=t->recv_offset[z];RSBlock*b=rs_block+id;for(int k=0;k<b->na*b->nb;k++)rs_sum[rs_offset[q]+k]+=rs_recv[base+k];}
 }
 for(size_t i=0;i<rs_v;i++)rs_sum[i]*=.5;rs_exchange(rs_sum,rs_ng,rs_sum,1.);rs_pack(field,rs_sum,1);
}
int CpuRSym_RepresentativeCertificate(int n,int maxn,const double *ko){
 rp_setup();return CpuNSym_RepresentativeCertificate(n,maxn,ko,rs_ng,rs_nat,(const int*const*)rs_perm,rs_type,rs_norb,(const double*const*)rs_d);
}

#include "native_panel_metadata.inc"

#include "cpu_nsoc_rep_tail.h"
int CpuRSym_RepresentativeTailCandidate(int n,int maxn,const double*ko,int*count,const int**ids,double*amplification){
 rp_setup();if(!CpuNSym_RepresentativeTailCandidate(n,maxn,ko,rs_ng,rs_nat,(const int*const*)rs_perm,rs_type,rs_norb,(const double*const*)rs_d,count,ids))return 0;
 static double cached_amplification=0.;
 if(!cached_amplification){double row;cached_amplification=CpuRepTailAmplification(rs_ng,rs_ns,rs_norb,(const double*const*)rs_d,&row);if(!isfinite(cached_amplification))rs_die("invalid stored-AO amplification bound");if(rs_rank==Host_ID){printf("CPU_NSOC_REP_TAIL_AMPLIFICATION max_AO_row_L1=%.15g bound=%.15g source=actual_stored_AO_action upward_inflated=1\n",row,cached_amplification);fflush(stdout);}}
 *amplification=cached_amplification;return 1;
}
