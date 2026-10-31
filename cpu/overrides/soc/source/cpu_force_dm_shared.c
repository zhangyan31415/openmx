/* Node-local force density views. The six native CDM/iDM components are
 * copied bit-for-bit; only node leaders request off-node owned atom blocks. */
#include "cpu_force_dm_shared.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#ifndef CPU_FORCE_DM_TILE_ELEMENTS
#define CPU_FORCE_DM_TILE_ELEMENTS 131072
#endif
struct CpuForceDM {
 MPI_Comm comm,node; MPI_Win win;
 int visible; int *rows,*neighbors;
 double *****cdm,*****idm; double dummy[6];
 size_t bytes;
};
static void fail(MPI_Comm c,const char*s){int r;MPI_Comm_rank(c,&r);if(!r){fprintf(stderr,"[CPU_FORCE_DM_SHARED] %s\n",s);fflush(stderr);}MPI_Abort(c,71);abort();}
static void require(MPI_Comm c,int ok,const char*s){int all;MPI_Allreduce(&ok,&all,1,MPI_INT,MPI_MIN,c);if(!all)fail(c,s);}
static void *alloc(MPI_Comm c,size_t n,size_t sz){void*p;if(sz&&n>SIZE_MAX/sz)fail(c,"allocation size overflow");p=calloc(n?n:1,sz);if(!p)fail(c,"metadata allocation failed");return p;}
int CpuForceDM_AddBytes(size_t r,size_t c,size_t prev,size_t*next){size_t v;if(!next||(c&&r>SIZE_MAX/c))return 0;v=r*c;if(v>SIZE_MAX/sizeof(double)||prev>SIZE_MAX-v*sizeof(double))return 0;v=prev+v*sizeof(double);if(v>(size_t)PTRDIFF_MAX||(MPI_Aint)v<0||(size_t)(MPI_Aint)v!=v)return 0;*next=v;return 1;}
int CpuForceDM_Requested(MPI_Comm c){const char*s=getenv("OPENMX_CPU_FORCE_DM_SHARED");int v=!s||!strcmp(s,"0")?0:!strcmp(s,"1")?1:-1,lo,hi;MPI_Allreduce(&v,&lo,1,MPI_INT,MPI_MIN,c);MPI_Allreduce(&v,&hi,1,MPI_INT,MPI_MAX,c);if(lo<0||lo!=hi)fail(c,"flag must be collective-consistent 0 or 1");return v;}
static uint64_t hash_int(uint64_t h,int v){return (h^(uint32_t)v)*UINT64_C(1099511628211);}
static double *****view(CpuForceDM*x,int comps,int comp0,const int*fm,const int*n,const int*f,int*const*nn,const size_t*elements,double**atom){
 int c,m,h,i;double *****v=alloc(x->comm,comps,sizeof(*v));
 for(c=0;c<comps;c++){
  v[c]=alloc(x->comm,(size_t)x->visible+1,sizeof(**v));
  for(m=0;m<=x->visible;m++){
   int g=m?fm[m]:0,nr=m?n[g]:1,nh=m?f[g]+1:1;size_t off=m?(size_t)(c+comp0)*elements[g]:0;
   v[c][m]=alloc(x->comm,nh,sizeof(***v));
   for(h=0;h<nh;h++){
    int nc=m?n[nn[g][h]]:1;
    v[c][m][h]=alloc(x->comm,nr,sizeof(****v));
    for(i=0;i<nr;i++){v[c][m][h][i]=m?atom[g]+off:&x->dummy[c+comp0];off+=(size_t)nc;}
   }
  }
 }
 return v;
}
static void free_view(CpuForceDM*x,double*****v,int comps){int c,m,h;if(!v)return;for(c=0;c<comps;c++){for(m=0;m<=x->visible;m++){for(h=0;h<x->neighbors[m];h++)free(v[c][m][h]);free(v[c][m]);}free(v[c]);}free(v);}
/* Source rows are ragged native arrays. This cursor packs arbitrary bounded
 * contiguous slices without creating a second complete atom/density copy. */
typedef struct {int a,c,h,i,j;} Cursor;
static size_t pack(Cursor*t,const int*list,int count,const int*g2m,const int*n,const int*f,int*const*nn,double*****cdm,double*****idm,double*out,size_t capacity){
 size_t used=0;
 while(t->a<count&&used<capacity){
  int g=list[t->a],m=g2m[g],cols=n[nn[g][t->h]];size_t take=(size_t)(cols-t->j);
  double*row=t->c<4?cdm[t->c][m][t->h][t->i]:idm[t->c-4][m][t->h][t->i];
  if(take>capacity-used)take=capacity-used;
  memcpy(out+used,row+t->j,take*sizeof(double));used+=take;t->j+=(int)take;
  if(t->j==cols){t->j=0;if(++t->i==n[g]){t->i=0;if(++t->h>f[g]){t->h=0;if(++t->c==6){t->c=0;t->a++;}}}}
 }
 return used;
}
static size_t unpack(int*a,size_t*offset,const int*list,int count,const size_t*elements,double**atom,const double*in,size_t length){
 size_t used=0;
 while(*a<count&&used<length){int g=list[*a];size_t total=6*elements[g],take=total-*offset;if(take>length-used)take=length-used;memcpy(atom[g]+*offset,in+used,take*sizeof(double));used+=take;*offset+=take;if(*offset==total){*offset=0;(*a)++;}}
 return used;
}
CpuForceDM *CpuForceDM_Create(MPI_Comm comm,int atoms,int owned,int ghosts,const int*owner,const int*m2g,const int*fm2g,const int*norb,const int*fnan,int*const*natn,double*****cdm,double*****idm,double******cv,double******iv){
 CpuForceDM*x;int rank,np,nnode,nrank,g,h,m,c,i,ok=1,lo,hi,rc,k,needed_atoms=0;
 size_t *segment_bytes,*off,*elements,unique=0,legacy=0,nodebytes=0;double **base,**atom,*local;
 int *seen,*global,*need,*needed,*members,*world_to_local,*g2m;
 int *sc,*sd,*rcount,*rd,*request,*requested,*pos,total_send=0,total_recv=0;
 uint64_t hash=UINT64_C(1469598103934665603),hlo,hhi;
 unsigned long long legacy_all,legacy_local,node_all,node_local,max_rounds=0,rounds=0,round;
 const size_t cap=CPU_FORCE_DM_TILE_ELEMENTS; /* one MiB per buffer, independent of world size */
 double *sendbuf,*recvbuf;
 MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&np);
 require(comm,cap>0&&cap<=(size_t)INT_MAX,"invalid bounded force tile capacity");
 MPI_Allreduce(&atoms,&lo,1,MPI_INT,MPI_MIN,comm);MPI_Allreduce(&atoms,&hi,1,MPI_INT,MPI_MAX,comm);
 require(comm,atoms>0&&atoms==lo&&lo==hi&&atoms<INT_MAX&&owned>=0&&ghosts>=0&&owned<=atoms&&ghosts<=atoms-owned&&owner&&m2g&&fm2g&&norb&&fnan&&natn&&cdm&&idm&&cv&&iv,"invalid atom metadata");
 x=alloc(comm,1,sizeof(*x));x->comm=comm;x->node=MPI_COMM_NULL;x->win=MPI_WIN_NULL;x->visible=owned+ghosts;
 rc=MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&x->node);require(comm,rc==MPI_SUCCESS&&x->node!=MPI_COMM_NULL,"shared communicator creation failed");MPI_Comm_size(x->node,&nnode);MPI_Comm_rank(x->node,&nrank);
 segment_bytes=alloc(comm,nnode,sizeof(*segment_bytes));off=alloc(comm,(size_t)atoms+1,sizeof(*off));elements=alloc(comm,(size_t)atoms+1,sizeof(*elements));base=alloc(comm,nnode,sizeof(*base));atom=alloc(comm,(size_t)atoms+1,sizeof(*atom));seen=alloc(comm,(size_t)atoms+1,sizeof(*seen));global=alloc(comm,(size_t)atoms+1,sizeof(*global));need=alloc(comm,(size_t)atoms+1,sizeof(*need));needed=alloc(comm,(size_t)atoms+1,sizeof(*needed));members=alloc(comm,nnode,sizeof(*members));world_to_local=alloc(comm,np,sizeof(*world_to_local));g2m=alloc(comm,(size_t)atoms+1,sizeof(*g2m));
 MPI_Allgather(&rank,1,MPI_INT,members,1,MPI_INT,x->node);
 for(g=0;g<np;g++)world_to_local[g]=-1;for(g=0;g<nnode;g++)world_to_local[members[g]]=g;
 for(g=1;g<=atoms;g++){
  size_t next=0;int r=owner[g];
  if(r<0||r>=np||norb[g]<=0||fnan[g]<0||fnan[g]>=INT_MAX||!natn[g]){ok=0;continue;}
  hash=hash_int(hash,r);hash=hash_int(hash,norb[g]);hash=hash_int(hash,fnan[g]);
  for(h=0;h<=fnan[g];h++){int ng=natn[g][h];hash=hash_int(hash,ng);if(ng<1||ng>atoms||norb[ng]<=0){ok=0;continue;}if(!CpuForceDM_AddBytes(norb[g],norb[ng],next,&next))ok=0;}
  elements[g]=next/sizeof(double);if(!CpuForceDM_AddBytes(6,elements[g],unique,&unique))ok=0;
 }
 require(comm,ok,"invalid neighbor metadata or payload overflow");
 MPI_Allreduce(&hash,&hlo,1,MPI_UINT64_T,MPI_MIN,comm);MPI_Allreduce(&hash,&hhi,1,MPI_UINT64_T,MPI_MAX,comm);require(comm,hlo==hhi,"global atom metadata differs across ranks");
 for(m=1;m<=owned;m++){g=m2g[m];if(g<1||g>atoms||owner[g]!=rank||fm2g[m]!=g){ok=0;continue;}seen[g]++;g2m[g]=m;}
 MPI_Allreduce(seen,global,atoms+1,MPI_INT,MPI_SUM,comm);for(g=1;g<=atoms;g++)if(global[g]!=1)ok=0;
 memset(seen,0,((size_t)atoms+1)*sizeof(*seen));
 for(m=1;m<=x->visible;m++){g=fm2g[m];if(g<1||g>atoms){ok=0;continue;}if(seen[g]++)ok=0;need[g]=1;if(!CpuForceDM_AddBytes(6,elements[g],legacy,&legacy))ok=0;}
 require(comm,ok,"owned/ghost mapping must be bijective and owner-complete");
 MPI_Allreduce(need,needed,atoms+1,MPI_INT,MPI_MAX,x->node);
 /* Local owners first-touch their own rows. Off-node rows are held once by
  * the node leader, and viewed by every local force reader that needs them. */
 for(g=1;g<=atoms;g++)if(needed[g]){int lr=world_to_local[owner[g]];needed_atoms++;if(lr<0)lr=0;off[g]=segment_bytes[lr]/sizeof(double);if(!CpuForceDM_AddBytes(6,elements[g],segment_bytes[lr],&segment_bytes[lr]))ok=0;}
 require(comm,ok,"needed-node shared segment overflow");
 for(g=0;g<nnode;g++){if(nodebytes>SIZE_MAX-segment_bytes[g])fail(comm,"node payload overflow");nodebytes+=segment_bytes[g];}
 MPI_Comm_set_errhandler(x->node,MPI_ERRORS_RETURN);
 rc=MPI_Win_allocate_shared((MPI_Aint)segment_bytes[nrank],sizeof(double),MPI_INFO_NULL,x->node,&local,&x->win);require(comm,rc==MPI_SUCCESS&&x->win!=MPI_WIN_NULL,"shared window allocation failed");MPI_Win_set_errhandler(x->win,MPI_ERRORS_RETURN);
 for(g=0;g<nnode;g++){MPI_Aint sz;int disp;rc=MPI_Win_shared_query(x->win,g,&sz,&disp,&base[g]);if(rc!=MPI_SUCCESS||sz!=(MPI_Aint)segment_bytes[g]||disp!=(int)sizeof(double)||(sz&&!base[g]))ok=0;}
 require(comm,ok,"node-local owner segment query failed");
 rc=MPI_Win_lock_all(0,x->win);require(comm,rc==MPI_SUCCESS,"shared window epoch failed");
 for(g=1;g<=atoms;g++)if(needed[g]){int lr=world_to_local[owner[g]];if(lr<0)lr=0;atom[g]=base[lr]+off[g];}
 for(m=1;m<=owned;m++){size_t at=0;g=m2g[m];for(c=0;c<6;c++)for(h=0;h<=fnan[g];h++)for(i=0;i<norb[g];i++){int cols=norb[natn[g][h]];double*src=c<4?cdm[c][m][h][i]:idm[c-4][m][h][i];memcpy(atom[g]+at,src,(size_t)cols*sizeof(double));at+=(size_t)cols;}}
 /* Requests contain only off-node atoms, once per requesting physical node.
  * Metadata is bounded by atom count, payload MPI counts by cap<=INT_MAX. */
 sc=alloc(comm,np,sizeof(*sc));sd=alloc(comm,np,sizeof(*sd));rcount=alloc(comm,np,sizeof(*rcount));rd=alloc(comm,np,sizeof(*rd));pos=alloc(comm,np,sizeof(*pos));
 if(!nrank)for(g=1;g<=atoms;g++)if(needed[g]&&world_to_local[owner[g]]<0)sc[owner[g]]++;
 MPI_Alltoall(sc,1,MPI_INT,rcount,1,MPI_INT,comm);
 for(g=0;g<np;g++){sd[g]=total_send;rd[g]=total_recv;if(sc[g]>INT_MAX-total_send||rcount[g]>INT_MAX-total_recv){ok=0;break;}total_send+=sc[g];total_recv+=rcount[g];}
 require(comm,ok,"request metadata MPI count overflow");
 request=alloc(comm,total_send,sizeof(*request));requested=alloc(comm,total_recv,sizeof(*requested));memcpy(pos,sd,(size_t)np*sizeof(*pos));
 if(!nrank)for(g=1;g<=atoms;g++)if(needed[g]&&world_to_local[owner[g]]<0)request[pos[owner[g]]++]=g;
 MPI_Alltoallv(request,sc,sd,MPI_INT,requested,rcount,rd,MPI_INT,comm);
 for(g=0;g<total_recv;g++){int a=requested[g];if(a<1||a>atoms||owner[a]!=rank||!g2m[a])ok=0;}
 require(comm,ok,"received force density request has wrong owner");
 for(g=0;g<np;g++){size_t count=0;for(k=0;k<rcount[g];k++){int a=requested[rd[g]+k];if(count>SIZE_MAX-6*elements[a])fail(comm,"send payload overflow");count+=6*elements[a];}if((unsigned long long)(count/cap+(count%cap!=0))>rounds)rounds=(unsigned long long)(count/cap+(count%cap!=0));}
 MPI_Allreduce(&rounds,&max_rounds,1,MPI_UNSIGNED_LONG_LONG,MPI_MAX,comm);
 sendbuf=alloc(comm,cap,sizeof(*sendbuf));recvbuf=alloc(comm,cap,sizeof(*recvbuf));
 /* Directed ring ordering gives every source/destination pair exactly once;
  * fixed global round count prevents empty-owner rank schedule divergence. */
 for(k=1;k<np;k++){
  int dest=(rank+k)%np,src=(rank-k+np)%np,a=0;size_t at=0;Cursor cursor={0,0,0,0,0};
  for(round=0;round<max_rounds;round++){
   size_t sent=pack(&cursor,requested+rd[dest],rcount[dest],g2m,norb,fnan,natn,cdm,idm,sendbuf,cap);int received;
   MPI_Status status;
   rc=MPI_Sendrecv(sendbuf,(int)sent,MPI_DOUBLE,dest,712,recvbuf,(int)cap,MPI_DOUBLE,src,712,comm,&status);if(rc!=MPI_SUCCESS)fail(comm,"force density node exchange failed");
   MPI_Get_count(&status,MPI_DOUBLE,&received);
   if(received<0||unpack(&a,&at,request+sd[src],sc[src],elements,atom,recvbuf,(size_t)received)!=(size_t)received)fail(comm,"force density receive stream mismatch");
  }
  if(cursor.a!=rcount[dest]||cursor.c||cursor.h||cursor.i||cursor.j||a!=sc[src]||at)fail(comm,"force density stream ended before atom completion");
 }
 free(sendbuf);free(recvbuf);
 if(MPI_Win_sync(x->win)!=MPI_SUCCESS)fail(comm,"writer sync failed");MPI_Barrier(x->node);if(MPI_Win_sync(x->win)!=MPI_SUCCESS)fail(comm,"reader sync failed");
 x->rows=alloc(comm,(size_t)x->visible+1,sizeof(*x->rows));x->neighbors=alloc(comm,(size_t)x->visible+1,sizeof(*x->neighbors));for(m=0;m<=x->visible;m++){x->rows[m]=m?norb[fm2g[m]]:1;x->neighbors[m]=m?fnan[fm2g[m]]+1:1;}for(c=0;c<6;c++)x->dummy[c]=NAN;
 x->cdm=view(x,4,0,fm2g,norb,fnan,natn,elements,atom);x->idm=view(x,2,4,fm2g,norb,fnan,natn,elements,atom);x->bytes=unique;*cv=x->cdm;*iv=x->idm;
 legacy_local=(unsigned long long)legacy;MPI_Allreduce(&legacy_local,&legacy_all,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);
 if(!nrank){printf("[CPU_FORCE_DM_NODE] leader=%d local_ranks=%d needed_atoms=%d node_payload_bytes=%zu remote_requested_atoms=%d\n",rank,nnode,needed_atoms,nodebytes,total_send);fflush(stdout);}
 node_local=nrank?0:(unsigned long long)nodebytes;MPI_Allreduce(&node_local,&node_all,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,comm);
 if(!rank){printf("[CPU_FORCE_DM_SHARED] active=1 unique_payload_bytes=%zu legacy_payload_bytes_sum=%llu node_payload_bytes_sum=%llu shared_pages_counted_once=1 cross_node_needed_only=1 tile_bytes=%zu\n",unique,legacy_all,node_all,cap*sizeof(double));fflush(stdout);}
 free(segment_bytes);free(off);free(elements);free(base);free(atom);free(seen);free(global);free(need);free(needed);free(members);free(world_to_local);free(g2m);free(sc);free(sd);free(rcount);free(rd);free(request);free(requested);free(pos);return x;
}
size_t CpuForceDM_Bytes(const CpuForceDM*x){return x?x->bytes:0;}
void CpuForceDM_Destroy(CpuForceDM**p){CpuForceDM*x;if(!p||!*p)return;x=*p;
 if(MPI_Win_sync(x->win)!=MPI_SUCCESS)fail(x->comm,"release sync failed");MPI_Barrier(x->node);if(MPI_Win_sync(x->win)!=MPI_SUCCESS)fail(x->comm,"release reader sync failed");
 free_view(x,x->cdm,4);free_view(x,x->idm,2);free(x->rows);free(x->neighbors);
 if(MPI_Win_unlock_all(x->win)!=MPI_SUCCESS)fail(x->comm,"shared epoch close failed");if(MPI_Win_free(&x->win)!=MPI_SUCCESS)fail(x->comm,"shared window release failed");MPI_Comm_free(&x->node);free(x);*p=NULL;
}
