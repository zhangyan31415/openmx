/* Exact-shape thread scratch cache; no MPI and no shared writable arrays. */
#include "cpu_orbital_workspace.h"
#include <limits.h>
#include <float.h>
#include <math.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define CACHE_LIMIT_BYTES ((size_t)65536)
static _Thread_local int mode=-1;
static _Thread_local CpuOrbitalWorkspace *cache;
static pthread_key_t key;
static pthread_once_t once=PTHREAD_ONCE_INIT;
static void destroy(void *value){
 CpuOrbitalWorkspace *w=value;if(!w)return;free(w->values);free(w->rows);free(w);cache=NULL;
}
void CpuOrbitalWorkspace_ReleaseThread(void){
 if(cache){CpuOrbitalWorkspace *w=cache;cache=NULL;pthread_setspecific(key,NULL);destroy(w);}
}
static void create_key(void){
 if(pthread_key_create(&key,destroy)!=0){fprintf(stderr,"ORBITAL_WORKSPACE: pthread_key_create failed\n");abort();}
 if(atexit(CpuOrbitalWorkspace_ReleaseThread)!=0){fprintf(stderr,"ORBITAL_WORKSPACE: atexit failed\n");abort();}
}
static int enabled(void){
 if(mode<0){const char *p=getenv("OPENMX_CPU_ORBITAL_WORKSPACE");
  if(!p||!strcmp(p,"0"))mode=0;
  else if(!strcmp(p,"1"))mode=1;
  else {fprintf(stderr,"OPENMX_CPU_ORBITAL_WORKSPACE must be0or1\n");abort();}
 }
 return mode;
}
CpuOrbitalWorkspace *CpuOrbitalWorkspace_Acquire(int maxl,int radial,double x,double y,double z,double cutoff){
 size_t l,m,n,ptr_bytes,value_bytes,i;double radius,*p;CpuOrbitalWorkspace *w;
 if(!enabled())return NULL;
 /* Stock's outside-cutoff branch leaves angular scratch/derivatives undefined.
    Preserve that legacy allocation path; only fully initialized calls reuse. */
 radius=sqrt((x*x+y*y)+z*z);
 if(!isfinite(radius)||!isfinite(cutoff)||!(radius>1.e-12&&radius<cutoff*(1.-64.*DBL_EPSILON)))return NULL;
 if(maxl<0||maxl>(INT_MAX-3)/2||radial<=0){fprintf(stderr,"ORBITAL_WORKSPACE: invalid scratch dimensions\n");abort();}
 l=(size_t)maxl+1;m=2*l+1;
 if(l>SIZE_MAX/5/sizeof(double*)||l>SIZE_MAX/(size_t)radial/2)return NULL;
 n=2*l*(size_t)radial;
 if(l>SIZE_MAX/m/3||n>SIZE_MAX-3*l*m)return NULL;n+=3*l*m;
 if(n>SIZE_MAX/sizeof(double))return NULL;
 ptr_bytes=5*l*sizeof(double*);value_bytes=n*sizeof(double);
 if(ptr_bytes>CACHE_LIMIT_BYTES||value_bytes>CACHE_LIMIT_BYTES-ptr_bytes)return NULL;
 w=cache;if(w&&w->active)return NULL; /* Reentrant same-thread call stays legacy. */
 if(!w){
  if(pthread_once(&once,create_key)!=0)abort();
  w=calloc(1,sizeof(*w));if(!w)abort();w->allocations=1;cache=w;
  if(pthread_setspecific(key,w)!=0){destroy(w);abort();}
 }
 if(!w->rows||w->maxl!=maxl||w->radial!=radial){
  free(w->values);free(w->rows);w->rows=NULL;w->values=NULL;
  w->rows=malloc(ptr_bytes);w->values=malloc(value_bytes);w->allocations+=2;
  if(!w->rows||!w->values){destroy(w);pthread_setspecific(key,NULL);abort();}
  w->RF=w->rows;w->dRF=w->rows+l;w->AF=w->rows+2*l;w->dAFQ=w->rows+3*l;w->dAFP=w->rows+4*l;p=w->values;
  for(i=0;i<l;i++){w->RF[i]=p;p+=radial;}
  for(i=0;i<l;i++){w->dRF[i]=p;p+=radial;}
  for(i=0;i<l;i++){w->AF[i]=p;p+=m;}
  for(i=0;i<l;i++){w->dAFQ[i]=p;p+=m;}
  for(i=0;i<l;i++){w->dAFP[i]=p;p+=m;}
  w->maxl=maxl;w->radial=radial;w->bytes=ptr_bytes+value_bytes;
 }
 w->active=1;return w;
}
void CpuOrbitalWorkspace_Release(CpuOrbitalWorkspace *w){if(w)w->active=0;}
size_t CpuOrbitalWorkspace_ThreadAllocationCount(void){return cache?cache->allocations:0;}
