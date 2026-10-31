#ifndef OPENMX_HNL_WORKSPACE_H
#define OPENMX_HNL_WORKSPACE_H
#include <stddef.h>
#include <stdint.h>
#include <limits.h>
typedef struct {dcomplex***m[6];dcomplex***planes;dcomplex**rows;dcomplex*data;size_t bytes;unsigned long long atoms,pairs,calls;double seconds;} CpuHNLWorkspace;
static int CpuHNLWorkspaceAllocate(CpuHNLWorkspace*w,int width){
 if(width<1||width>128)return -1;size_t n=(size_t)width,values=18*n*n,rowptrs=18*n;
 memset(w,0,sizeof(*w));w->data=malloc(values*sizeof(dcomplex));w->rows=malloc(rowptrs*sizeof(dcomplex*));w->planes=malloc(18*sizeof(dcomplex**));
 if(!w->data||!w->rows||!w->planes){free(w->data);free(w->rows);free(w->planes);memset(w,0,sizeof(*w));return -1;}
 for(size_t s=0;s<18;s++){w->planes[s]=w->rows+s*n;for(size_t i=0;i<n;i++)w->planes[s][i]=w->data+(s*n+i)*n;}
 for(int k=0;k<6;k++)w->m[k]=w->planes+3*k;w->bytes=values*sizeof(dcomplex)+rowptrs*sizeof(dcomplex*)+18*sizeof(dcomplex**);return 0;
}
static void CpuHNLWorkspaceFree(CpuHNLWorkspace*w){free(w->data);free(w->rows);free(w->planes);memset(w,0,sizeof(*w));}
#endif
