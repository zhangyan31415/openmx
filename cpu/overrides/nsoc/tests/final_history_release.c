#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "cpu_sparse_lifetime.h"
static double *****slot(int atoms,const int*map,const int*species,const int*norb,const int*fan){
 double *****p=malloc(sizeof(*p));p[0]=malloc((atoms+1)*sizeof(*p[0]));
 for(int a=0;a<=atoms;a++){int g=a?map[a]:0,n=a?norb[species[g]]:1,f=a?fan[g]:0;p[0][a]=malloc((f+1)*sizeof(*p[0][a]));for(int h=0;h<=f;h++){p[0][a][h]=malloc(n*sizeof(*p[0][a][h]));for(int i=0;i<n;i++){p[0][a][h][i]=malloc(7*sizeof(double));for(int j=0;j<7;j++)p[0][a][h][i][j]=1000*a+100*h+10*i+j;}}}return p;
}
int main(void){
 int map[]={0,2,1},species[]={0,1,0},norb[]={2,3},fan[]={0,2,1};double ******p=malloc(4*sizeof(*p));
 for(int s=0;s<4;s++)p[s]=slot(2,map,species,norb,fan);
 double *****keep=p[0];double value=p[0][0][2][2][2][6];
 size_t released=CpuSparseHistoryRelease(p,4,1,1,2,map,species,norb,fan);
 assert(released>0&&p[0]==keep&&p[0][0][2][2][2][6]==value);
 for(int s=1;s<4;s++)assert(p[s]==NULL);
 assert(CpuSparseHistoryRelease(p,4,1,1,2,map,species,norb,fan)==0);
 assert(CpuSparseHistoryRelease(p,4,0,1,2,map,species,norb,fan)>0);
 free(p);puts("final history release: preserve current density, mapped species, idempotent cleanup PASS");return 0;
}
