#include "rep_panel.h"
#include <assert.h>
#include <stdio.h>
int main(int argc,char**argv){
 MPI_Init(&argc,&argv);int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
 int64_t cp[5]={0,1,2,3,4};int ri[4]={0,1,2,3},types[4]={0,0,0,0},norb[1]={1},perm[4]={0,1,2,3};double cv[4]={1,1,1,1},rot[1]={1};const int*perms[6]={perm,perm,perm,perm,perm,perm};const double*rots[6]={rot,rot,rot,rot,rot,rot};
 int coords[256];assert(np<=128);for(int i=0;i<np;i++){coords[2*i]=0;coords[2*i+1]=i;}
 RpaBasis b={0};b.n=4;b.natoms=4;b.species=1;b.ng=6;for(int k=0;k<4;k++)b.sectors[k]=1;b.cp=cp;b.ri=ri;b.cv=cv;b.types=types;b.norb=norb;b.permutations=perms;b.rotations=rots;b.context=(RpaContext){1,2,3,4};
 char error[256];int kind=-1;RpaPlan*p=Rpa_CreateChecked(MPI_COMM_WORLD,&b,0,NULL,32,1,np,coords,rank==0?0:1024*1024,error,sizeof(error),&kind);
 assert(!p&&kind==RPA_CREATE_RESOURCE); /* one rank exhausts; all must agree */
 b.n=3;kind=-1;p=Rpa_CreateChecked(MPI_COMM_WORLD,&b,0,NULL,32,1,np,coords,1024*1024,error,sizeof(error),&kind);assert(!p&&kind==RPA_CREATE_INVALID);
 b.n=4;kind=-1;p=Rpa_CreateChecked(MPI_COMM_WORLD,&b,0,NULL,32,1,np,coords,1024*1024,error,sizeof(error),&kind);assert(!p&&kind==RPA_CREATE_INVALID); /* false D3 representation stays fatal */
 if(!rank)puts("collective resource vs invalid classification PASS");MPI_Finalize();return 0;
}
