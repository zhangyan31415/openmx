#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
static int Matomnum=1,order_max=2,HS_fileout=1,Host_ID=0,mpi_comm_level1=0;
static int atoms_map[]={0,1},species[]={0,0},norbs[]={2},fans[]={0,1};
static int nbr0[]={0},nbr1[]={1,1};static int *neighbors[]={nbr0,nbr1};
#define M2G atoms_map
#define WhatSpecies species
#define Spe_Total_CNO norbs
#define FNAN fans
#define natn neighbors
static void MPI_Comm_rank(int comm,int*r){(void)comm;*r=0;}
static int grid_alive=1,position_calls=0,momentum_calls=0;
static void Calc_OLPpo(double ******p){assert(grid_alive);position_calls++;for(int d=0;d<3;d++)for(int o=0;o<2;o++)for(int a=0;a<=1;a++)for(int h=0;h<=(a?1:0);h++)for(int i=0;i<(a?2:1);i++)for(int j=0;j<(a?2:1);j++)p[d][o][a][h][i][j]=100*d+10*o+3*a+h+i+.1*j;}
static void Calc_OLPmo(double *****p){assert(grid_alive);momentum_calls++;for(int d=0;d<3;d++)for(int a=0;a<=1;a++)for(int h=0;h<=(a?1:0);h++)for(int i=0;i<(a?2:1);i++)for(int j=0;j<(a?2:1);j++)p[d][a][h][i][j]=200*d+3*a+h+i+.1*j;}
#include "cpu_scf_operator_cache.inc"
int main(void){
 double ******p=NULL;double *****m=NULL;
 assert(!CpuSCFFileCacheTake(&p,&m));
 CpuSCFFileCacheOperators();assert(position_calls==1&&momentum_calls==1);
 grid_alive=0;assert(CpuSCFFileCacheTake(&p,&m));assert(!CpuSCFFileCacheTake(&p,&m));
 assert(p[2][1][1][1][1][1]==215.1);assert(m[2][1][1][1][1]==405.1);
 for(int d=0;d<3;d++){for(int o=0;o<2;o++){for(int a=0;a<=1;a++){for(int h=0;h<=(a?1:0);h++){for(int i=0;i<(a?2:1);i++)free(p[d][o][a][h][i]);free(p[d][o][a][h]);}free(p[d][o][a]);}free(p[d][o]);}free(p[d]);for(int a=0;a<=1;a++){for(int h=0;h<=(a?1:0);h++){for(int i=0;i<(a?2:1);i++)free(m[d][a][h][i]);free(m[d][a][h]);}free(m[d][a]);}free(m[d]);}free(p);free(m);
 HS_fileout=0;CpuSCFFileCacheOperators();assert(position_calls==1&&momentum_calls==1);
 puts("operator cache: native shape, single computation, lifetime transfer, export-off PASS");return 0;
}
