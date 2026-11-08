#ifndef CPU_NSOC_OWNER_H
#define CPU_NSOC_OWNER_H
#include <stddef.h>
#include <mpi.h>
typedef struct {int row,rows,col,cols,owner;} CpuNOBlock;
typedef struct CpuNOPlan CpuNOPlan;
typedef struct {void*context;int(*rows)(void*context,int first_row,int rows,double*output,int leading_dimension);} CpuNORowProvider;
CpuNOPlan *CpuNOCreate(MPI_Comm comm,int nao,int blocks,const CpuNOBlock *b,size_t tile_bytes);
int CpuNOCompute(CpuNOPlan *p,int local_states,const double *coeff,const double *occupation,const double *energy);
int CpuNOComputeProvided(CpuNOPlan*p,int local_states,const CpuNORowProvider*provider,const double*occupation,const double*energy);
size_t CpuNOEntries(const CpuNOPlan *p);
size_t CpuNOTemporaryBytes(const CpuNOPlan *p);
const double *CpuNODensity(const CpuNOPlan *p);
const double *CpuNOEnergyDensity(const CpuNOPlan *p);
void CpuNODestroy(CpuNOPlan *p);
#endif
