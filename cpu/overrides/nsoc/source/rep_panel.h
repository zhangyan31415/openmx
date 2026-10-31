#ifndef OPENMX_REP_PANEL_H
#define OPENMX_REP_PANEL_H
#include <stddef.h>
#include <stdint.h>
#include <mpi.h>
typedef struct RpaPlan RpaPlan;
typedef struct {int row,rows,col,cols,multiplicity;size_t offset;} RpaBlock;
typedef struct {uint64_t geometry,basis,distribution,projector_map;} RpaContext;
typedef struct {RpaContext context;uint64_t revision,value_hash;int kind,step;} RpaStamp;
typedef struct {
 int n,natoms,species,ng,sectors[4];const int64_t *cp;const int *ri;const double *cv;
 const int *types,*norb;const int *const *permutations;
 const double *const *rotations;RpaContext context;
} RpaBasis;
typedef struct {size_t retained_bytes,peak_bytes,local_representatives,
 endpoint_count,panels,routes,send_values,recv_values,max_panel_values;
 double intertwiner_residual;int value_collectives;} RpaStats;
/* blocks and packed values are local representative-owner data. Actual
 * periodic orbit validation belongs to the trusted RS metadata provider. */
RpaPlan *Rpa_Create(MPI_Comm comm,const RpaBasis *basis,int blocks,
 const RpaBlock *block,int block_size,int grid_rows,int grid_cols,
 const int *actual_rank_coordinates,size_t cap,char *error,size_t error_size);
uint64_t Rpa_ValueHash(const double *values,size_t count);
size_t Rpa_InputSize(const RpaPlan *p);
size_t Rpa_LocalSize(const RpaPlan *p,int sector); /* E,A+,A- */
int Rpa_Apply(RpaPlan *p,const RpaStamp *stamp,const double *values,
 double *const *out,const size_t *capacities,char *error,size_t error_size);
void Rpa_GetStats(const RpaPlan *p,RpaStats *stats);
void Rpa_Destroy(RpaPlan *p);
#endif
