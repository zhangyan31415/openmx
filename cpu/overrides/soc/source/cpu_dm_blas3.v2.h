#ifndef OPENMX_CPU_DM_BLAS3_H
#define OPENMX_CPU_DM_BLAS3_H
/* State-major complex spinors: [band][up spatial, down spatial].
 * Each block is {row_start,row_count,col_start,col_count,output_offset}.
 * Weights are already occupation-filtered; signed weights are supported.
 * Output matches native OpenMX components, including its imaginary signs.
 * This local helper performs no MPI and does not allocate a global dense D. */
int cpu_dm_blocks(int n, int bands, const void *vectors, const double *weights,
                  int block_count, const int *blocks, int entries,
                  int energy_density, int tile, double *output);
int cpu_dm_compare(int count, const double *reference, const double *candidate, double *max_abs, double *scale);
#endif
