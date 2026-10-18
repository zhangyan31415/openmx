#ifndef OPENMX_NC_LOCAL_DM_GPU_H
#define OPENMX_NC_LOCAL_DM_GPU_H

#ifdef __cplusplus
extern "C" {
#endif

/* State-major local spinors, interleaved real/imaginary doubles. Weights
   already include the occupation cutoff and, for EDM, the eigenvalue.
   Sparse indices address spatial orbitals [0,n). No global dense matrix is
   formed. The caller retains MPI reduction and atom-block storage.
   Return 1: computed (including trivial zero inputs); 0: safely declined;
   -1: device execution/transfer failed. On failure output may be partial.
   Called by the MPI master thread only. Allocations are released per call. */
int openmx_nc_local_dm_gpu(int n, int states, int entries, int flag,
                          const double *vectors, const double *weights,
                          const int *basis0, const int *basis1, double *output);
const char *openmx_nc_local_dm_gpu_error(void);

#ifdef __cplusplus
}
#endif
#endif
