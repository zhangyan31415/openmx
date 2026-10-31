#ifndef OPENMX_GSOC_RUNTIME_GUARD_H
#define OPENMX_GSOC_RUNTIME_GUARD_H
#include <stddef.h>
#include <stdint.h>
#include "mpi.h"
/* This guard must be called by every rank, even when its local GSOC flag is off,
 * before native dense placeholders are allocated. Restart-read success remains
 * checked at the solve after RestartFileDFT has actually read the checkpoint. */
typedef struct {
 int solver,spin,soc,md,contracted,scf_steps,c2n,restart,restart_spin;
 int hubbard,constraint,zeeman_spin,zeeman_orbital;
 int cwf,export_dm,energy_window_dm,dos,dos_gauss,optical,mo;
 int atoms;
 /* Other state-based postprocessors listed in the native dense-release guard. */
 int dense_postprocessing;
 int level_fileout,partial_charge,atomic_analysis;
} GsocRuntimeInput;
typedef struct {
 GsocRuntimeInput input;
 int enabled,compact,tile,probes,project,allow_fresh,block_buckets;
 double covariance_tol,geometry_tol;
 uint64_t environment_digest[2],basis_digest[2],basis_bytes;
 int32_t header[6],sector_size[48],sector_type[48],sector_partner[48];
 int64_t basis_nnz;
} GsocRuntimeConfig;
/* Capture never aborts or enters MPI: local errors must still participate in
 * gsoc_runtime_agree, preventing a missing file on one rank from hanging peers. */
int gsoc_runtime_capture(GsocRuntimeConfig *,const GsocRuntimeInput *,char *,size_t);
void gsoc_runtime_agree(MPI_Comm,const GsocRuntimeConfig *,int,const char *);
void CpuGSOC_ValidateRuntime(void);
#endif
