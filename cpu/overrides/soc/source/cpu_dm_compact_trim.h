#ifndef CPU_DM_COMPACT_TRIM_H
#define CPU_DM_COMPACT_TRIM_H
#include <stddef.h>
#include "cpu_dm_compact.h"
/* No environment parsing or MPI in the kernel. Collective startup validation
 * and these setters are the caller's responsibility. Both defaults are zero. */
void CpuDMCompactSetTrim(int enabled);
void CpuDMCompactSetProfile(int enabled);
typedef struct {
 int logical_bands,effective_bands,nonzero_weights;
 size_t pack_calls,gemm_calls,scratch_bytes;
 double kernel_seconds,pack_seconds,gemm_seconds,scatter_seconds;
} CpuDMCompactProfile;
/* Last invocation on this process; valid when profile was enabled. */
void CpuDMCompactGetProfile(CpuDMCompactProfile *out);
#endif
