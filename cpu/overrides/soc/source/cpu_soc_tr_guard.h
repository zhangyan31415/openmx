#ifndef CPU_SOC_TR_GUARD_H
#define CPU_SOC_TR_GUARD_H
/* Pure parser: unset/0 -> 0, 1 -> 1, anything else -> -1. No MPI call. */
int CpuSOCNative_TRRequested(void);
/* Call collectively at DFT's common entry. The MarkSolve hook also validates
 * collectively so enabling the guard cannot silently bypass native export. */
void CpuSOCNative_TRValidateConfig(void);
/* A configuration request alone is never proof of actual operator symmetry. */
int CpuSOCNative_TRPassed(void);
#endif
