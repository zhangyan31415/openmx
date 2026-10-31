#ifndef CPU_HASSEMBLY_WORKSPACE_H
#define CPU_HASSEMBLY_WORKSPACE_H
#include <stddef.h>
int CpuHAssembly_Compact(void);
void CpuHAssembly_Validate(void);
/* Borrow storage only for this synchronous assembly. No cached H values,
 * ownership transfer, basis or communicator assumptions are introduced. */
void CpuHAssembly_WithScratch(double ****RH,double *Hs,int *MP,int spin,int myworld1,double *workspace,size_t capacity);
#endif
