/* Native grids/structures are emitted before deferred projected SOC bands. */
#define OutData CpuOriginal_OutData
#include "../../../vendor/openmx4.0/source/OutData.c"
#undef OutData
extern void CpuTargetBands_RunPending(void);
double OutData(char *inputfile){
 double native=CpuOriginal_OutData(inputfile),start=MPI_Wtime();
 CpuTargetBands_RunPending();return native+MPI_Wtime()-start;
}
