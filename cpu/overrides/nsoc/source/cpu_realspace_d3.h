#ifndef CPU_REALSPACE_D3_H
#define CPU_REALSPACE_D3_H
int CpuRSym_Enabled(void);
void CpuRSym_SetOccupationReference(double electrons,double band_energy);
void CpuRSym_BeforeSolve(int step, double *****hamiltonian, double ****overlap);
void CpuRSym_AfterDensity(int step, double *****density, double *****energy_density);
void CpuRSym_BeforeExport(void);
#endif
