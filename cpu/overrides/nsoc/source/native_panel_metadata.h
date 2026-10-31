#ifndef OPENMX_NATIVE_PANEL_METADATA_H
#define OPENMX_NATIVE_PANEL_METADATA_H
#include "rep_panel.h"
/* Set collectively once before BeforeSolve(step1), default0. */
int CpuRSym_PanelEnable(int enabled);
/* Borrowed immutable metadata, valid while the fixed native geometry lives. */
int CpuRSym_PanelAOMetadata(int *nat,int *species,int *ng,
 const int **types,const int **norb,const int *const **permutations,
 const double *const **rotations);
/* Local original-owner representatives, same order as rp_offsets, caller free. */
RpaBlock *CpuRSym_PanelBlocks(int *count);
RpaContext CpuRSym_PanelContext(uint64_t basis,uint64_t distribution);
/* Reads only the field sealed immediately by BeforeSolve. A changed pointer,
 * value, geometry, step or field type is refused; no fresh stamp is invented. */
int CpuRSym_PanelPack(int kind,int step,RpaContext context,double *values,
 size_t entries,RpaStamp *stamp);
int CpuRSym_PanelFields(double ****h,double ****s);
#endif
