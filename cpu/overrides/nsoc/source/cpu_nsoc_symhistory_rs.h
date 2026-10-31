#ifndef CPU_NSOC_SYMHISTORY_RS_H
#define CPU_NSOC_SYMHISTORY_RS_H
#include <stddef.h>
size_t CpuRSym_HistorySize(void);
const double *CpuRSym_HistoryWeights(void);
void CpuRSym_HistoryPack(double ****field,double *representative);
void CpuRSym_HistoryPackChecked(double ****field,double *representative,int step,const char *label);
void CpuRSym_HistoryExpand(const double *representative,double ****field);
void CpuRSym_HistoryProbe(int step,double ******history,double ******residual);
#endif
