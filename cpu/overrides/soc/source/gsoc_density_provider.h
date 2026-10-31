#ifndef OPENMX_GSOC_DENSITY_PROVIDER_H
#define OPENMX_GSOC_DENSITY_PROVIDER_H
/* Coefficient-provider adapter for the existing owner-tiled BLAS3 DM/EDM
 * consumer. It changes neither weights nor spin/output conventions and owns
 * no vectors, caches or MPI collectives.
 *
 * Contract shared by C3 and GSOC providers:
 *   LocalBands: all logical physical states held by this rank, including TR
 *               partners; local index is zero based.
 *   SortedBand: one-based global physical state for native FermiFunc_NC.
 *   Pack:       contiguous zero-based full-spinor rows, local band interval;
 *               interleaved FP64 complex [column*ld+row], ld >= rows.
 * Providers must remain alive through the last requested DM/EDM component.
 *
 * GSOC_DENSITY_REMAP_C3_API is an opt-in compatibility shim for CONSUMERS
 * only (Cluster_DFT_NonCol.c and cpu_dm_compact_trim.c). Never set it when
 * compiling cpu_c3.c or cpu_gsoc.c: those implement the provider symbols.
 * No whole-matrix DM/EDM allocation is introduced by this adapter.
 */
#include <stdio.h>
#include <stdlib.h>

int CpuC3_CompactEnabled(void);
int CpuC3_CompactLocalBands(void);
int CpuC3_CompactSortedBand(int local_band);
int CpuC3_CompactPack(int row0,int rows,int band0,int bands,void *out,int ld);
int CpuGSOC_CompactEnabled(void);
int CpuGSOC_CompactLocalBands(void);
int CpuGSOC_CompactSortedBand(int local_band);
int CpuGSOC_CompactPack(int row0,int rows,int band0,int bands,void *out,int ld);

static inline int CpuDMCoefficient_Kind(void)
{
  int gsoc=CpuGSOC_CompactEnabled(),c3=CpuC3_CompactEnabled();
  if(gsoc&&c3){
    fprintf(stderr,"GSOC density provider: simultaneous C3 and GSOC compact providers are unsupported\n");
    abort();
  }
  return gsoc?2:(c3?1:0);
}
static inline int CpuDMCoefficient_Enabled(void)
{
  return CpuDMCoefficient_Kind()!=0;
}
static inline int CpuDMCoefficient_LocalBands(void)
{
  int kind=CpuDMCoefficient_Kind();
  return kind==2?CpuGSOC_CompactLocalBands():(kind==1?CpuC3_CompactLocalBands():0);
}
static inline int CpuDMCoefficient_SortedBand(int local_band)
{
  int kind=CpuDMCoefficient_Kind();
  return kind==2?CpuGSOC_CompactSortedBand(local_band):(kind==1?CpuC3_CompactSortedBand(local_band):-1);
}
static inline int CpuDMCoefficient_Pack(int row0,int rows,int band0,int bands,void *out,int ld)
{
  int kind=CpuDMCoefficient_Kind();
  return kind==2?CpuGSOC_CompactPack(row0,rows,band0,bands,out,ld):
         (kind==1?CpuC3_CompactPack(row0,rows,band0,bands,out,ld):-1);
}

#ifdef GSOC_DENSITY_REMAP_C3_API
#define CpuC3_CompactEnabled CpuDMCoefficient_Enabled
#define CpuC3_CompactLocalBands CpuDMCoefficient_LocalBands
#define CpuC3_CompactSortedBand CpuDMCoefficient_SortedBand
#define CpuC3_CompactPack CpuDMCoefficient_Pack
#endif
#endif
