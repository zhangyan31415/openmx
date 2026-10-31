#ifndef OPENMX_CPU_DM_COMPACT_H
#define OPENMX_CPU_DM_COMPACT_H
/* Same block/output convention as cpu_dm_blocks. Vectors are supplied by
 * CpuC3_CompactPack; the vectors argument is unused and may be NULL. */
int cpu_dm_blocks_compact(int n,int bands,const void *vectors,const double *weights,
 int block_count,const int *blocks,int entries,int energy_density,int tile,double *output);
#endif
