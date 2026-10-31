#ifndef CPU_MEMORY_INVENTORY_H
#define CPU_MEMORY_INVENTORY_H

/* Caller includes openmx_common.h and malloc.h. This is diagnostic-only and
 * must run while the native history trees and geometry maps are still alive.
 * Slots: DM, iDM, ResidualDM, iResidualDM, HisH1, HisH2.
 *
 * malloc_usable_size includes allocator padding but is not resident memory.
 * usable_bytes + 8 * allocation_count is a glibc chunk-footprint estimate,
 * not a PSS measurement (large mmap chunks/arena fragmentation differ). */
static int CpuMemoryInventoryAdd(const void *pointer, size_t minimum,
                                unsigned long long *bytes,
                                unsigned long long *count)
{
  size_t usable;
  if (!pointer) return 0;
  usable = malloc_usable_size((void *)pointer);
  if (usable < minimum) return 0;
  *bytes += (unsigned long long)usable;
  ++*count;
  return 1;
}

static int CpuMemoryHistoryTree(double ******tree, int histories, int channels,
                               const int *atom_map,
                               unsigned long long *bytes,
                               unsigned long long *count)
{
  int m, k, a, h, i;
  if (!tree || histories <= 0 || channels <= 0 || !atom_map ||
      !CpuMemoryInventoryAdd(tree, (size_t)histories * sizeof(*tree),
                             bytes, count)) return 0;
  for (m = 0; m < histories; ++m) {
    if (!CpuMemoryInventoryAdd(tree[m], (size_t)channels * sizeof(*tree[m]),
                               bytes, count)) return 0;
    for (k = 0; k < channels; ++k) {
      if (!CpuMemoryInventoryAdd(tree[m][k],
              ((size_t)Matomnum + 1) * sizeof(*tree[m][k]), bytes, count)) return 0;
      for (a = 0; a <= Matomnum; ++a) {
        int global_atom = a ? atom_map[a] : 0;
        int neighbors = 0, rows = 1;
        if (a) {
          int species;
          if (global_atom < 1 || global_atom > atomnum) return 0;
          species = WhatSpecies[global_atom];
          if (species < 0 || species >= SpeciesNum) return 0;
          neighbors = FNAN[global_atom];
          rows = Spe_Total_NO[species];
          if (neighbors < 0 || rows < 1 || !natn[global_atom]) return 0;
        }
        if (!CpuMemoryInventoryAdd(tree[m][k][a],
                ((size_t)neighbors + 1) * sizeof(*tree[m][k][a]),
                bytes, count)) return 0;
        for (h = 0; h <= neighbors; ++h) {
          int columns = 1;
          if (a) {
            int neighbor_atom = natn[global_atom][h], neighbor_species;
            if (neighbor_atom < 1 || neighbor_atom > atomnum) return 0;
            neighbor_species = WhatSpecies[neighbor_atom];
            if (neighbor_species < 0 || neighbor_species >= SpeciesNum) return 0;
            columns = Spe_Total_NO[neighbor_species];
            if (columns < 1) return 0;
          }
          if (!CpuMemoryInventoryAdd(tree[m][k][a][h],
                  (size_t)rows * sizeof(*tree[m][k][a][h]), bytes, count)) return 0;
          for (i = 0; i < rows; ++i)
            if (!CpuMemoryInventoryAdd(tree[m][k][a][h][i],
                    (size_t)columns * sizeof(*tree[m][k][a][h][i]),
                    bytes, count)) return 0;
        }
      }
    }
  }
  return 1;
}

static void CpuMemoryHistoryInventory(unsigned long long usable_bytes[6],
                                     unsigned long long allocation_count[6])
{
  unsigned long long bytes[6] = {0,0,0,0,0,0};
  unsigned long long count[6] = {0,0,0,0,0,0};
  int j;
  for (j = 0; j < 6; ++j) usable_bytes[j] = allocation_count[j] = 0;
  /* In other mixing/spin modes some native globals may point to dummy trees
   * with different dimensions; never attempt to traverse those layouts. */
  if (Mixing_switch != 1 || SpinP_switch != 3 || SO_switch != 1 || Cnt_switch ||
      Matomnum < 0 || atomnum < 1 || SpeciesNum < 1 ||
      List_YOUSO[16] < 1 || List_YOUSO[39] < 1 || !M2G || !S_M2G ||
      !FNAN || !WhatSpecies || !Spe_Total_NO || !natn ||
      !DM || !iDM || !ResidualDM || !iResidualDM || !HisH1 || !HisH2) return;
  if (!CpuMemoryHistoryTree(DM, List_YOUSO[16], 4, M2G, bytes+0, count+0) ||
      !CpuMemoryHistoryTree(iDM, List_YOUSO[16], 2, M2G, bytes+1, count+1) ||
      !CpuMemoryHistoryTree(ResidualDM, List_YOUSO[16], 4, M2G, bytes+2, count+2) ||
      !CpuMemoryHistoryTree(iResidualDM, List_YOUSO[16], 2, M2G, bytes+3, count+3) ||
      !CpuMemoryHistoryTree(HisH1, List_YOUSO[39], 4, S_M2G, bytes+4, count+4) ||
      !CpuMemoryHistoryTree(HisH2, List_YOUSO[39], 3, S_M2G, bytes+5, count+5)) return;
  for (j = 0; j < 6; ++j) {
    usable_bytes[j] = bytes[j];
    allocation_count[j] = count[j];
  }
}

#endif
