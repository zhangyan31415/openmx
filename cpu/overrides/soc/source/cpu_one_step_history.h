#ifndef CPU_ONE_STEP_HISTORY_H
#define CPU_ONE_STEP_HISTORY_H

#include <string.h>

/* Capacity planning only: the input history, mixing weights, and convergence
 * controls are deliberately not modified. Keep three H slots as well as three
 * DM slots: before restart success is known, first-step Simple_Mixing_H can
 * still read/write H history slots 0 and 1. */
typedef struct {
  int md_switch, md_iter, solver, spin, soc, contraction;
  int max_scf, mixing, control_temp, pulay_start, symmetry_mode;
  int dm_capacity, h_capacity;
} CpuOneStepHistoryInput;

static int CpuOneStepHistoryPlan(const char *flag,
                                const CpuOneStepHistoryInput *input,
                                int *dm_capacity, int *h_capacity)
{
  *dm_capacity = input->dm_capacity;
  *h_capacity = input->h_capacity;
  if (!flag || !strcmp(flag, "0")) return 0;
  if (strcmp(flag, "1")) return -1;
  if (input->md_switch != 0 || input->md_iter != 1 || input->solver != 2 ||
      input->spin != 3 || input->soc != 1 || input->contraction != 0 ||
      input->max_scf != 1 || input->mixing != 1 || input->control_temp != 0 ||
      input->pulay_start < 4 || (input->symmetry_mode != 0 && input->symmetry_mode != 3) ||
      input->dm_capacity < 3 || input->h_capacity < 3) return -2;
  *dm_capacity = 3;
  *h_capacity = 3;
  return 1;
}

#endif
