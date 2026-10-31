#ifndef CPU_NSOC_REP_TAIL_H
#define CPU_NSOC_REP_TAIL_H
#include <mpi.h>
typedef struct {int admitted,cut_states,covered_states,invalid;double density_bound,energy_density_bound,amplification,tolerance;}CpuRepTailEvidence;
/* C is native AO-major, local-state contiguous; f and energy are actual local
 * weights, not estimates. Cut IDs are 1-based global returned-state indices. */
int CpuRepTailAdmit(MPI_Comm comm,int nao,int first,int last,const double*coeff,
 const double*f,const double*energy,int cut_count,const int*cut_ids,
 double amplification,double tolerance,CpuRepTailEvidence*out);
/* Upper bound for actual stored AO block transforms, including finite
 * row-sum/product roundoff; independent of ideal orthogonality. */
double CpuRepTailAmplification(int operations,int species,const int*norb,
 const double*const*rotations,double*max_row_l1);
#endif
