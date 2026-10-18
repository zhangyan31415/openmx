/* Standalone NVHPC 24.1 OpenMP private pointer-array miscompilation.
 * No OpenMX, MPI or GPU calls. On EPYC 7282, 24.1 fails even with one
 * thread and -Mnovect; 24.7 succeeds. -DUSE_SCALAR tests the workaround.
 *
 * nvc -O3 -Mnovect -fopenmp nvc_omp_pointer_array_bug.c -o repro
 * OMP_NUM_THREADS=1 ./repro 3
 * nvc -O3 -Mnovect -fopenmp -DUSE_SCALAR nvc_omp_pointer_array_bug.c -o scalar
 * OMP_NUM_THREADS=1 ./scalar 3
 *
 * GCC: use -O3 -fopenmp. Without OpenMP, all pragmas may be ignored.
 */
#include <stdio.h>
#include <stdlib.h>

#ifdef USE_SCALAR
#define D0 dst0
#define D1 dst1
#define D2 dst2
#define D3 dst3
#else
#define D0 dst[0]
#define D1 dst[1]
#define D2 dst[2]
#define D3 dst[3]
#endif

int main(int argc, char **argv)
{
    int spin = argc > 1 ? atoi(argv[1]) : 3;
    const int n = 4096;
    double *grid[4];
    if (spin != 0 && spin != 1 && spin != 3) return 2;
    for (int s = 0; s < 4; ++s) {
        grid[s] = calloc((size_t)n, sizeof(double));
        if (!grid[s]) return 2;
    }
#pragma omp parallel
    {
#ifdef USE_SCALAR
        double *dst0 = grid[0];
        double *dst1 = spin >= 1 ? grid[1] : NULL;
        double *dst2 = spin == 3 ? grid[2] : NULL;
        double *dst3 = spin == 3 ? grid[3] : NULL;
#else
        double *dst[4];
        for (int s = 0; s <= spin; ++s) dst[s] = grid[s];
#endif
#pragma omp for
        for (int k = 0; k < n; ++k) {
            if (spin == 0) {
                D0[k] = 4.0 * (k + 1);
            } else if (spin == 1) {
                D0[k] = 2.0 * (k + 1);
                D1[k] = 2.0 * (k + 2);
            } else {
                D0[k] = 2.0 * (k + 1);
                D1[k] = 2.0 * (k + 2);
                D2[k] = 4.0 * (k + 3);
                D3[k] = -4.0 * (k + 4);
            }
        }
    }
    int errors = 0;
    for (int s = 0; s < 4; ++s) {
        for (int k = 0; k < n; ++k) {
            double expected = 0.0;
            if (s <= spin) {
                double weight = spin == 0 || s == 2 ? 4.0 : (s == 3 ? -4.0 : 2.0);
                expected = weight * (k + s + 1);
            }
            errors += grid[s][k] != expected;
        }
        free(grid[s]);
    }
    printf("spin=%d errors=%d\n", spin, errors);
    return errors != 0;
}
