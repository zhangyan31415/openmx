#include "openmx_nc_local_dm_gpu.h"
#include <cuda_runtime_api.h>
#include <openacc.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static char dm_error[256];

const char *openmx_nc_local_dm_gpu_error(void)
{
    return dm_error;
}

static int decline(const char *message)
{
    snprintf(dm_error, sizeof(dm_error), "%s", message);
    return 0;
}

int openmx_nc_local_dm_gpu(int n, int states, int entries, int flag,
                          const double *vectors, const double *weights,
                          const int *basis0, const int *basis1, double *output)
{
    double *dv = NULL, *dw = NULL, *dr = NULL;
    int *da = NULL, *db = NULL;
    int status = 0, device = 0;
    dm_error[0] = '\0';
    if (n <= 0 || n > INT_MAX/2 || states < 0 || entries < 0)
        return decline("invalid local DM dimensions");
    if (flag < 1 || flag > 10)
        return decline("unsupported local DM component");
    const int n2 = 2*n;
    const size_t state_bytes = (size_t)n2 * 2 * sizeof(double);
    if ((size_t)states > SIZE_MAX/state_bytes)
        return decline("local spinor size overflow");
    if (entries == 0) return 1;
    if (output == NULL) return decline("missing local DM output");
    if (states == 0) {
        memset(output, 0, (size_t)entries*sizeof(double));
        return 1;
    }
    if (vectors == NULL || weights == NULL || basis0 == NULL || basis1 == NULL)
        return decline("missing local DM input");
    for (int p = 0; p < entries; ++p)
        if (basis0[p] < 0 || basis0[p] >= n || basis1[p] < 0 || basis1[p] >= n)
            return decline("sparse DM orbital index out of range");
    int nonzero = 0;
    for (int k = 0; k < states; ++k) {
        if (!isfinite(weights[k])) return decline("nonfinite local DM weight");
        if (weights[k] != 0.0) nonzero = 1;
    }
    if (!nonzero) {
        memset(output, 0, (size_t)entries*sizeof(double));
        return 1;
    }
    cudaError_t initial = cudaGetDevice(&device);
    if (initial != cudaSuccess) {
        snprintf(dm_error, sizeof(dm_error), "CUDA unavailable: %s", cudaGetErrorString(initial));
        (void)cudaGetLastError();
        return 0;
    }
    acc_set_device_num(device, acc_device_nvidia);

#define DM_CUDA(expression, failure_status) do { \
    cudaError_t dm_status = (expression); \
    if (dm_status != cudaSuccess) { \
        snprintf(dm_error, sizeof(dm_error), "%s: %s", #expression, cudaGetErrorString(dm_status)); \
        status = (failure_status); \
        goto cleanup; \
    } \
} while (0)

    const int capacity = entries < 65536 ? entries : 65536;
    DM_CUDA(cudaMalloc((void **)&dv, (size_t)states*state_bytes), 0);
    DM_CUDA(cudaMalloc((void **)&dw, (size_t)states*sizeof(double)), 0);
    DM_CUDA(cudaMalloc((void **)&da, (size_t)capacity*sizeof(int)), 0);
    DM_CUDA(cudaMalloc((void **)&db, (size_t)capacity*sizeof(int)), 0);
    DM_CUDA(cudaMalloc((void **)&dr, (size_t)capacity*sizeof(double)), 0);
    DM_CUDA(cudaMemcpy(dv, vectors, (size_t)states*state_bytes, cudaMemcpyHostToDevice), -1);
    DM_CUDA(cudaMemcpy(dw, weights, (size_t)states*sizeof(double), cudaMemcpyHostToDevice), -1);

    const int a_offset = (flag == 2 || flag == 6 || flag == 8) ? n : 0;
    const int b_offset = (flag == 1 || flag == 5 || flag == 7) ? 0 : n;
    const int imaginary = (flag == 4 || flag == 5 || flag == 6 || flag == 10);
    for (int offset = 0; offset < entries; ) {
        const int count = entries-offset < capacity ? entries-offset : capacity;
        DM_CUDA(cudaMemcpy(da, basis0+offset, (size_t)count*sizeof(int), cudaMemcpyHostToDevice), -1);
        DM_CUDA(cudaMemcpy(db, basis1+offset, (size_t)count*sizeof(int), cudaMemcpyHostToDevice), -1);
#pragma acc parallel loop gang vector vector_length(128) deviceptr(dv,dw,da,db,dr)
        for (int p = 0; p < count; ++p) {
            const int ia = da[p] + a_offset;
            const int ib = db[p] + b_offset;
            double sum = 0.0;
            if (imaginary) {
#pragma acc loop seq
                for (int k = 0; k < states; ++k) {
                    const double w = dw[k];
                    if (w != 0.0) {
                        const size_t ai = 2*((size_t)k*n2 + ia);
                        const size_t bi = 2*((size_t)k*n2 + ib);
                        sum += w*(dv[ai]*dv[bi+1] - dv[ai+1]*dv[bi]);
                    }
                }
            }
            else {
#pragma acc loop seq
                for (int k = 0; k < states; ++k) {
                    const double w = dw[k];
                    if (w != 0.0) {
                        const size_t ai = 2*((size_t)k*n2 + ia);
                        const size_t bi = 2*((size_t)k*n2 + ib);
                        sum += w*(dv[ai]*dv[bi] + dv[ai+1]*dv[bi+1]);
                    }
                }
            }
            dr[p] = sum;
        }
        DM_CUDA(cudaGetLastError(), -1);
        DM_CUDA(cudaDeviceSynchronize(), -1);
        DM_CUDA(cudaMemcpy(output+offset, dr, (size_t)count*sizeof(double), cudaMemcpyDeviceToHost), -1);
        offset += count;
    }
    status = 1;
cleanup:
    if (dr != NULL) (void)cudaFree(dr);
    if (db != NULL) (void)cudaFree(db);
    if (da != NULL) (void)cudaFree(da);
    if (dw != NULL) (void)cudaFree(dw);
    if (dv != NULL) (void)cudaFree(dv);
    (void)cudaGetLastError();
    return status;
#undef DM_CUDA
}
