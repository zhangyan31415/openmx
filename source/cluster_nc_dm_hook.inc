/* Rank-local sparse DM/EDM: retain the established band distribution and
   host MPI reduction. A root-dense eigenvector copy is never required. */
static int ClusterNonCol_TryLocalGpuDM(int flag, int myid, int entries,
                                      const int *is2, const int *ie2, int *MP,
                                      int n, int n2, const double *ko,
                                      const dcomplex *vectors, double *out)
{
    const char *setting = getenv("OPENMX_DISTRIBUTED_DM_GPU");
    const char *gpu_request = getenv("OPENMX_GPU");
    const int mode = setting ? atoi(setting) : 1;
    if (mode == 0 || (gpu_request && atoi(gpu_request) == 0) ||
        scf_eigen_lib_flag != GPUSOLVER || flag < 1 || flag > 10 ||
        n2 < Cluster_DFT_NonCol_GpuSwitchNum() || entries <= 0 || entries == INT_MAX)
        return 0;
    _Static_assert(sizeof(dcomplex) == 2*sizeof(double), "unexpected complex layout");
    ClusterNonCol_DMEntryCache_Ensure(MP, n);
    ClusterNonColDMEntryCache *cache = &ClusterNonCol_dm_entry_cache;
    if (cache->entry_count != entries)
        ClusterNonCol_AbortWithMessage("Rank-local GPU DM sparse entry count mismatch");
    int states = ie2[myid] - is2[myid] + 1;
    if (states < 0) states = 0;
    double *weights = malloc((size_t)(states ? states : 1)*sizeof(double));
    int active = 0, result = 0;
    const char *reason = "host occupation allocation failed";
    double start, end;
    dtime(&start);
    if (weights != NULL) {
        for (int local = 0; local < states; ++local) {
            const int k = is2[myid] + local;
            double x = (ko[k]-ChemP)*Beta;
            if (x <= -60.0) x = -60.0;
            if (x >= 60.0) x = 60.0;
            const double occupation = FermiFunc_NC(x, k);
            /* Apply exactly the legacy occupation test before multiplying
               by a possibly negative energy. */
            weights[local] = (occupation > 1.0e-13)
                ? (flag >= 7 ? occupation*ko[k] : occupation) : 0.0;
            if (weights[local] != 0.0) ++active;
        }
        result = openmx_nc_local_dm_gpu(n, states, entries, flag,
                    (const double *)vectors, weights, cache->basis0, cache->basis1, out);
        reason = openmx_nc_local_dm_gpu_error();
        free(weights);
    }
    dtime(&end);
    if (result < 0 || (mode >= 2 && result == 0)) {
        char message[512];
        snprintf(message, sizeof(message), "rank %d GPU DM flag %d failed: %s", myid, flag, reason);
        ClusterNonCol_AbortWithMessage(message);
    }
    const char *verbose = getenv("OPENMX_DIST_DM_VERBOSE");
    if (result == 0 || (verbose && atoi(verbose))) {
        int device = -1;
        (void)cudaGetDevice(&device);
        fprintf(stderr, "<DIST_DM> rank=%d flag=%d states=%d active=%d entries=%d device=%d mode=%s compute_seconds=%.9f%s%s\n",
                myid, flag, states, active, entries, device,
                result ? (active ? "GPU" : "zero") : "CPU-fallback", end-start,
                result ? "" : " reason=", result ? "" : reason);
    }
    /* A declined attempt may have filled early tiles. Reset them before
       the legacy += loop so allocation fallback cannot double count. */
    if (result == 0) memset(out, 0, (size_t)entries*sizeof(double));
    return result == 1;
}
