# Public provenance and local path redaction

Published source manifests and saved validation records omit personal filesystem
paths. Source content, numerical results, original artifact hashes, software
versions, and scientific authorship remain recorded.

The following names are provenance labels, not dependencies that a user must
create or paths that the build reads:

- `snapshots/results/...`: archived experiment and verification artifacts.
- `snapshots/openmx4.0/source/...`: the original installed source snapshot.
- `snapshots/project/...` and `snapshots/moirekp-release/...`: historical source
  locations outside this standalone CPU package.
- `snapshots/scratch/...`: historical temporary build artifacts.
- `repo/...`: a path within the repository at the time of a saved command.

An existing source or artifact SHA256 continues to identify the original bytes
described by that record. It is not a checksum of a public JSON record whose path
strings have been redacted. Where a distributed source file changes, its original
snapshot hash is retained separately from the current distributed file hash.
Archived command records containing these labels are evidence, not runnable
configuration files. Use `tools/examples/` and the current workflow instructions
to prepare a new calculation.

User instructions use `$CONDA_PREFIX` for an activated Python environment and
`/path/to/python-env` where an explicit prefix is required in JSON. The build
uses `ONEAPI_ROOT` when set, otherwise `/opt/intel/oneapi`; `--oneapi-prefix`
overrides both. The supported classic Intel compiler, MPI and MKL versioned
subdirectory layout is documented in `BUILD.zh-CN.md`. No shell-variable
expansion is performed inside JSON configuration values.
