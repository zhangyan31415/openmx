# Source-only CPU dependencies

`openmx4.0/source` is a source/header/template snapshot of the installed CPU OpenMX4.0 tree at `snapshots/openmx4.0/source`, including bundled `elpa-2018.05.001`. It has not been independently established to be a pristine upstream tarball. Source copyright and license notices are retained. No objects, compiler module files, executables, or prebuilt libraries are vendored.

`ORIGINS.json` records origin labels, exact SHA256 values, and source-only reconstruction provenance for optimized files under `../overrides/soc`. Personal paths have been replaced by snapshot-relative labels; see [public provenance rules](../docs/PROVENANCE_PATHS.md). These labels are provenance, never build inputs. The HNL batch Force translation was recovered with the historical source-only transformation and verified against its saved translation hash.

The CPU target selection lives in `../overrides/soc/SOURCES.json`; only those translation units are compiled. The other stock source/header files provide the original package context and include closure. Declared external dependencies are Intel classic C/Fortran compiler, Intel MPI, and MKL/FFTW compatibility interfaces. See `../docs/BUILD.zh-CN.md` for the supported build command and profile limitations.
