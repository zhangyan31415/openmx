#!/usr/bin/env python3
"""Full-application Force3 regression against a saved CPU MoTe2 SOC reference.

Run inside a GPU allocation after sourcing the intended compiler environment.
The input must export mote2_<atoms>.out and .scfout; it is copied unchanged.
This is a short SCF correctness check, not a convergence or performance test.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import socket
import subprocess
import time

FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"


def numbers(folder, atoms):
    text = (folder / f"mote2_{atoms}.out").read_text()
    energy = float(re.findall(r"\bUtot\.?\s*=*\s*(" + FLOAT + r")", text)[-1])
    section = text.split("Eigenvalues (Hartree) for non-collinear KS-eq.", 1)[1]
    eig = []
    for line in section.splitlines():
        fields = line.split()
        if len(fields) == 2 and fields[0].isdigit() and int(fields[0]) == len(eig) + 1:
            eig.append(float(fields[1]))
        elif eig and line.strip():
            break
    section = text.split("<coordinates.forces", 1)[1].split("coordinates.forces>", 1)[0]
    forces = []
    for line in section.splitlines():
        fields = line.split()
        if len(fields) == 8 and fields[0].isdigit() and fields[1] in ("Mo", "Te"):
            forces.extend(map(float, fields[-3:]))
    return energy, eig, forces


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--binary", required=True, type=Path)
    ap.add_argument("--input", required=True, type=Path)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--work", required=True, type=Path)
    ap.add_argument("--atoms", type=int, choices=(24, 96), required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--force3-gpu", choices=("0", "1"), default="0")
    ap.add_argument("--all-force-cpu", action="store_true")
    args = ap.parse_args()
    if args.threads < 1:
        ap.error("threads must be positive")
    binary, inp, reference, work = [p.resolve() for p in (
        args.binary, args.input, args.reference, args.work)]
    assert binary.is_file() and inp.is_file()
    expected = numbers(reference, args.atoms)
    assert expected[1] and len(expected[2]) == 3 * args.atoms
    assert all(map(math.isfinite, [expected[0]] + expected[1] + expected[2]))
    work.mkdir(parents=True, exist_ok=False)
    (work / "input.dat").write_bytes(inp.read_bytes())
    env = dict(os.environ, OMP_NUM_THREADS=str(args.threads), MKL_NUM_THREADS="1",
               OPENBLAS_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1", OMP_PROC_BIND="close",
               OMP_PLACES="cores", OMP_DISPLAY_ENV="VERBOSE", OPENMX_GPU="1",
               OPENMX_FORCE3_GPU=args.force3_gpu, OPENMX_FORCE_PROFILE="1",
               NVCOMPILER_ACC_NOTIFY="1")
    if args.all_force_cpu:
        assert args.force3_gpu == "0"
        env.update(OPENMX_FORCE_GPU="0", OPENMX_FORCE4B_GPU="0", OPENMX_FORCEHNL_GPU="0")
    mpi = Path(env["OPENMX_MPI_ROOT"]) / "bin/.bin/orterun"
    # Preserve Slurm's allowed CPU set; OpenMP places/binding select its cores.
    # OpenMPI's topology can otherwise select CPUs excluded by the job cgroup.
    cmd = [str(mpi), "--mca", "pml", "ob1", "--mca", "btl", "self,vader,tcp", "-np", "1",
           "--bind-to", "none", "--report-bindings",
           str(binary), "./input.dat", "-nt", str(args.threads)]
    metadata = dict(node=socket.gethostname(), allocation=env.get("SLURM_JOB_ID"),
                    environment=env.get("OPENMX_GPU_ENV"), runtime_prefix="snapshots/home/mambaforge",
                    command=cmd, binary=str(binary), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                    source_input=str(inp), input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),
                    reference=str(reference), atoms=args.atoms, MPI=1, OMP=args.threads,
                    force3_gpu=args.force3_gpu, all_force_cpu=args.all_force_cpu,
                    CUDA_VISIBLE_DEVICES=env.get("CUDA_VISIBLE_DEVICES"),
                    cpu_affinity=sorted(os.sched_getaffinity(0)),
                    test_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    scope="short-SCF MoTe2 SOC Gamma smoke; not production convergence")
    (work / "launch.json").write_text(json.dumps(metadata, indent=2) + "\n")
    start = time.monotonic()
    loaded, gpu_samples, thread_samples = {}, [], []
    with (work / "stdout.log").open("w") as stdout, (work / "stderr.log").open("w") as stderr:
        child = subprocess.Popen(cmd, cwd=work, env=env, stdin=subprocess.DEVNULL,
                                 stdout=stdout, stderr=stderr, start_new_session=True)
        while child.poll() is None:
            # Record only this run's OpenMX process and runtime libraries.
            for entry in Path("/proc").iterdir():
                if not entry.name.isdigit():
                    continue
                try:
                    if (entry / "exe").resolve() != binary or (entry / "cwd").resolve() != work:
                        continue
                    maps = (entry / "maps").read_text()
                    loaded[entry.name] = sorted({line.split()[-1] for line in maps.splitlines()
                        if "/" in line and any(s in line for s in ("libmpi", "libnvomp", "libacc", "libcudart"))})
                    masks = {}
                    for task in (entry / "task").iterdir():
                        try:
                            status = (task / "status").read_text()
                            masks[task.name] = re.search(r"^Cpus_allowed_list:\s*(.*)$", status, re.M)[1]
                        except (OSError, TypeError):
                            pass
                    thread_samples.append(dict(seconds=time.monotonic()-start, pid=entry.name,
                        threads=len(masks), thread_cpu_masks=masks))
                except (OSError, RuntimeError):
                    pass
            query = subprocess.run(["nvidia-smi", "--query-gpu=uuid,memory.used,utilization.gpu",
                                    "--format=csv,noheader"], text=True, capture_output=True)
            gpu_samples.append(dict(seconds=time.monotonic()-start, gpu=query.stdout.strip()))
            if time.monotonic() - start > 900:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                break
            time.sleep(1)
        code = child.wait()
    elapsed = time.monotonic() - start
    out, err = [(work / name).read_text(errors="replace") for name in ("stdout.log", "stderr.log")]
    teams = re.findall(r"The number of threads in each node for OpenMP parallelization is (\d+)", out)
    checks = dict(exit_zero=code == 0, normal_finish="The calculation was normally finished." in out,
                  requested_threads_reported=teams == [str(args.threads)],
                  force3_finished=bool(re.search(r"\[Force profile\]\s+force#3\s+max=", out)))
    errors = {}
    if checks["normal_finish"]:
        energy, eig, force = numbers(work, args.atoms)
        checks.update(finite=all(map(math.isfinite, [energy] + eig + force)),
                      eigen_count=len(eig) == len(expected[1]),
                      force_count=len(force) == len(expected[2]) == 3 * args.atoms)
        errors = dict(energy_Ha=abs(energy-expected[0]),
                      eigen_Ha=max((abs(x-y) for x, y in zip(eig, expected[1])), default=math.inf),
                      force_Ha_per_bohr=max((abs(x-y) for x, y in zip(force, expected[2])), default=math.inf))
        checks.update({k: value <= 1e-7 for k, value in errors.items()})
        hs = work / f"mote2_{args.atoms}.scfout"
        checks["HS_export_exists"] = hs.is_file() and hs.stat().st_size > 1000
    force_kernels = []
    for line in err.splitlines():
        source = re.search(r"\bfile=(\S+)", line)
        if "launch CUDA kernel" in line and source and Path(source[1]).name == "Force.c":
            force_kernels.append(line)
    if args.all_force_cpu:
        checks["no_force_device_kernel"] = not force_kernels
    elif args.force3_gpu == "1":
        checks["force3_device_kernel"] = "Force3_GpuTrace_Part" in err
    else:
        checks["no_force3_device_kernel"] = not any(
            "function=Force3_GpuTrace_Part" in line for line in force_kernels)
    result = dict(**metadata, exit_code=code, seconds=elapsed, checks=checks, errors=errors,
                  passed=all(checks.values()), force_kernel_lines=force_kernels)
    for name, value in (("result.json", result), ("loaded_runtime.json", loaded),
                        ("gpu_samples.json", gpu_samples), ("thread_samples.json", thread_samples)):
        (work / name).write_text(json.dumps(value, indent=2) + "\n")
    print(work.name, "PASS" if result["passed"] else "FAIL", "exit", code,
          "seconds", round(elapsed, 3), errors, flush=True)
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
