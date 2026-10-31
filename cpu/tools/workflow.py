#!/usr/bin/env python3
"""Explicit, reproducible single/multi-node NSOC -> own-checkpoint -> one SOC.

`plan` is read-only. `run` submits a supervised controller on the first named
SSH host. No scheduler operations, inherited optimization flags or input edits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import socket
import signal
import subprocess
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
THREAD_KEYS = ('OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS')
SAFE_NAME = re.compile(r'[A-Za-z0-9_.-]+\Z')
HOST_NAME = re.compile(r'[A-Za-z0-9][A-Za-z0-9_.-]*\Z')


def loads_unique(text):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('duplicate JSON key: '+key)
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError('nonfinite JSON number: '+value)
    return json.loads(text, object_pairs_hook=pairs, parse_constant=nonfinite)


def load(path):
    return loads_unique(Path(path).read_text())


def write_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name+'.writing-'+str(os.getpid()))
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    tmp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def absolute(value, label, *, exists=True, file=False):
    p = Path(value)
    if not p.is_absolute():
        raise ValueError(label+' must be an absolute path')
    p = p.resolve()
    if exists and not p.exists():
        raise ValueError(label+' does not exist: '+str(p))
    if file and not p.is_file():
        raise ValueError(label+' must be a regular file')
    return p


def positive_int(value, label):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(label+' must be a positive integer')
    return value


def input_values(text):
    values, blocks, active = {}, {}, None
    for raw in text.splitlines():
        line = raw.split('#',1)[0].strip()
        if not line:
            continue
        if line.startswith('<'):
            key = line[1:].strip().lower()
            if active is not None or key in blocks:
                raise ValueError('duplicate or nested input block: '+key)
            active = key; blocks[key] = []
        elif active is not None:
            if line.lower() == active+'>':
                active = None
            else:
                blocks[active].append(line.split())
        else:
            words = line.split(None,1)
            if len(words) != 2:
                raise ValueError('invalid input line: '+line)
            key, value = words[0].lower(), words[1]
            if key in values:
                raise ValueError('duplicate input keyword: '+key)
            values[key] = value.strip()
    if active is not None:
        raise ValueError('unterminated input block: '+active)
    return values, blocks


def _numeric(value):
    try:
        f = float(value)
    except ValueError:
        return value.lower()
    if not math.isfinite(f):
        raise ValueError('nonfinite physical input')
    return f


def structure(text):
    values, blocks = input_values(text)
    result = {'coordinate_unit':values.get('atoms.speciesandcoordinates.unit','au').lower(),
              'cell_unit':values.get('atoms.unitvectors.unit','au').lower()}
    for name in ('atoms.unitvectors','definition.of.atomic.species','atoms.speciesandcoordinates'):
        if name not in blocks:
            raise ValueError('missing structure/basis block: '+name)
        result[name] = [[_numeric(x) for x in (row[:5] if name=='atoms.speciesandcoordinates' else row)] for row in blocks[name]]
    return result


def validate_inputs(nsoc_text, soc_text):
    n, _ = input_values(nsoc_text); s, _ = input_values(soc_text)
    if structure(nsoc_text) != structure(soc_text):
        raise ValueError('NSOC/SOC structure, atomic order and basis must match')
    for v in (n,s):
        if not SAFE_NAME.fullmatch(v.get('system.name','')):
            raise ValueError('simple explicit System.Name required')
        if v.get('system.currrentdirectory','./') not in ('.','./'):
            raise ValueError('System.CurrrentDirectory must be ./ for isolated staging')
        if v.get('scf.symmetry.mode','off').lower() != 'off':
            raise ValueError('input parser Symmetry Auto/User is not wired: use explicit prepared profiles')
        if v.get('scf.eigenvaluesolver','').lower() not in ('cluster','band') or v.get('scf.kgrid','1 1 1').split()!=['1','1','1']:
            raise ValueError('only fixed nonmagnetic Gamma Cluster/Band inputs supported')
        if v.get('scf.eigen.lib','').lower() not in ('elpa1','elpa2'):
            raise ValueError('explicit elpa1/elpa2 required')
        if v.get('md.type','nomd').lower() not in ('nomd','none'):
            raise ValueError('fixed geometry/no-MD required')
        if v.get('hs.fileout','off').lower() != 'on':
            raise ValueError('HS.fileout on required for native artifact acceptance')
        if not absolute(v.get('data.path',''), 'DATA.PATH').is_dir():
            raise ValueError('DATA.PATH must be a directory')
    if n.get('scf.spinpolarization','off').lower()!='off' or n.get('scf.spinorbit.coupling','off').lower()!='off':
        raise ValueError('NSOC input must be nonmagnetic scalar')
    if n.get('scf.restart','off').lower() not in ('off','no'):
        raise ValueError('NSOC must start independently with scf.restart off')
    if int(n.get('scf.maxiter','200')) <= 1:
        raise ValueError('NSOC must allow complete convergence, not one-step smoke')
    for key, expected in {'scf.spinpolarization':'nc','scf.spinorbit.coupling':'on','scf.restart':'c2n','scf.maxiter':'1'}.items():
        if s.get(key,'').lower()!=expected:
            raise ValueError('one-step SOC requires '+key+' '+expected)
    if s.get('scf.restart.filename') != n['system.name']:
        raise ValueError('SOC restart must name this workflow own NSOC checkpoint')
    # Reject implicit/unequal important physical defaults instead of guessing them.
    for key in ('data.path','scf.xctype','scf.energycutoff','scf.electronictemperature'):
        if key not in n or key not in s:
            raise ValueError('explicit physical setting required in both inputs: '+key)
        different = str(Path(n[key]).resolve()) != str(Path(s[key]).resolve()) if key=='data.path' else _numeric(n[key]) != _numeric(s[key])
        if different:
            raise ValueError('NSOC/SOC physical setting differs: '+key)
    if float(n.get('scf.system.charge','0')) != float(s.get('scf.system.charge','0')):
        raise ValueError('NSOC/SOC effective scf.system.charge differs')
    # The SOC criterion is irrelevant to a one-step solve, but never rewrite it.
    return n,s


def clean_environment(inherited, additions):
    blocked = set(THREAD_KEYS) | {'MKL_DOMAIN_NUM_THREADS','OMP_DYNAMIC','MKL_DYNAMIC','OMP_PROC_BIND','OMP_PLACES','OMP_MAX_ACTIVE_LEVELS','PYTHONPATH','PYTHONHOME','GOMP_CPU_AFFINITY'}
    result = {k:v for k,v in inherited.items() if not k.startswith(('OPENMX_CPU_','I_MPI_','FI_','UCX_','OMP_','MKL_','OPENBLAS_','NUMEXPR_','BLIS_','KMP_')) and k not in blocked}
    result.update(additions)
    return result


def enabled(env,key):
    return str(env.get(key,'0')).strip().lower() not in ('0','','off','false','no')


def validate_profile(stage, profile, symmetry):
    if not isinstance(profile,dict) or not profile.get('name') or not isinstance(profile.get('environment'),dict):
        raise ValueError('explicit named '+stage+' profile/environment required')
    env = dict(profile['environment'])
    for key, value in env.items():
        if not re.fullmatch(r'[A-Z][A-Z0-9_]*',key) or not isinstance(value,str) or '\x00' in value:
            raise ValueError('environment must contain uppercase names and string values')
        if not key.startswith(('OPENMX_CPU_','OMP_','MKL_','OPENBLAS_','NUMEXPR_')):
            raise ValueError('runtime transport/binding is controlled by runtime/layout, not profile: '+key)
    if stage=='soc':
        generic = enabled(env,'OPENMX_CPU_GSOC')
        if generic and any(enabled(env,k) for k in ('OPENMX_CPU_C3','OPENMX_CPU_D3','OPENMX_CPU_SOC_COMPACT','OPENMX_CPU_SOC_TR','OPENMX_CPU_SOC_REAL')):
            raise ValueError('conflicting generic and legacy SOC flags')
        project = bool(symmetry.get('project',False)); compact = bool(symmetry.get('compact',False))
        if enabled(env,'OPENMX_CPU_GSOC_PROJECT') != project:
            raise ValueError('SOC project flag must agree with explicit symmetry.project')
        if enabled(env,'OPENMX_CPU_GSOC_COMPACT') != compact:
            raise ValueError('SOC compact flag must agree with explicit symmetry.compact')
        if (project or compact) and not generic:
            raise ValueError('generic SOC must be enabled for project/compact')
        if compact and (env.get('OPENMX_CPU_DM_BLAS3')!='1' or any(enabled(env,k) for k in ('OPENMX_CPU_DM_SORT','OPENMX_CPU_DM_TRIM','OPENMX_CPU_DM_CYCLIC'))):
            raise ValueError('compact SOC requires DM_BLAS3=1 and old SORT/TRIM/CYCLIC disabled')
    if stage=='soc':
        policy=symmetry.get('soc_export_policy','raw')
        if policy not in ('raw','symmetrize'):
            raise ValueError('soc_export_policy must be raw or symmetrize')
        generic=enabled(env,'OPENMX_CPU_GSOC')
        native_group=env.get('OPENMX_CPU_SOC_NATIVE_GROUP')
        native_selector=env.get('OPENMX_CPU_SOC_NATIVE_D3')
        if native_selector is not None and native_selector not in ('0','1'):
            raise ValueError('legacy SOC_NATIVE_D3 selector must be 0 or 1')
        if native_group is not None:
            if native_group not in ('C3','D3'):
                raise ValueError('legacy native group must be C3 or D3')
            if native_selector is not None and (native_selector=='0' or native_group!='D3'):
                raise ValueError('legacy native group and SOC_NATIVE_D3 selectors conflict')
        native_active=native_group is not None or native_selector=='1'
        if generic and policy!='raw':
            raise ValueError('generic GSOC currently supports raw H/S export only')
        if native_active and policy!='symmetrize':
            raise ValueError('legacy native export projects H/S: explicit soc_export_policy=symmetrize required')
        if policy=='symmetrize' and not native_active:
            raise ValueError('symmetrized export requires explicit legacy native-map activation')
        if native_active:
            group=native_group or 'D3'
            if enabled(env,'OPENMX_CPU_D3')!=(group=='D3'):
                raise ValueError('legacy native group and D3 solver selectors disagree')
            if env.get('OPENMX_CPU_C3')!='3' or not env.get('OPENMX_CPU_SOC_NATIVE_MAP'):
                raise ValueError('legacy native export requires C3 mode3 and explicit native map')
        if any(enabled(env,k) for k in ('OPENMX_CPU_SOC_TR','OPENMX_CPU_SOC_REAL')):
            if not native_active or (native_group or 'D3')!='D3' or not enabled(env,'OPENMX_CPU_D3') or not env.get('OPENMX_CPU_SOC_NATIVE_MAP'):
                raise ValueError('legacy SOC TR/REAL requires D3 with activated native map; C3 must keep TR/REAL off')
    if stage=='nsoc' and enabled(env,'OPENMX_CPU_NSOC_SYM') and not env.get('OPENMX_CPU_NSOC_BASIS'):
        raise ValueError('NSOC symmetry requires explicit OPENMX_CPU_NSOC_BASIS')
    if stage=='soc' and enabled(env,'OPENMX_CPU_GSOC') and not env.get('OPENMX_CPU_GSOC_BASIS'):
        raise ValueError('generic SOC requires explicit OPENMX_CPU_GSOC_BASIS')
    if stage=='soc' and any(enabled(env,k) for k in ('OPENMX_CPU_C3','OPENMX_CPU_SOC_COMPACT')) and not env.get('OPENMX_CPU_C3_BASIS'):
        raise ValueError('legacy SOC requires explicit OPENMX_CPU_C3_BASIS')
    if stage=='soc' and any(enabled(env,k) for k in ('OPENMX_CPU_NSOC_SYM','OPENMX_CPU_NSOC_FULL_PROJECT')):
        raise ValueError('conflicting scalar NSOC flags in SOC profile')
    return env


def plan_config(config):
    if config.get('schema') != 1:
        raise ValueError('workflow schema must be 1')
    binary = absolute(config['binary'],'binary',file=True)
    if not os.access(binary,os.X_OK):
        raise ValueError('binary is not executable')
    sources = {stage:absolute(config['inputs'][stage],stage+' input',file=True) for stage in ('nsoc','soc')}
    n,s = validate_inputs(sources['nsoc'].read_text(),sources['soc'].read_text())
    output = absolute(config['output'],'output',exists=False)
    for source in sources.values():
        if output.is_relative_to(source.parent) or source.is_relative_to(output):
            raise ValueError('output must be outside input source directories, including symlink aliases')
    if output.exists():
        raise ValueError('output already exists; use a fresh run directory')
    layout = config['layout']; hosts = layout['hosts']
    if not isinstance(hosts,list) or len(hosts) not in (1,2,4) or len(set(hosts))!=len(hosts) or any(not isinstance(h,str) or not HOST_NAME.fullmatch(h) for h in hosts):
        raise ValueError('hosts must be 1, 2 or 4 distinct valid SSH host names')
    per_host = positive_int(layout['ranks_per_host'],'ranks_per_host'); threads = positive_int(layout['threads'],'threads')
    ids = layout.get('cpu_ids_per_host',{})
    if set(ids)!=set(hosts):
        raise ValueError('explicit CPU ownership cpu_ids_per_host required for every host')
    for h in hosts:
        cpus = ids[h]
        if not isinstance(cpus,list) or len(cpus)!=per_host*threads or any(isinstance(x,bool) or not isinstance(x,int) or x<0 for x in cpus) or len(set(cpus))!=len(cpus):
            raise ValueError('CPU ids must be distinct and have ranks_per_host * threads entries')
    runtime = config['runtime']
    mpiexec = absolute(runtime['mpiexec'],'mpiexec',file=True)
    prefix = absolute(runtime['runtime_prefix'],'runtime_prefix')
    python_prefix = absolute(runtime['environment_prefix'],'environment_prefix')
    python = python_prefix/'bin/python'
    if not python.is_file():
        raise ValueError('environment_prefix/bin/python missing')
    transport = runtime.get('transport','shm:ofi' if len(hosts)>1 else 'shm')
    if transport not in ('shm','shm:ofi','ofi') or len(hosts)>1 and transport=='shm':
        raise ValueError('multi-node requires explicit OFI transport')
    symmetry = config.get('symmetry',{})
    for flag in ('project','compact'):
        if flag in symmetry and not isinstance(symmetry[flag],bool):
            raise ValueError('symmetry.'+flag+' must be boolean')
    for name in ('geometry_tolerance_ang','covariance_tolerance'):
        if name in symmetry:
            value=symmetry[name]
            if isinstance(value,bool) or not isinstance(value,(float,int)) or not math.isfinite(value) or value<=0:
                raise ValueError('symmetry.'+name+' must be a finite positive number')
    soc_env=config['profiles']['soc'].get('environment',{})
    if enabled(soc_env,'OPENMX_CPU_GSOC'):
        for field,key in (('geometry_tolerance_ang','OPENMX_CPU_GSOC_GEOMETRY_TOL_ANG'),('covariance_tolerance','OPENMX_CPU_GSOC_COV_TOL')):
            if field in symmetry and (key not in soc_env or float(soc_env[key])!=symmetry[field]):
                raise ValueError('explicit symmetry.'+field+' must match SOC environment '+key)
    profiles = {stage:{'name':config['profiles'][stage]['name'],'environment':validate_profile(stage,config['profiles'][stage],symmetry)} for stage in ('nsoc','soc')}
    if threads>1:
        for stage,profile in profiles.items():
            if enabled(profile['environment'],'OPENMX_CPU_VNA_SHARED'):
                raise ValueError(stage+' OPENMX_CPU_VNA_SHARED requires exactly one thread; explicitly disable it for threaded layouts')
            if enabled(profile['environment'],'OPENMX_CPU_GRID_NUMA'):
                raise ValueError(stage+' OPENMX_CPU_GRID_NUMA must be0 in this driver mixed-thread profile until full-team NUMA locality is audited')
    scalar_env=profiles['nsoc']['environment']
    realspace_mode=scalar_env.get('OPENMX_CPU_REALSPACE_D3','off')
    if realspace_mode not in ('0','off','export','scf','1'):
        raise ValueError('OPENMX_CPU_REALSPACE_D3 must be off,export,scf (or native aliases0/1)')
    if enabled(scalar_env,'OPENMX_CPU_NSOC_SECTOR_PARALLEL'):
        world=len(hosts)*per_host
        if len(hosts)!=4 or not ((world==256 and threads==1) or (world==64 and threads==4)):
            raise ValueError('NSOC sector parallel layout requires4hosts with256MPIx1 or64MPIx4')
        for host in hosts:
            if set(ids[host])!=set(range(64)):
                raise ValueError('NSOC sector parallel CPU IDs must cover exactly0..63 on every host')
    if enabled(scalar_env,'OPENMX_CPU_NSOC_SECTOR_PARALLEL'):
        parts=scalar_env.get('OPENMX_CPU_NSOC_SECTOR_RANKS','').split(',')
        if len(parts)!=3 or any(not x.strip().isdigit() or int(x)<=0 for x in parts) or sum(map(int,parts))!=len(hosts)*per_host:
            raise ValueError('explicit NSOC sector rank partition requires three positive entries summing to current MPI rank count')
    elpa_threads=scalar_env.get('OPENMX_CPU_NSOC_ELPA_OMP_THREADS','0')
    if elpa_threads not in ('0','1'):
        raise ValueError('NSOC ELPA_OMP_THREADS must be0 or1')
    if elpa_threads=='1' and (not enabled(scalar_env,'OPENMX_CPU_NSOC_SECTOR_PARALLEL') or len(hosts)*per_host!=64 or threads!=4):
        raise ValueError('NSOC ELPA_OMP_THREADS=1 requires enabled sector parallel with64MPIx4')
    dependencies = {str(binary):sha256(binary), str(mpiexec):sha256(mpiexec), **{str(p):sha256(p) for p in sources.values()}}
    # Every explicitly referenced native basis/map is frozen; no relative metadata.
    for stage in profiles:
        for key,value in profiles[stage]['environment'].items():
            required=key.endswith(('_BASIS','_MAP')) and value not in ('','0','none')
            other_input_file=Path(value).is_absolute() and Path(value).is_file()
            if required or other_input_file:
                path = absolute(value,key,file=True)
                # Freeze the path consumed by native code, not only the hash
                # target: retargeting the original symlink cannot redirect it.
                profiles[stage]['environment'][key]=str(path)
                dependencies[str(path)] = sha256(path)
    manifest = None
    if config.get('build_manifest'):
        manifest_path = absolute(config['build_manifest'],'build_manifest',file=True)
        manifest = load(manifest_path)
        if manifest.get('binary',{}).get('sha256')!=sha256(binary):
            raise ValueError('build manifest binary hash differs')
        if manifest.get('external_openmx_objects') or manifest.get('external_elpa_objects'):
            raise ValueError('standalone build cannot require external OpenMX/ELPA objects')
        caps = manifest.get('capabilities',{})
        for stage,env in ((x,profiles[x]['environment']) for x in profiles):
            required = []
            if enabled(env,'OPENMX_CPU_NSOC_SYM'): required.append('nsoc_symmetry')
            if enabled(env,'OPENMX_CPU_GSOC'): required.append('generic_soc')
            if stage=='soc' and any(enabled(env,k) for k in ('OPENMX_CPU_C3','OPENMX_CPU_D3','OPENMX_CPU_SOC_COMPACT')): required.append('legacy_soc')
            if enabled(env,'OPENMX_CPU_GSOC_COMPACT'): required.append('compact_soc')
            if len(hosts)>1: required.append('multinode')
            for cap in required:
                if not caps.get(cap):
                    raise ValueError('source build lacks capability '+cap)
        dependencies[str(manifest_path)] = sha256(manifest_path)
    for tool in (Path(__file__).resolve(), HERE/'workflow_node.py'):
        dependencies[str(tool)] = sha256(tool)
    interval = float(config.get('monitor_interval_seconds',0.25))
    if not 0.05<=interval<=5:
        raise ValueError('monitor interval must be in [0.05,5] seconds')
    return {'schema':1,'category':'full converged NSOC then own exactly-one-step SOC; native H/S and original physics retained',
        'binary':str(binary),'binary_sha256':sha256(binary),'build_manifest':str(Path(config['build_manifest']).resolve()) if config.get('build_manifest') else None,
        'inputs':{k:str(v) for k,v in sources.items()},'input_hashes':{k:sha256(v) for k,v in sources.items()},
        'output':str(output),'scratch_root':str(Path('/data/work/tmp')/os.environ.get('USER','zy')/'codex/openmx-unified'/hashlib.sha256(str(output).encode()).hexdigest()[:20]),'layout':layout,'hosts':hosts,'ranks_per_host':per_host,'ranks':len(hosts)*per_host,'threads':threads,
        'total_cores':len(hosts)*per_host*threads,'runtime':{**runtime,'mpiexec':str(mpiexec),'runtime_prefix':str(prefix),'environment_prefix':str(python_prefix),'transport':transport},
        'profiles':profiles,'symmetry':{**symmetry,'soc_export_policy':symmetry.get('soc_export_policy','raw')},'soc_export_policy':symmetry.get('soc_export_policy','raw'),'dependencies_sha256':dependencies,'system_names':{'nsoc':n['system.name'],'soc':s['system.name']},
        'monitor_interval_seconds':interval,'memory_limit_bytes_per_host':positive_int(config['memory_limit_bytes_per_host'],'memory_limit_bytes_per_host') if config.get('memory_limit_bytes_per_host') is not None else None,'timeout_seconds':positive_int(config.get('timeout_seconds',86400),'timeout_seconds'),
        'environment_python':str(python),'physics_input_modified':False,'input_parser_auto_enabled':False,
        'nsoc_density_projection_explicit':False,
        'deprecated_fields':{'nsoc_density_projection_explicit':'Deprecated misnomer: REALSPACE_D3 projects H/S; AfterDensity checks DM/EDM without writing them. This field is always false.'},
        'nsoc_realspace_HS_projection_explicit':realspace_mode in ('export','scf','1'),
        'nsoc_realspace_HS_projection_scope':'each_scf_and_export' if realspace_mode in ('scf','1') else 'export_only' if realspace_mode=='export' else 'none',
        'nsoc_export_HS':'PROJECTED_REALSPACE_H_S' if realspace_mode in ('export','scf','1') else 'RAW_REALSPACE_H_S',
        'nsoc_sector_solver_all_steps':enabled(profiles['nsoc']['environment'],'OPENMX_CPU_NSOC_FULL_PROJECT'),
        'source_commit':subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip() if (REPO/'.git').exists() else None}


def parse_native(stdout, native_out):
    def last(pattern,text,conv=float):
        hits = re.findall(pattern,text,re.MULTILINE)
        return conv(hits[-1]) if hits else None
    return {'normal_finish':'The calculation was normally finished.' in stdout,
            'native_converged':bool(re.search(r'The SCF was achieved at MD=\s*1',stdout)),
            'SCF_steps':last(r'MD=\s*1\s+SCF=\s*(\d+)',stdout,int),
            'energy_Ha':last(r'^\s*Utot\.\s+([-+\d.eE]+)',native_out),
            'native_wall_seconds':last(r'^\s*Elapsed\.Time\.\s+([-+\d.eE]+)',native_out),
            'restart_read_success':bool(re.search(r'<Restart>\s+Found restart files',stdout))}


def check_checkpoint(directory, prefix, ranks):
    directory = Path(directory)
    if not directory.is_dir() or directory.is_symlink():
        raise ValueError('own NSOC checkpoint missing or symlinked')
    check = directory/(prefix+'.crst_check')
    if not check.is_file() or int(check.read_text().split()[0])!=ranks:
        raise ValueError('checkpoint rank count does not match SOC layout')
    required = [check] + [directory/(prefix+'.adrst'+str(i)) for i in range(ranks)] + [directory/(prefix+'.crst0_'+str(i)+'_0') for i in range(ranks)]
    if any(not p.is_file() or p.is_symlink() or p.stat().st_size==0 for p in required):
        raise ValueError('missing or invalid own NSOC restart checkpoint files')
    if any(p.is_symlink() for p in directory.rglob('*')):
        raise ValueError('checkpoint tree contains symlinks')
    return directory


def ssh(host, command, *, check=True, timeout=35):
    return subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8','-o','ServerAliveInterval=5','-o','ServerAliveCountMax=3',host,shlex.join([str(x) for x in command])],stdin=subprocess.DEVNULL,capture_output=True,text=True,check=check,timeout=timeout)


def _library_environment(plan,stage):
    env = dict(plan['profiles'][stage]['environment'])
    if plan.get('scratch_root'):
        scratch=Path(plan['scratch_root'])/stage
        scratch.mkdir(parents=True,exist_ok=False)
        env['TMPDIR']=str(scratch)
    env.update({k:str(plan['threads']) for k in THREAD_KEYS})
    env.update(OMP_DYNAMIC='FALSE',MKL_DYNAMIC='FALSE',OMP_MAX_ACTIVE_LEVELS='1',I_MPI_PIN='0',
               I_MPI_FABRICS=plan['runtime']['transport'],I_MPI_ROOT=str(Path(plan['runtime']['mpiexec']).parent.parent),I_MPI_DEBUG='5',I_MPI_HYDRA_BOOTSTRAP='ssh',
               OPENMX_WORKFLOW_STAGE=stage)
    if plan['runtime'].get('provider'):
        env['FI_PROVIDER']=plan['runtime']['provider']
    prefix = Path(plan['runtime']['runtime_prefix'])
    dirs = [prefix/'mkl/2022.2.1/lib/intel64',prefix/'compiler/2022.2.1/linux/compiler/lib/intel64_lin',prefix/'mpi/2021.7.1/lib/release',prefix/'mpi/2021.7.1/lib']
    extra = plan['runtime'].get('library_paths',[])
    if any(not Path(p).is_absolute() for p in extra):
        raise ValueError('runtime.library_paths must be absolute')
    # Do not retain unrelated inherited MPI/BLAS libraries.
    env['LD_LIBRARY_PATH']=':'.join([str(p) for p in dirs if p.is_dir()]+extra)
    env['PATH']=str(Path(plan['runtime']['mpiexec']).parent)+':/usr/bin:/bin:/usr/sbin:/sbin'
    return clean_environment(os.environ,env)


def _units(plan,stage):
    token = hashlib.sha256(plan['output'].encode()).hexdigest()[:16]+stage
    return 'omxwf'+token+'.slice','omxwf'+token+'monitor.service'


def start_monitors(plan, stage, directory):
    slice_name, unit = _units(plan,stage)
    units = {}
    logs=directory/'monitor_logs'; logs.mkdir()
    for host in plan['hosts']:
        command = ['systemd-run','--user','--unit='+unit,'--slice='+slice_name,'--property=MemoryAccounting=yes',
                   '--property=RuntimeMaxSec='+str(plan['timeout_seconds']+120),'--property=TasksMax=infinity','--property=StandardOutput=file:'+str(logs/(host+'.stdout.log')),
                   '--property=StandardError=file:'+str(logs/(host+'.stderr.log')),
                   plan['environment_python'],str(HERE/'workflow_node.py'),'monitor','--directory',str(directory),
                   '--host',host,'--interval',str(plan['monitor_interval_seconds'])]
        q = ssh(host,command)
        # Register immediately, before any follow-up can fail, so partial
        # startup is still cleaned up by the stage finally block.
        units[host]={'unit':unit,'slice':slice_name,'command':command,'submission_stdout':q.stdout,'submission_stderr':q.stderr}
        write_json(directory/'monitor_units.json',units)
        if plan.get('memory_limit_bytes_per_host'):
            ssh(host,['systemctl','--user','set-property','--runtime',slice_name,'MemoryMax='+str(plan['memory_limit_bytes_per_host']),'MemorySwapMax=0'])
    deadline = time.monotonic()+30
    while time.monotonic()<deadline:
        if all((directory/'nodes'/host/'ready.json').exists() for host in plan['hosts']):
            return units
        time.sleep(.1)
    raise RuntimeError('per-host physical-memory monitor did not become ready')


def stop_monitors(plan,stage,directory,units):
    (directory/'MONITOR_STOP').touch()
    deadline = time.monotonic()+15
    while time.monotonic()<deadline and not all((directory/'nodes'/h/'final.json').exists() for h in units):
        time.sleep(.1)
    for host,u in units.items():
        q=ssh(host,['systemctl','--user','show',u['unit'],'-p','ActiveState','-p','SubState','-p','Result','-p','ExecMainStatus'],check=False)
        u['final_state']=q.stdout; u['show_exit_code']=q.returncode
        # Native scopes already exited (or failed-stage cleanup stopped them).
        # Remove this task's empty slice after saving kernel peak accounting.
        stopped=ssh(host,['systemctl','--user','stop',u['slice']],check=False)
        u['slice_stop_exit_code']=stopped.returncode
        final=ssh(host,['systemctl','--user','show',u['slice'],'-p','ActiveState','-p','SubState'],check=False)
        u['final_slice_state']=final.stdout
    write_json(directory/'monitor_units.json',units)


def memory_summary(directory,hosts,interval):
    records={}; peaks={}
    for host in hosts:
        d=directory/'nodes'/host
        final=load(d/'final.json')
        peaks[host]=final['slice_peak_bytes']
        records[host]=[json.loads(line) for line in (d/'samples.jsonl').read_text().splitlines() if line.strip()]
    # A sample bundle is simultaneous within one common wall-clock interval. Summing each
    # host's independent maximum is a separate conservative upper bound.
    indexed={h:{r['sample_id']:r for r in records[h]} for h in hosts}
    common=set.intersection(*(set(indexed[h]) for h in hosts))
    sums=[]; skews=[]
    for sample_id in sorted(common):
        selected=[indexed[h][sample_id] for h in hosts]
        skew=max(x['time_unix'] for x in selected)-min(x['time_unix'] for x in selected)
        if skew<=interval:
            sums.append(sum(x['slice_current_bytes'] for x in selected)); skews.append(skew)
    return {'metric':'per-stage per-host cgroup v2 slice memory.current / memory.peak, shared pages counted once per host; includes monitor, excludes remote MPI bootstrap outside slice',
            'per_node_peak_bytes':peaks,'max_node_peak_bytes':max(peaks.values()),
            'synchronized_cluster_peak_bytes':max(sums) if sums else None,
            'sum_independent_node_peaks_bytes':sum(peaks.values()),'synchronized_samples':len(sums),
            'sample_ids_intersection':len(common),
            'maximum_sample_skew_seconds':max(skews) if skews else None,'sample_interval_seconds':interval}


def _verify_rank_records(plan,directory):
    files=list((directory/'ranks').glob('rank*.json'))
    if len(files)!=plan['ranks']:
        raise ValueError('actual MPI rank witness count differs from requested layout')
    seen=set()
    for file in files:
        r=load(file); rank=r['rank']; host=r['host']; local=r['local_rank']
        if rank in seen or rank<0 or rank>=plan['ranks'] or host not in plan['hosts']:
            raise ValueError('invalid rank/host witness')
        seen.add(rank)
        expected=plan['layout']['cpu_ids_per_host'][host][local*plan['threads']:(local+1)*plan['threads']]
        if sorted(r['cpu_affinity'])!=sorted(expected):
            raise ValueError('rank actual binding differs from explicit CPU list')
        if not r['inside_task_slice']:
            raise ValueError('MPI rank escaped monitored task slice')
    return [str(p) for p in sorted(files)]


def verify_dependencies(plan):
    for path,want in plan['dependencies_sha256'].items():
        if sha256(path)!=want:
            raise ValueError('frozen planned dependency changed: '+path)


def start_anchors(plan):
    directory=Path(plan['output'])/'anchors'; directory.mkdir()
    anchors={}
    try:
        for host in plan['hosts']:
            command=[plan['environment_python'],str(HERE/'workflow_node.py'),'anchor','--directory',str(directory),'--host',host]
            log=(directory/(host+'.stderr.log')).open('w')
            q=subprocess.Popen(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8','-o','ServerAliveInterval=5','-o','ServerAliveCountMax=3',host,shlex.join(command)],
                               stdin=subprocess.PIPE,stdout=subprocess.DEVNULL,stderr=log)
            anchors[host]=(q,log)
        deadline=time.monotonic()+20
        while time.monotonic()<deadline:
            if all((directory/(host+'.ready.json')).exists() for host in anchors): return anchors
            if any(q.poll() is not None for q,log in anchors.values()):
                raise RuntimeError('persistent SSH login anchor failed; inspect anchors/*.stderr.log')
            time.sleep(.1)
        raise RuntimeError('persistent SSH login anchors did not become ready')
    except BaseException:
        stop_anchors(plan,anchors)
        raise


def stop_anchors(plan,anchors):
    directory=Path(plan['output'])/'anchors'
    (directory/'STOP').touch()
    for host,(q,log) in anchors.items():
        if q.stdin: q.stdin.close()
        try: q.wait(timeout=15)
        except subprocess.TimeoutExpired:
            q.terminate()
            try:q.wait(timeout=5)
            except subprocess.TimeoutExpired:q.kill();q.wait()
        log.close()


def run_stage(plan,stage):
    verify_dependencies(plan)
    directory=Path(plan['output'])/stage
    directory.mkdir(); (directory/'ranks').mkdir(); (directory/'nodes').mkdir()
    shutil.copy2(plan['inputs'][stage],directory/'input.dat')
    prefix=plan['system_names'][stage]
    parent_checkpoint=None
    if stage=='soc':
        parent=Path(plan['output'])/'nsoc'; parent_meta=load(parent/'run.json')
        if parent_meta['state']!='complete' or not parent_meta['native_converged'] or parent_meta['ranks']!=plan['ranks'] or parent_meta['threads']!=plan['threads']:
            raise ValueError('SOC requires own complete converged NSOC with same MPI/thread layout')
        name=plan['system_names']['nsoc']
        rst=check_checkpoint(parent/(name+'_rst'),name,plan['ranks'])
        files={str(f.relative_to(rst)):sha256(f) for f in rst.rglob('*') if f.is_file()}
        shutil.copytree(rst,directory/rst.name)
        for name,want in files.items():
            if sha256(rst/name)!=want or sha256(directory/rst.name/name)!=want:
                raise ValueError('NSOC checkpoint changed during staging: '+name)
        verify_dependencies(plan)
        parent_checkpoint={'source':str(rst),'own_nsoc':str(parent),'ranks':plan['ranks'],'threads':plan['threads'],'files_sha256':files}
        write_json(directory/'RESTART_PROVENANCE.json',parent_checkpoint)
    meta={'state':'starting','stage':stage,'ranks':plan['ranks'],'threads':plan['threads'],'hosts':plan['hosts'],
          'binary':plan['binary'],'binary_sha256':plan['binary_sha256'],'profile':plan['profiles'][stage],
          'input_sha256':sha256(directory/'input.dat'),'physics_input_modified':False,
          'soc_export_policy':plan['soc_export_policy'] if stage=='soc' else None,
          'nsoc_density_projection_explicit':False if stage=='nsoc' else None,
          'deprecated_fields':plan['deprecated_fields'] if stage=='nsoc' else {},
          'nsoc_realspace_HS_projection_explicit':plan['nsoc_realspace_HS_projection_explicit'] if stage=='nsoc' else None,
          'nsoc_realspace_HS_projection_scope':plan['nsoc_realspace_HS_projection_scope'] if stage=='nsoc' else None,
          'nsoc_export_HS':plan['nsoc_export_HS'] if stage=='nsoc' else None,
          'nsoc_sector_solver_all_steps':plan['nsoc_sector_solver_all_steps'] if stage=='nsoc' else None,'own_checkpoint':parent_checkpoint,
          'environment_prefix':plan['runtime']['environment_prefix'],'started_unix':time.time()}
    write_json(directory/'run.json',meta)
    units={}; q=None
    try:
        units=start_monitors(plan,stage,directory)
        env=_library_environment(plan,stage)
        write_json(directory/'environment.json',{k:v for k,v in env.items() if k.startswith(('OPENMX_','OMP_','MKL_','I_MPI_','FI_')) or k in THREAD_KEYS or k=='LD_LIBRARY_PATH'})
        cmd=[plan['runtime']['mpiexec'],'-hosts',','.join(plan['hosts']),'-ppn',str(plan['ranks_per_host']),'-np',str(plan['ranks']),
             plan['environment_python'],str(HERE/'workflow_node.py'),'rank','--plan',str(Path(plan['output'])/'plan.json'),'--stage',stage]
        meta.update(state='running',command=cmd); write_json(directory/'run.json',meta)
        verify_dependencies(plan)
        if sha256(directory/'input.dat')!=plan['input_hashes'][stage]:
            raise ValueError('staged input changed before MPI launch')
        start=time.monotonic()
        with (directory/'stdout.log').open('wb') as out,(directory/'stderr.log').open('wb') as err:
            q=subprocess.run(cmd,cwd=directory,env=env,stdin=subprocess.DEVNULL,stdout=out,stderr=err,timeout=plan['timeout_seconds'])
        meta['launch_wall_seconds']=time.monotonic()-start
        combined_log=(directory/'stdout.log').read_text(errors='replace')+'\n'+(directory/'stderr.log').read_text(errors='replace')
        providers=re.findall(r'libfabric provider:\s*(\S+)',combined_log)
        versions=re.findall(r'Intel\(R\) MPI Library, Version ([^\n]+)',combined_log)
        meta['runtime_observed']={'MPI_versions':sorted(set(versions)),'libfabric_providers':sorted(set(providers)),
                                  'library_resolution':str(directory/'LOADED_LIBRARY_RESOLUTION.txt')}
        out=(directory/(prefix+'.out'))
        metrics=parse_native((directory/'stdout.log').read_text(errors='replace'),out.read_text(errors='replace') if out.exists() else '')
        hs=directory/(prefix+'.scfout')
        success=q.returncode==0 and metrics['normal_finish'] and metrics['energy_Ha'] is not None and metrics['native_wall_seconds'] is not None and hs.is_file() and hs.stat().st_size>32
        success=success and (metrics['native_converged'] if stage=='nsoc' else metrics['restart_read_success'] and metrics['SCF_steps']==1)
        if stage=='nsoc' and success:
            check_checkpoint(directory/(prefix+'_rst'),prefix,plan['ranks'])
        rank_records=_verify_rank_records(plan,directory)
        meta.update(exit_code=q.returncode,normal_finish=metrics['normal_finish'],native_converged=metrics['native_converged'],
                    metrics=metrics,wall_seconds=metrics['native_wall_seconds'],HS_output=str(hs) if hs.exists() else None,
                    rank_records=rank_records,artifacts_inspected=bool(success),state='complete' if success else 'failed')
    except BaseException as exc:
        meta.update(state='failed',error=str(exc),traceback=traceback.format_exc())
    finally:
        if not units and (directory/'monitor_units.json').exists():
            units=load(directory/'monitor_units.json')
        if units:
            try:
                # Stop only this task's rank scopes if a timeout/native failure left
                # children. This never scans or targets unrelated OpenMX processes.
                if meta['state']=='failed':
                    slice_name,_=_units(plan,stage)
                    token=slice_name[:-6]
                    for hi,host in enumerate(plan['hosts']):
                        rank_units=[token+'rank'+str(rank)+'.scope' for rank in range(hi*plan['ranks_per_host'],(hi+1)*plan['ranks_per_host'])]
                        ssh(host,['systemctl','--user','stop',*rank_units],check=False)
                stop_monitors(plan,stage,directory,units)
                meta['memory']=memory_summary(directory,plan['hosts'],plan['monitor_interval_seconds'])
                if not meta['memory']['synchronized_samples']:
                    raise ValueError('no simultaneous physical-memory samples')
                for host in plan['hosts']:
                    node_final=load(directory/'nodes'/host/'final.json')
                    if node_final['memory_events'].get('oom_kill',0):
                        raise ValueError('task cgroup recorded OOM kill on '+host)
                    if node_final.get('remaining_descendant_pids'):
                        raise ValueError('task cgroup retains native processes on '+host)
            except BaseException as exc:
                meta.update(state='failed',monitor_error=str(exc))
        else:
            meta.update(state='failed',monitor_error='no task memory monitors started')
        meta['ended_unix']=time.time(); write_json(directory/'run.json',meta)
    return meta


def controller(plan_path):
    plan_path=Path(plan_path); plan=load(plan_path); output=Path(plan['output'])
    if socket.gethostname().split('.')[0]!=plan['hosts'][0].split('.')[0]:
        raise ValueError('controller must run on first configured authorized compute host')
    verify_dependencies(plan)
    lock=output/'CONTROLLER_STARTED'
    with lock.open('x') as f: f.write(str(os.getpid())+'\n')
    def terminate(signum,frame):
        raise RuntimeError('controller received termination signal '+str(signum))
    signal.signal(signal.SIGTERM,terminate)
    signal.signal(signal.SIGINT,terminate)
    state={'schema':1,'state':'running','hosts':plan['hosts'],'ranks':plan['ranks'],'threads':plan['threads'],
           'total_cores':plan['total_cores'],'binary_sha256':plan['binary_sha256'],'stages':{},'started_unix':time.time()}
    write_json(output/'workflow.json',state)
    anchors={}
    try:
        anchors=start_anchors(plan)
        for stage in ('nsoc','soc'):
            result=run_stage(plan,stage); state['stages'][stage]=result; write_json(output/'workflow.json',state)
            if result['state']!='complete':
                raise RuntimeError(stage+' did not pass native/lifecycle acceptance')
        stages=list(state['stages'].values())
        state.update(state='complete',full_native_wall_seconds=sum(s['wall_seconds'] for s in stages),
                     full_launch_wall_seconds=sum(s['launch_wall_seconds'] for s in stages),
                     core_hours=sum(s['wall_seconds'] for s in stages)*plan['total_cores']/3600,
                     synchronized_cluster_peak_bytes=max(s['memory']['synchronized_cluster_peak_bytes'] for s in stages),
                     max_node_peak_bytes=max(s['memory']['max_node_peak_bytes'] for s in stages),
                     science_comparison='not inferred from normal exit; compare matching H/S and density diagnostics to reference separately')
    except BaseException as exc:
        state.update(state='failed',error=str(exc),traceback=traceback.format_exc())
    finally:
        if anchors: stop_anchors(plan,anchors)
        state['ended_unix']=time.time(); write_json(output/'workflow.json',state)
    return 0 if state['state']=='complete' else 1


def submit(config, *, wait=False):
    if not wait:
        raise ValueError('run requires --wait: foreground SSH keeps non-lingering user managers alive; disconnect is not durable')
    plan=plan_config(config)
    if not plan['build_manifest']:
        raise ValueError('run requires source-built build_manifest; plan may inspect legacy candidates read-only')
    directory=Path(plan['output']); directory.mkdir(parents=True,exist_ok=False)
    write_json(directory/'config.json',config); write_json(directory/'plan.json',plan)
    unit='omxwfcontroller'+hashlib.sha256(str(directory).encode()).hexdigest()[:16]+'.service'
    command=['systemd-run','--user','--wait','--unit='+unit,'--property=MemoryAccounting=yes','--property=TimeoutStopSec=120','--property=RuntimeMaxSec='+str(2*plan['timeout_seconds']+600),
             '--property=StandardOutput=file:'+str(directory/'controller.stdout.log'),'--property=StandardError=file:'+str(directory/'controller.stderr.log'),
             '--working-directory='+str(directory),plan['environment_python'],str(Path(__file__).resolve()),'_controller',str(directory/'plan.json')]
    submission={'job_id':None,'scheduler_submission':False,'authorization':'user explicit direct SSH to named hosts','host':plan['hosts'][0],
                'unit':unit,'command':command,'environment_prefix':plan['runtime']['environment_prefix'],'ranks':plan['ranks'],'threads':plan['threads'],
                'logs':[str(directory/'controller.stdout.log'),str(directory/'controller.stderr.log')],'result_path':str(directory),'submitted_unix':time.time()}
    write_json(directory/'SUBMISSION.json',submission)
    q=ssh(plan['hosts'][0],command,check=False,timeout=2*plan['timeout_seconds']+900)
    submission.update(submission_stdout=q.stdout,submission_stderr=q.stderr,foreground_ssh_exit_code=q.returncode); write_json(directory/'SUBMISSION.json',submission)
    if q.returncode!=0:
        final=load(directory/'workflow.json') if (directory/'workflow.json').exists() else {}
        if final.get('state') in ('complete','failed'):
            return final
        raise ValueError('foreground SSH failed with exit '+str(q.returncode)+'; inspect SUBMISSION.json; no completion inferred')
    if wait:
        while True:
            state=inspect(directory)
            if state.get('state') in ('complete','failed'):
                return state
            time.sleep(1)
    return {'state':'submitted','output':str(directory),'controller_host':plan['hosts'][0],'unit':unit}


def inspect(directory):
    directory=Path(directory)
    result=load(directory/'workflow.json') if (directory/'workflow.json').exists() else {'state':'submitted','output':str(directory)}
    sub=load(directory/'SUBMISSION.json') if (directory/'SUBMISSION.json').exists() else None
    if sub:
        q=ssh(sub['host'],['systemctl','--user','show',sub['unit'],'-p','ActiveState','-p','SubState','-p','Result','-p','ExecMainStatus'],check=False)
        result['controller_service_state']=q.stdout
        if result['state'] in ('running','submitted') and ('ActiveState=failed' in q.stdout or 'ActiveState=inactive' in q.stdout or q.returncode):
            result.update(state='failed',error='controller service ended before workflow completion')
    return result


def preparation_command(args,nsoc,soc,dest):
    python=str(absolute(args.environment_prefix,'environment_prefix')/'bin/python')
    command=[python,str(REPO/'cpu/geometry/prepare.py'),'--nsoc-input',str(nsoc),'--soc-input',str(soc),
             '--output',str(dest),'--group',args.group,'--geometry-tol',str(args.geometry_tol)]
    if args.symprec_ladder:
        command+=['--symprec-ladder',*[str(x) for x in args.symprec_ladder]]
    if args.project_nsoc:
        command+=['--project-nsoc']
    return command


def prepare(args):
    """Delegate to the repository-local common selected-subgroup generator."""
    nsoc=absolute(args.nsoc_input,'NSOC input',file=True); soc=absolute(args.soc_input,'SOC input',file=True)
    validate_inputs(nsoc.read_text(),soc.read_text())
    dest=absolute(args.output,'metadata output',exists=False)
    if dest.exists(): raise ValueError('metadata output already exists')
    command=preparation_command(args,nsoc,soc,dest)
    # The generator owns creation/fail-closed PREPARED.json. Do not pre-create
    # its output or leave old enabled settings from a failed preparation.
    q=subprocess.run(command,check=True,stdin=subprocess.DEVNULL,capture_output=True,text=True)
    result=load(dest/'PREPARED.json')
    if result.get('status')!='prepared':
        raise ValueError('common metadata generator did not report prepared')
    write_json(dest/'WORKFLOW_PREPARATION.json',{'command':command,'stdout':q.stdout,'stderr':q.stderr,
                'scope':'metadata only; full runtime profiles, covariance and physical validation remain explicit'})
    return result


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='action',required=True)
    for action in ('plan','run'):
        p=sub.add_parser(action); p.add_argument('config',type=Path)
        if action=='run': p.add_argument('--wait',action='store_true')
    p=sub.add_parser('inspect'); p.add_argument('output',type=Path)
    p=sub.add_parser('_controller',help=argparse.SUPPRESS); p.add_argument('plan',type=Path)
    p=sub.add_parser('prepare'); p.add_argument('--nsoc-input',required=True,type=Path); p.add_argument('--soc-input',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path); p.add_argument('--group','--point-group',dest='group',default='auto',choices=['auto','C2','D2','C3','D3','C4','D4','C6','D6','Ci','Cs'])
    p.add_argument('--symprec-ladder',type=float,nargs='+'); p.add_argument('--geometry-tol',type=float,default=5e-5)
    p.add_argument('--project-nsoc',action='store_true',help='explicit NSOC real-space H/S projection; DM/EDM are checked, not projected by this hook'); p.add_argument('--environment-prefix',default='snapshots/home/mambaforge')
    args=parser.parse_args(argv)
    try:
        if args.action=='_controller': return controller(args.plan)
        if args.action=='plan': result=plan_config(load(args.config))
        elif args.action=='run': result=submit(load(args.config),wait=args.wait)
        elif args.action=='inspect': result=inspect(args.output)
        else: result=prepare(args)
        print(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False))
        return 1 if result.get('state')=='failed' else 0
    except (ValueError,KeyError,OSError,subprocess.SubprocessError) as exc:
        print('openmx workflow: '+str(exc),file=sys.stderr); return 2


if __name__=='__main__': sys.exit(main())
