#!/usr/bin/env python3
"""Internal supervised node memory monitor and MPI-preserving rank wrapper."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import socket
import select
import subprocess
import sys
import time


def atomic_json(path,data):
    path=Path(path); tmp=path.with_name(path.name+'.writing-'+str(os.getpid()))
    tmp.write_text(json.dumps(data,indent=2)+'\n'); tmp.replace(path)


def rank_identity(env):
    try:
        rank=int(env['PMI_RANK'])
        if 'MPI_LOCALRANKID' in env:
            local=int(env['MPI_LOCALRANKID'])
        else:
            local=int(env['PMI_LOCAL_RANK'])
        if rank<0 or local<0: raise ValueError()
        return rank,local
    except (KeyError,ValueError):
        raise ValueError('Intel MPI rank/local-rank identity missing or invalid') from None


def rank_cpus(plan,host,local):
    if local<0 or local>=plan['ranks_per_host']:
        raise ValueError('local rank outside declared host layout')
    t=plan['threads']; ids=plan['layout']['cpu_ids_per_host'][host][local*t:(local+1)*t]
    if len(ids)!=t: raise ValueError('rank CPU binding unavailable')
    return ids


def task_slice_path(text):
    path=next((line.split('::',1)[1] for line in text.splitlines() if line.startswith('0::')),None)
    if not path or not path.startswith('/user.slice/'):
        raise ValueError('task must be in delegated user cgroup v2 slice')
    parts=Path(path).parts
    for i,part in enumerate(parts):
        if part.startswith('omxwf') and part.endswith('.slice'):
            return Path('/sys/fs/cgroup').joinpath(*parts[1:i+1])
    raise ValueError('dedicated workflow slice not found')


def cpu_topology(cpus, root=Path('/sys/devices/system/cpu')):
    result={}
    seen=set()
    for cpu in cpus:
        p=root/('cpu'+str(cpu))/'topology'
        pair=(int((p/'physical_package_id').read_text()),int((p/'core_id').read_text()))
        if pair in seen:
            raise ValueError('CPU list repeats a physical core through SMT threads')
        seen.add(pair)
        result[str(cpu)]={'socket':pair[0],'core':pair[1]}
    return result


def monitor(args):
    actual=socket.gethostname().split('.')[0]
    if actual!=args.host.split('.')[0]: raise ValueError('monitor started on wrong host')
    directory=args.directory; dest=directory/'nodes'/args.host; dest.mkdir(parents=True,exist_ok=False)
    cg=task_slice_path(Path('/proc/self/cgroup').read_text())
    # Require kernel physical accounting, never substitute summed process RSS.
    for name in ('memory.current','memory.peak','memory.events'):
        if not (cg/name).is_file(): raise ValueError('kernel cgroup accounting unavailable: '+name)
    plan_path=directory.parent/'plan.json'
    plan=json.loads(plan_path.read_text())
    topology=cpu_topology(plan['layout']['cpu_ids_per_host'][args.host])
    ready={'cpu_topology':topology,'host':actual,'requested_host':args.host,'pid':os.getpid(),'cgroup':str(cg),'time_unix':time.time(),
           'category':'task slice physical memory includes all own MPI scopes and monitor; native shared allocations charged once'}
    atomic_json(dest/'ready.json',ready)
    peaks=0; count=0
    def sample():
        nonlocal peaks,count
        now=time.time(); cur=int((cg/'memory.current').read_text()); peak=int((cg/'memory.peak').read_text())
        peaks=max(peaks,peak,cur); count+=1
        mem={}
        for line in Path('/proc/meminfo').read_text().splitlines():
            if line.startswith(('MemTotal:','MemAvailable:')):
                key,value,*_=line.split(); mem[key.rstrip(':')]=int(value)*1024
        events={line.split()[0]:int(line.split()[1]) for line in (cg/'memory.events').read_text().splitlines()}
        return {'sample_id':math.floor(now/args.interval),'time_unix':now,'host':args.host,'slice_current_bytes':cur,
                'slice_peak_bytes':peak,'node_mem_total_bytes':mem.get('MemTotal'),'node_mem_available_bytes':mem.get('MemAvailable'),
                'memory_events':events}
    with (dest/'samples.jsonl').open('w',buffering=1) as out:
        while not (directory/'MONITOR_STOP').exists():
            out.write(json.dumps(sample())+'\n')
            # Align all hosts to the next shared wall-clock sample boundary.
            now=time.time(); time.sleep(max(.001,(math.floor(now/args.interval)+1)*args.interval-now))
        final=sample(); out.write(json.dumps(final)+'\n')
    remaining=set()
    for procs in cg.rglob('cgroup.procs'):
        try: remaining.update(int(x) for x in procs.read_text().split())
        except FileNotFoundError: pass
    remaining.discard(os.getpid())
    atomic_json(dest/'final.json',{**ready,**final,'slice_peak_bytes':peaks,'samples':count,
                                 'remaining_descendant_pids':sorted(remaining),'state':'complete'})
    return 0


def run_scope_command(command):
    # Hydra passes a live PMI_FD; closing descriptors here breaks MPI_Init.
    return subprocess.run(command, stdin=subprocess.DEVNULL, close_fds=False, timeout=None).returncode


def rank(args):
    plan=json.loads(args.plan.read_text()); rankno,local=rank_identity(os.environ)
    token=hashlib.sha256(plan['output'].encode()).hexdigest()[:16]+args.stage
    command=['systemd-run','--user','--scope','--quiet','--unit=omxwf'+token+'rank'+str(rankno),
             '--slice=omxwf'+token+'.slice','--property=MemoryAccounting=yes','--property=TasksMax=infinity',
             '--property=RuntimeMaxSec=infinity','--property=TimeoutStopSec=infinity',
             plan['environment_python'],str(Path(__file__).resolve()),'exec','--plan',str(args.plan),'--stage',args.stage]
    # scope mode inherits the original PMI environment and descriptors; service
    # mode would detach from Hydra and is deliberately not used for MPI ranks.
    return run_scope_command(command)


def binding_environment(inherited,cpus):
    env={k:v for k,v in inherited.items() if not k.startswith('KMP_') and k!='GOMP_CPU_AFFINITY'}
    # One explicitly selected physical CPU per OpenMP worker. Binding happens
    # before the native executable initializes Intel OpenMP/MKL, and native
    # sector guards still independently verify every actual worker mask.
    env['OMP_PLACES']=','.join('{'+str(cpu)+'}' for cpu in cpus)
    env['OMP_PROC_BIND']='close'
    return env


def execute(args):
    plan=json.loads(args.plan.read_text()); rankno,local=rank_identity(os.environ)
    actual=socket.gethostname().split('.')[0]
    candidates=[h for h in plan['hosts'] if h.split('.')[0]==actual]
    if len(candidates)!=1: raise ValueError('MPI rank appeared on unexpected/ambiguous host')
    host=candidates[0]; cpus=rank_cpus(plan,host,local)
    os.sched_setaffinity(0,cpus)
    env=binding_environment(os.environ,cpus)
    cg=task_slice_path(Path('/proc/self/cgroup').read_text())
    token=hashlib.sha256(plan['output'].encode()).hexdigest()[:16]+args.stage
    if cg.name!='omxwf'+token+'.slice': raise ValueError('rank outside this task memory accounting slice')
    directory=Path(plan['output'])/args.stage
    witness={'rank':rankno,'local_rank':local,'host':host,'actual_hostname':socket.gethostname(),'pid':os.getpid(),
             'cpu_affinity':sorted(os.sched_getaffinity(0)),'inside_task_slice':True,'cgroup':str(cg),
             'threads':plan['threads'],'cwd':str(Path.cwd()),'time_unix':time.time(),
             'thread_environment':{k:env.get(k) for k in ('OMP_NUM_THREADS','MKL_NUM_THREADS','OMP_DYNAMIC','MKL_DYNAMIC','OMP_PLACES','OMP_PROC_BIND')},
             'MPI_environment':{k:v for k,v in os.environ.items() if k.startswith(('PMI_','MPI_LOCAL','I_MPI_','FI_'))}}
    atomic_json(directory/'ranks'/('rank'+str(rankno)+'.json'),witness)
    if rankno==0:
        q=subprocess.run(['ldd',plan['binary']],stdin=subprocess.DEVNULL,text=True,capture_output=True)
        (directory/'LOADED_LIBRARY_RESOLUTION.txt').write_text(q.stdout+q.stderr)
    os.chdir(directory)
    os.execve(plan['binary'],[plan['binary'],'input.dat','-nt',str(plan['threads'])],env)


def anchor(args):
    directory=args.directory
    host=socket.gethostname().split('.')[0]
    if host!=args.host.split('.')[0]: raise ValueError('SSH anchor on unexpected host')
    atomic_json(directory/(args.host+'.ready.json'),{'host':host,'pid':os.getpid(),'state':'SSH login held','time_unix':time.time()})
    while not (directory/'STOP').exists():
        ready,_,_=select.select([sys.stdin],[],[],.25)
        if ready and not os.read(sys.stdin.fileno(),1): break
    atomic_json(directory/(args.host+'.final.json'),{'host':host,'pid':os.getpid(),'state':'closed','time_unix':time.time()})
    return 0


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('anchor');a.add_argument('--directory',type=Path,required=True);a.add_argument('--host',required=True)
    m=sub.add_parser('monitor');m.add_argument('--directory',type=Path,required=True);m.add_argument('--host',required=True);m.add_argument('--interval',type=float,required=True)
    for action in ('rank','exec'):
        r=sub.add_parser(action);r.add_argument('--plan',type=Path,required=True);r.add_argument('--stage',choices=['nsoc','soc'],required=True)
    args=p.parse_args()
    try:
        if args.action=='anchor': return anchor(args)
        if args.action=='monitor': return monitor(args)
        if args.action=='rank': return rank(args)
        return execute(args)
    except (ValueError,OSError,KeyError) as exc:
        print('openmx workflow node: '+str(exc),file=sys.stderr);return 2

if __name__=='__main__':sys.exit(main())
