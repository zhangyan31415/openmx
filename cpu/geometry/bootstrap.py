#!/usr/bin/env python3
"""Internal rank-zero preparation for a native symmetry-enabled launch.

The actual OpenMX input is never rewritten. The paired canonical files below
are metadata-only views for the existing checked preparation implementation;
neither file is submitted to OpenMX. The final environment file is an atomic
success marker, with literal KEY=VALUE records and no shell interpretation.
"""
from __future__ import annotations

import argparse
import contextlib
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import socket
import sys

# This helper is serial even when its parent MPI rank owns a threaded BLAS
# environment. Set these before importing the geometry/numerical modules.
for _name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
              'NUMEXPR_NUM_THREADS', 'BLIS_NUM_THREADS'):
    os.environ[_name] = '1'

HERE = Path(__file__).resolve().parent
CPU = HERE.parent


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_text(path, text):
    temporary = path.with_name(path.name+'.tmp-'+str(os.getpid()))
    try:
        with temporary.open('x', encoding='utf-8') as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def canonical_text(text, replacements):
    """Retain geometry, basis and physical settings; normalize bookkeeping."""
    result = ['# METADATA ONLY: never run this file with OpenMX.',
              '# Actual native input is unchanged; see BOOTSTRAP.json.']
    active = None
    for raw in text.splitlines():
        line = raw.split('#', 1)[0].strip()
        if active is not None:
            result.append(raw)
            if line.lower() == active+'>':
                active = None
            continue
        if line.startswith('<'):
            active = line[1:].strip().lower()
            result.append(raw)
            continue
        key = line.split(None, 1)[0].lower() if line else ''
        if key in replacements or key == 'scf.symmetry' or key.startswith('scf.symmetry.'):
            continue
        result.append(raw)
    result.extend(key+' '+value for key, value in replacements.items())
    return '\n'.join(result)+'\n'


def native_logical(value):
    """Match the stock Inputtools.c logical spelling contract."""
    value = value.lower()
    if value in ('on', 'yes', 'true', '.true.', 'ok'):
        return True
    if value in ('off', 'no', 'false', '.false.', 'ng'):
        return False
    raise ValueError('invalid OpenMX logical value: '+value)


def check_scope(values, blocks, stage):
    if stage not in ('nsoc', 'soc'):
        raise ValueError('stage must be nsoc or soc')
    if (values.get('scf.eigenvaluesolver', '').lower() not in ('cluster', 'band') or
            values.get('scf.kgrid', '1 1 1').split() != ['1', '1', '1']):
        raise ValueError('native symmetry currently requires Gamma Cluster/Band')
    if values.get('scf.eigen.lib', '').lower() not in ('elpa1', 'elpa2'):
        raise ValueError('native symmetry requires scf.eigen.lib elpa1 or elpa2')
    if values.get('md.type', 'nomd').lower() not in ('nomd', 'none'):
        raise ValueError('native symmetry requires fixed geometry/no MD')
    if (values.get('orbitalopt.method','off').lower() != 'off' or
            native_logical(values.get('scf.mixing.control.temp','off'))):
        raise ValueError('native symmetry requires uncontracted orbitals and fixed electronic temperature')
    spin = values.get('scf.spinpolarization', 'off').lower()
    spinorbit = native_logical(values.get('scf.spinorbit.coupling', 'off'))
    if stage == 'nsoc' and (spin != 'off' or spinorbit):
        raise ValueError('NSOC preparation requires nonmagnetic scalar input')
    if stage == 'nsoc' and values.get('scf.restart','off').lower() not in ('off','no','false','.false.','ng'):
        raise ValueError('native symmetry NSOC requires a fresh calculation with restart off')
    if stage == 'soc' and (spin != 'nc' or not spinorbit or
            values.get('scf.restart', 'off').lower() != 'c2n' or
            int(values.get('scf.maxiter', '40')) != 1):
        raise ValueError('SOC symmetry currently supports only nonmagnetic one-step C2N SOC')
    for row in blocks.get('atoms.speciesandcoordinates', []):
        if len(row) < 7:
            raise ValueError('explicit equal nonmagnetic initial atom populations required')
        up, down = float(row[5]), float(row[6])
        if not (abs(up-down) <= 1e-12):
            raise ValueError('native symmetry requires nonmagnetic initial atom populations')


def fast_memory_scope(values):
    """Conservative intersection of the verified sparse/cache native guards."""
    off = ('off','no','false','.false.','ng','0','none')
    if (values.get('scf.eigenvaluesolver', '').lower() not in ('cluster','band') or
            values.get('scf.mixing.type', 'simple').lower() != 'rmm-diis' or
            values.get('orbitalopt.method','off').lower() != 'off' or
            not native_logical(values.get('scf.proexpn.vna','on')) or
            int(values.get('level.of.fileout','1')) >= 2):
        return False
    blocked = ('scf.hubbard.u','scf.constraint.nc.spin','scf.nc.zeeman.spin',
               'scf.nc.zeeman.orbital','scf.core.hole','xmcd.calc','partial.charge',
               'cddf.start','lnao.calc.flag','lnbo.calc.flag','nbo.switch',
               'voronoi.charge','voronoi.orbital.moment','dos.fileout','dosgauss.fileout',
               'lno.flag','cwf.calc','mo.fileout','xanes.calc','dm.export',
               'scf.stress.tensor','scf.mixing.control.temp')
    if any(values.get(key,'off').lower() not in off for key in blocked):
        return False
    return all(float(x.replace('D','E').replace('d','e')) == 0
               for x in values.get('scf.electric.field','0 0 0').split())


def preferred_soc_route(values):
    # Explicit advanced covariance requests retain the generic route and its
    # declared diagnostics/guards rather than switching to the legacy route.
    if (values.get('scf.eigenvaluesolver', '').lower() != 'cluster' or not fast_memory_scope(values) or
            'scf.symmetry.covariancetolerance' in values):
        return 'generic'
    return 'auto'


def profiles(prepared, q, api, *, threads=1, values=None, ranks_single_cpu_pinned=False):
    group = prepared['selected_group']
    templates = json.loads((CPU/'tools/examples/profile_templates.json').read_text())['templates']
    scalar_template = 'D3_legacy_fast_candidate' if group == 'D3' else 'C3_gsoc'
    def defaults(template, stage):
        return {k:v for k,v in templates[template]['profiles'][stage]['environment'].items()
                if '${' not in v}
    common = {'OPENMX_CPU_INPUT_SYMMETRY_MODE': q['mode'],
              'OPENMX_CPU_INPUT_SYMMETRY_GROUP': group,
              'OPENMX_CPU_INPUT_SYMMETRY_TOL': format(q['tolerance'], '.17g'),
              'OPENMX_CPU_INPUT_SYMMETRY_EXPORT': q['export_policy'],
              'OPENMX_CPU_INPUT_SYMMETRY_COV_TOL': format(q['covariance_tolerance'], '.17g'),
              'OPENMX_CPU_DM_BLAS3': '1'}
    nenv = {**defaults(scalar_template, 'nsoc'), **prepared['NSOC_environment_additions'], **common,
            'OPENMX_CPU_NSOC_GENERIC_COV_TOL': format(q['covariance_tolerance'], '.17g'),
            'OPENMX_CPU_NSOC_WARM_STEPS': '0', 'OPENMX_CPU_NSOC_RESET_DIIS': '0',
            'OPENMX_CPU_NSOC_RAW_RESET_AT': '0',
            'OPENMX_CPU_NSOC_PHASE': '0',
            'OPENMX_CPU_NSOC_SPARSE_U': '1', 'OPENMX_CPU_NSOC_SPARSE_BACK': '1',
            **api.return_optimizations(group)}
    # Python stays serial; select native caches from the actual parent -nt.
    nenv['OPENMX_CPU_GRID_NUMA'] = '1' if threads == 1 and ranks_single_cpu_pinned else '0'
    if threads > 1:
        nenv['OPENMX_CPU_VNA_SHARED'] = '0'
    values = values or {}
    fast_memory = fast_memory_scope(values)
    if not fast_memory:
        safe = set(prepared['NSOC_environment_additions']) | set(common) | set(api.return_optimizations(group))
        safe.update('OPENMX_CPU_NSOC_'+key for key in
                    ('GENERIC_COV_TOL','WARM_STEPS','RESET_DIIS','RAW_RESET_AT','PHASE','SPARSE_U','SPARSE_BACK'))
        nenv = {k:v for k,v in nenv.items() if k in safe}
    if int(values.get('scf.mixing.history','5')) != 30 or q['hamiltonian_policy'] != 'symmetrize':
        nenv['OPENMX_CPU_NSOC_SYM_HISTORY'] = '0'
    if q['hamiltonian_policy'] != 'symmetrize':
        nenv.update(OPENMX_CPU_NSOC_REPRESENTATIVE_DM='0', OPENMX_CPU_NSOC_REP_TAIL_BOUND='0')
    if q['export_policy'] == 'symmetrize' and q['hamiltonian_policy'] == 'free':
        nenv['OPENMX_CPU_REALSPACE_D3'] = 'export'
    if prepared['SOC_route'] == 'legacy':
        senv = {**defaults('D3_legacy_fast_candidate', 'soc'),
                **prepared['SOC_environment_additions'], **common}
        senv.update(OPENMX_CPU_GSOC='0', OPENMX_CPU_GSOC_COMPACT='0')
        if group == 'C3':
            senv.update(OPENMX_CPU_D3='0', OPENMX_CPU_SOC_TR='0', OPENMX_CPU_SOC_REAL='0')
    else:
        template = group+'_gsoc' if group in ('D2','C2') else 'C3_gsoc'
        senv = {**defaults(template, 'soc'), **prepared['SOC_environment_additions'], **common,
                'OPENMX_CPU_GSOC_PROJECT': '1' if q['export_policy'] == 'symmetrize' else '0',
                'OPENMX_CPU_GSOC_COV_TOL': format(q['covariance_tolerance'], '.17g'),
                'OPENMX_CPU_GSOC_GEOMETRY_TOL_ANG': format(q['tolerance'], '.17g')}
    senv['OPENMX_CPU_GRID_NUMA'] = '0'
    if not fast_memory:
        senv = {**prepared['SOC_environment_additions'], **common,
                'OPENMX_CPU_GSOC_PROJECT': '1' if q['export_policy']=='symmetrize' else '0',
                'OPENMX_CPU_GSOC_COV_TOL': format(q['covariance_tolerance'], '.17g'),
                'OPENMX_CPU_GSOC_GEOMETRY_TOL_ANG': format(q['tolerance'], '.17g')}
    if int(values.get('scf.mixing.startpulay','6')) < 4:
        senv['OPENMX_CPU_ONESTEP_HISTORY'] = '0'
    return {'nsoc': nenv, 'soc': senv}


def report_symmetry(prepared, stage, request):
    """Report the operations actually selected for this solver, once on rank 0."""
    selected = json.loads(Path(prepared['files']['selected_symmetry']).read_text())['selected']
    operations = selected['operations']
    print('CPU_SYMMETRY_SUMMARY stage='+stage+' selected_group='+prepared['selected_group']+
          ' detected_point_group='+str(prepared['detected_point_group']['symbol'])+
          ' operations='+str(len(operations))+
          ' geometry_max_Ang='+format(selected['max_geometry_residual_Ang'], '.12e')+
          ' H_policy='+request['hamiltonian_policy']+' HS_export='+request['export_policy'], flush=True)

    def vector(values):
        return ','.join(format(v, '.12g') for v in values)

    def matrix(values):
        return '['+';'.join(vector(row) for row in values)+']'

    print('CPU_SYMMETRY_CONVENTION fractional_column x_prime=W*x+t; R_cart acts on Cartesian columns; atom_maps_in=metadata/selected.symmetry.json', flush=True)
    for index, operation in enumerate(operations, 1):
        print('CPU_SYMMETRY_OPERATION index='+str(index)+'/'+str(len(operations))+
              ' W='+matrix(operation['W'])+' t_frac=['+vector(operation['t'])+']'+
              ' R_cart='+matrix(operation['R_cart_orthogonal'])+
              ' residual_Ang='+format(operation['residual_max_Ang'], '.12e'), flush=True)


def bootstrap(input_path, output, stage, environment_file, *, runtime_cwd=None, threads=1, ranks_single_cpu_pinned=False):
    if isinstance(threads, bool) or not isinstance(threads, int) or threads < 1:
        raise ValueError('native thread count must be a positive integer')
    input_path = Path(input_path).resolve(strict=True)
    out = Path(output).resolve(strict=True)
    envfile = Path(environment_file).absolute()
    cwd = Path.cwd().resolve() if runtime_cwd is None else Path(runtime_cwd).resolve(strict=True)
    if not out.is_dir() or envfile.parent.resolve() != out or envfile.is_symlink() or envfile.exists():
        raise ValueError('fresh environment file directly inside the task-owned output directory required')
    if any(any(c in str(p) for c in '\n\r\0') for p in (input_path, out, envfile, cwd)):
        raise ValueError('input/output paths cannot contain newline, CR or NUL')
    if (out/'BOOTSTRAP.json').exists() or (out/'canonical').exists() or (out/'metadata').exists():
        raise ValueError('bootstrap output already contains preparation artifacts; use a fresh directory')
    provenance = {'schema': 1, 'status': 'preparing', 'stage': stage,
                  'actual_input': str(input_path), 'actual_input_sha256': sha(input_path),
                  'actual_input_modified': False, 'runtime_cwd': str(cwd),
                  'hostname': socket.gethostname(), 'python': sys.executable,
                  'environment_prefix': sys.prefix, 'command': list(sys.argv),
                  'helper_sha256': sha(__file__), 'helper_threads': 1, 'parent_threads': threads,
                  'ranks_single_cpu_pinned': bool(ranks_single_cpu_pinned),
                  'native_calculation_started': False,
                  'canonical_inputs_are_metadata_only': True}
    try:
        workflow = module('bootstrap_workflow', CPU/'tools/workflow.py')
        api = module('bootstrap_input_symmetry', CPU/'tools/input_symmetry.py')
        text = input_path.read_text()
        values, blocks = workflow.input_values(text)
        q = api.request(values)
        if q['mode'] == 'off':
            raise ValueError('native preparation requires scf.Symmetry on or explicit Auto/User')
        check_scope(values, blocks, stage)
        root = Path(values.get('data.path', '../DFT_DATA19'))
        data = root.resolve() if root.is_absolute() else (cwd/root).resolve()
        if not data.is_dir():
            raise ValueError('DATA.PATH must resolve to a directory relative to the launch working directory')
        ops = q.get('operations_file')
        if ops:
            opath = Path(ops)
            ops = str((cwd/opath).resolve(strict=True) if not opath.is_absolute() else opath.resolve(strict=True))
        canonical = out/'canonical'
        canonical.mkdir(mode=0o700)
        # Strict paired-input validation is reused only for metadata. Native
        # defaults are made explicit here; they never change the actual input.
        common = {'data.path': str(data), 'system.currrentdirectory': './',
                  'hs.fileout': 'on', 'scf.xctype': values.get('scf.xctype', 'LDA'),
                  'scf.energycutoff': values.get('scf.energycutoff', '220'),
                  'scf.electronictemperature': values.get('scf.electronictemperature', '300')}
        n = {**common, 'system.name': 'metadata_nsoc', 'scf.spinpolarization': 'off',
             'scf.spinorbit.coupling': 'off', 'scf.restart': 'off',
             'scf.restart.filename': 'metadata_nsoc',
             'scf.maxiter': str(max(2, int(values.get('scf.maxiter', '40'))))}
        s = {**common, 'system.name': 'metadata_soc', 'scf.spinpolarization': 'NC',
             'scf.spinorbit.coupling': 'on', 'scf.restart': 'C2N',
             'scf.restart.filename': 'metadata_nsoc', 'scf.maxiter': '1'}
        npath, spath = canonical/'nsoc.dat', canonical/'soc.dat'
        npath.write_text(canonical_text(text, n))
        spath.write_text(canonical_text(text, s))
        prepare = module('bootstrap_prepare', HERE/'prepare.py')
        preparation_log = out/'preparation.log'
        with preparation_log.open('x') as stream, contextlib.redirect_stdout(stream):
            prepared = prepare.prepare(npath, spath, out/'metadata', q['group'],
                        atom_tol=q['tolerance'], project_nsoc=q['hamiltonian_policy'] == 'symmetrize',
                        soc_route=preferred_soc_route(values) if q['export_policy']=='symmetrize' else 'generic', export_policy=q['export_policy'],
                        detection_mode=q['mode'], operations_file=ops)
        env = profiles(prepared, q, api, threads=threads, values=values, ranks_single_cpu_pinned=ranks_single_cpu_pinned)[stage]
        export_requested=native_logical(values.get('hs.fileout','off'))
        if stage=='soc' and prepared['SOC_route']=='generic' and not export_requested:
            # Preserve the projected solve without activating a file-export
            # hook when the native input did not request an H/S file.
            env['OPENMX_CPU_GSOC_EXPORT_SYMMETRIZE']='0'
        api.validate_runtime_request(q, env, stage)
        for key, value in env.items():
            if (not re.fullmatch(r'OPENMX_CPU_[A-Z0-9_]+', key) or not isinstance(value, str) or
                    any(c in value for c in '\n\r\0')):
                raise ValueError('unsafe native environment record: '+key)
        if sha(input_path) != provenance['actual_input_sha256']:
            raise ValueError('actual input changed while symmetry metadata was being prepared')
        provenance.update(status='prepared', request=q, selected_group=prepared['selected_group'], SOC_route=prepared['SOC_route'],
                          actual_HS_export_requested=export_requested,
                          environment=env, environment_file=str(envfile),
                          canonical_inputs={k: {'path': str(p), 'sha256': sha(p)}
                                            for k, p in [('nsoc', npath), ('soc', spath)]},
                          artifacts_sha256=prepared['sha256'],
                          preparation_log=str(preparation_log),
                          prepared_manifest=str(out/'metadata/PREPARED.json'))
        atomic_text(out/'BOOTSTRAP.json', json.dumps(provenance, indent=2)+'\n')
        report_symmetry(prepared, stage, q)
        # This is deliberately the last successful write. C reads only this
        # final file after a successful helper exit; partial metadata is inert.
        atomic_text(envfile, ''.join(key+'='+env[key]+'\n' for key in sorted(env)))
        print('OpenMX symmetry prepared: stage='+stage+' subgroup='+prepared['selected_group']+
              ' H/S='+q['export_policy']+' metadata='+str(out), flush=True)
        return provenance
    except Exception as exc:
        provenance.update(status='failed', error=str(exc), native_launch_allowed=False)
        atomic_text(out/'BOOTSTRAP.json', json.dumps(provenance, indent=2)+'\n')
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=('nsoc', 'soc'), required=True)
    parser.add_argument('--environment-file', type=Path, required=True)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--rank-pinned', type=int, choices=(0,1), default=0)
    args = parser.parse_args()
    try:
        bootstrap(args.input, args.output, args.stage, args.environment_file, threads=args.threads,
                  ranks_single_cpu_pinned=bool(args.rank_pinned))
    except (OSError, ValueError, RuntimeError, ImportError) as exc:
        print('OpenMX symmetry preparation failed: '+str(exc), file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
