"""Translate explicit input requests into checked, complete CPU profiles.

No shell evaluation or material-specific geometry modification. Configuration
is prepared before MPI starts; native input code checks the same contract.
"""
import copy
import importlib.util
import math
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GROUPS = ('D3', 'C3', 'D2', 'C2', 'D4', 'C4', 'D6', 'C6', 'Ci', 'Cs')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def request(values):
    prefix = 'scf.symmetry.'
    allowed = {'mode', 'tolerance', 'pointgroup', 'operationsfile', 'geometrypolicy',
               'kpointpolicy', 'densitypolicy', 'forcepolicy', 'hsexport', 'covariancetolerance', 'hamiltonianpolicy'}
    for key in values:
        if key.startswith(prefix) and key[len(prefix):] not in allowed:
            raise ValueError('unknown symmetry keyword: '+key)
    mode = values.get(prefix+'mode', 'off').lower()
    if mode not in ('off', 'auto', 'user'):
        raise ValueError('scf.Symmetry.Mode must be Off, Auto or User')
    for key, default in [('geometrypolicy', 'strict'), ('kpointpolicy', 'auto'),
                         ('densitypolicy', 'free'), ('forcepolicy', 'free')]:
        if values.get(prefix+key, default).lower() != default:
            raise ValueError('unsupported scf.Symmetry.'+key+'; supported value: '+default)
    group = values.get(prefix+'pointgroup', 'none')
    group = next((g for g in GROUPS if g.lower() == group.lower()), group.lower())
    if mode == 'user' and group not in GROUPS:
        raise ValueError('User requires a supported scf.Symmetry.PointGroup')
    if mode != 'user' and group not in ('none', 'auto'):
        raise ValueError('scf.Symmetry.PointGroup requires Mode User')
    tol = float(values.get(prefix+'tolerance', '1e-5'))
    cov = float(values.get(prefix+'covariancetolerance', '1e-10'))
    if not math.isfinite(tol) or not 0 < tol <= .02:
        raise ValueError('scf.Symmetry.Tolerance must be finite in (0,0.02] Angstrom')
    if not math.isfinite(cov) or not 1e-14 <= cov <= 1e-3:
        raise ValueError('scf.Symmetry.CovarianceTolerance must be in [1e-14,1e-3]')
    policy = values.get(prefix+'hsexport', 'raw').lower()
    hpolicy = values.get(prefix+'hamiltonianpolicy','free').lower()
    if hpolicy not in ('free','symmetrize'):
        raise ValueError('scf.Symmetry.HamiltonianPolicy must be Free or Symmetrize')
    if policy not in ('raw', 'symmetrize'):
        raise ValueError('scf.Symmetry.HSExport must be raw or symmetrize')
    if policy == 'raw' and cov > 1e-10:
        raise ValueError('relaxed covariance requires explicit HSExport Symmetrize')
    ops = values.get(prefix+'operationsfile', 'none')
    if ops.lower() != 'none' and mode != 'user':
        raise ValueError('OperationsFile requires Mode User')
    if mode == 'off' and policy != 'raw':
        raise ValueError('HSExport Symmetrize requires Auto or User')
    if hpolicy=='symmetrize' and (mode=='off' or policy!='symmetrize'):
        raise ValueError('HamiltonianPolicy Symmetrize requires Auto/User and HSExport Symmetrize')
    return {'mode': mode, 'group': group if mode == 'user' else 'auto',
            'tolerance': tol, 'covariance_tolerance': cov, 'export_policy': policy,
            'operations_file': None if ops.lower() == 'none' else ops, 'hamiltonian_policy':hpolicy}


def validate_runtime_request(q, env, stage):
    if q['mode'] == 'off':
        return
    if env.get('OPENMX_CPU_INPUT_SYMMETRY_MODE') != q['mode']:
        raise ValueError('Auto/User requires matching prepared profiles: run workflow.py configure')
    group = env.get('OPENMX_CPU_INPUT_SYMMETRY_GROUP')
    if group not in GROUPS or (q['mode'] == 'user' and group != q['group']):
        raise ValueError('prepared subgroup disagrees with input PointGroup')
    if float(env.get('OPENMX_CPU_INPUT_SYMMETRY_TOL', 'nan')) != q['tolerance']:
        raise ValueError('prepared geometry tolerance disagrees with input')
    if env.get('OPENMX_CPU_INPUT_SYMMETRY_EXPORT') != q['export_policy']:
        raise ValueError('prepared H/S export disagrees with input')
    if stage=='nsoc':
        mode=env.get('OPENMX_CPU_REALSPACE_D3','off')
        if (mode in ('scf','1')) != (q['hamiltonian_policy']=='symmetrize'):
            raise ValueError('prepared NSOC Hamiltonian projection disagrees with input')
        expected='scf' if q['hamiltonian_policy']=='symmetrize' else 'export' if q['export_policy']=='symmetrize' else 'off'
        actual={'0':'off','1':'scf'}.get(mode,mode)
        if actual!=expected:
            raise ValueError('prepared NSOC export/output mode disagrees with input HSExport')
        if expected!='off' and not Path(env.get('OPENMX_CPU_REALSPACE_MAP','')).is_file():
            raise ValueError('prepared NSOC export requires a real-space map')
    covkey = 'OPENMX_CPU_NSOC_GENERIC_COV_TOL' if stage == 'nsoc' else 'OPENMX_CPU_GSOC_COV_TOL'
    if float(env.get(covkey, 'nan')) != q['covariance_tolerance']:
        raise ValueError('prepared covariance tolerance disagrees with input')
    key = 'OPENMX_CPU_NSOC_SYM' if stage == 'nsoc' else 'OPENMX_CPU_GSOC'
    basis = 'OPENMX_CPU_NSOC_BASIS' if stage == 'nsoc' else 'OPENMX_CPU_GSOC_BASIS'
    if env.get(key) != '1' or not Path(env.get(basis, '')).is_file():
        raise ValueError('prepared profile does not activate requested native '+stage+' symmetry')


def return_optimizations(group):
    # Native fused/selected return currently certifies only the real D3 v3
    # representation. Metadata generality does not extend that kernel's scope.
    return {'OPENMX_CPU_NSOC_FUSED_RETURN': '1' if group == 'D3' else '0',
            'OPENMX_CPU_NSOC_SELECTED_BACK': '1' if group == 'D3' else '0'}


def artifact_dependencies(config):
    result={}
    for path,expected in config.get('symmetry_artifacts_sha256',{}).items():
        p=Path(path)
        if not p.is_absolute() or not p.is_file():
            raise ValueError('prepared artifact missing: '+path)
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
            raise ValueError('prepared artifact changed since configure: '+path)
        result[str(p.resolve())]=expected
    return result


def configure(config, destination, wf):
    """Materialize metadata + full runnable config, without launching MPI."""
    config = copy.deepcopy(config)
    dest = wf.absolute(destination, 'configuration directory', exists=False)
    if dest.exists():
        raise ValueError('configuration directory already exists')
    paths = {s: wf.absolute(config['inputs'][s], s+' input', file=True) for s in ('nsoc', 'soc')}
    n, s = wf.validate_inputs(paths['nsoc'].read_text(), paths['soc'].read_text())
    q, qs = request(n), request(s)
    if q != qs:
        raise ValueError('NSOC/SOC symmetry input settings must match')
    if q['mode'] == 'off':
        raise ValueError('configure requires input scf.Symmetry.Mode Auto or User')
    if config.get('profiles') or config.get('symmetry'):
        raise ValueError('configure reads symmetry from inputs; omit manual profiles/symmetry')
    ops = q['operations_file']
    if ops:
        ops = str((paths['nsoc'].parent/ops).resolve()) if not Path(ops).is_absolute() else ops
        if not Path(ops).is_file():
            raise ValueError('OperationsFile does not exist')
        other = (paths['soc'].parent/qs['operations_file']).resolve()
        if Path(ops).resolve() != other:
            raise ValueError('NSOC/SOC OperationsFile must resolve to the same file')
    # All expensive geometry/basis work occurs on the caller's compute node.
    dest.mkdir(parents=True)
    prep = module('input_symmetry_preparation', ROOT/'geometry/prepare.py')
    data = prep.prepare(paths['nsoc'], paths['soc'], dest/'metadata', q['group'],
                        atom_tol=q['tolerance'], project_nsoc=q['hamiltonian_policy']=='symmetrize', soc_route='generic', export_policy=q['export_policy'],
                        detection_mode=q['mode'], operations_file=ops)
    group = data['selected_group']
    common = {'OPENMX_CPU_INPUT_SYMMETRY_MODE': q['mode'],
              'OPENMX_CPU_INPUT_SYMMETRY_GROUP': group,
              'OPENMX_CPU_INPUT_SYMMETRY_TOL': format(q['tolerance'], '.17g'),
              'OPENMX_CPU_INPUT_SYMMETRY_EXPORT': q['export_policy']}
    nenv = {**data['NSOC_environment_additions'], **common,
            'OPENMX_CPU_NSOC_GENERIC_COV_TOL': format(q['covariance_tolerance'], '.17g'),
            'OPENMX_CPU_NSOC_WARM_STEPS': '0', 'OPENMX_CPU_NSOC_RESET_DIIS': '0',
            'OPENMX_CPU_NSOC_RAW_RESET_AT': '0', 'OPENMX_CPU_DM_BLAS3': '1',
            'OPENMX_CPU_NSOC_SPARSE_U': '1', 'OPENMX_CPU_NSOC_SPARSE_BACK': '1',
            **return_optimizations(group)}
    # Output projection is separate from the SCF sector solve and does not
    # rewrite density or silently change convergence controls.
    if q['export_policy'] == 'symmetrize' and q['hamiltonian_policy']=='free':
        nenv['OPENMX_CPU_REALSPACE_D3'] = 'export'
    senv = {**data['SOC_environment_additions'], **common,
            'OPENMX_CPU_GSOC_PROJECT': '1' if q['export_policy'] == 'symmetrize' else '0',
            'OPENMX_CPU_GSOC_COV_TOL': format(q['covariance_tolerance'], '.17g'),
            'OPENMX_CPU_GSOC_GEOMETRY_TOL_ANG': format(q['tolerance'], '.17g'),
            'OPENMX_CPU_DM_BLAS3': '1'}
    config['profiles'] = {s: {'name': 'input-'+q['mode']+'-'+group+'-'+s, 'environment': e}
                          for s, e in [('nsoc', nenv), ('soc', senv)]}
    config['symmetry'] = {'project': q['export_policy'] == 'symmetrize', 'compact': False,
                          'soc_export_policy': q['export_policy'],
                          'geometry_tolerance_ang': q['tolerance'],
                          'covariance_tolerance': q['covariance_tolerance']}
    config['prepared_symmetry'] = str(dest/'metadata/PREPARED.json')
    config['symmetry_artifacts_sha256'] = {str(Path(path).resolve()):wf.sha256(path)
        for path in [*data['files'].values(),config['prepared_symmetry']]}
    config['symmetry_inputs_sha256'] = {stage:wf.sha256(path) for stage,path in paths.items()}
    config['symmetry_operations_sha256'] = {ops:wf.sha256(ops)} if ops else {}
    plan = wf.plan_config(config)
    wf.write_json(dest/'workflow.json', config)
    wf.write_json(dest/'plan.json', plan)
    return {'state': 'configured', 'config': str(dest/'workflow.json'), 'selected_group': group,
            'mode': q['mode'], 'geometry_modified': False, 'native_calculation_started': False}
