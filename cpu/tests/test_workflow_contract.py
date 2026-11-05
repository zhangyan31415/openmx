import sys
"""Safety contracts of the public NSOC -> own-checkpoint -> SOC workflow."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('workflow', Path(__file__).parents[1] / 'tools/workflow.py')
if SPEC and SPEC.loader and Path(SPEC.origin).exists():
    wf = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(wf)
else:
    wf = None

BASE = '''System.Name scalar
System.CurrrentDirectory ./
DATA.PATH {data}
scf.XcType GGA-PBE
scf.energycutoff 150
scf.ElectronicTemperature 300
scf.criterion 1e-7
scf.SpinPolarization off
scf.SpinOrbit.Coupling off
scf.EigenvalueSolver Cluster
scf.Eigen.lib elpa2
scf.restart off
scf.maxIter 200
HS.fileout on
Atoms.Number 1
Atoms.UnitVectors.Unit Ang
Atoms.SpeciesAndCoordinates.Unit Ang
<Atoms.UnitVectors
10 0 0
0 10 0
0 0 10
Atoms.UnitVectors>
<Definition.of.Atomic.Species
H H5.0-s1 H_PBE19
Definition.of.Atomic.Species>
<Atoms.SpeciesAndCoordinates
1 H 0 0 0 0.5 0.5
Atoms.SpeciesAndCoordinates>
'''

class WorkflowContract(unittest.TestCase):
    def setUp(self):
        self.assertIsNotNone(wf, 'workflow implementation is missing')
        self.tmp = tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        src = self.root / 'inputs'; src.mkdir()
        data = self.root / 'DFT_DATA19'; data.mkdir()
        self.ntext = BASE.format(data=data)
        self.stext = self.ntext.replace('System.Name scalar', 'System.Name spinor').replace('scf.SpinPolarization off','scf.SpinPolarization nc').replace('scf.SpinOrbit.Coupling off','scf.SpinOrbit.Coupling on').replace('scf.restart off','scf.restart c2n\nscf.restart.filename scalar').replace('scf.maxIter 200','scf.maxIter 1')
        self.n = src / 'nsoc.dat'; self.n.write_text(self.ntext)
        self.s = src / 'soc.dat'; self.s.write_text(self.stext)
        self.binary = self.root / 'openmx'; self.binary.write_text('#!/bin/sh\nexit 0\n'); self.binary.chmod(0o755)
        self.mpi = self.root / 'mpirun'; self.mpi.write_text('#!/bin/sh\nexit 0\n'); self.mpi.chmod(0o755)
        self.config = {'schema':1,'binary':str(self.binary), 'inputs':{'nsoc':str(self.n),'soc':str(self.s)},'output':str(self.root/'results'),
          'layout':{'hosts':['node040'],'ranks_per_host':2,'threads':1,'cpu_ids_per_host':{'node040':[0,1]}},
          'runtime':{'mpiexec':str(self.mpi),'runtime_prefix':str(self.root),'environment_prefix':sys.prefix,'transport':'shm'},
          'profiles':{'nsoc':{'name':'raw','environment':{}},'soc':{'name':'raw','environment':{}}},'symmetry':{'project':False,'compact':False}}

    def test_plan_is_read_only_and_preserves_inputs(self):
        before = self.n.read_bytes(), self.s.read_bytes()
        p = wf.plan_config(self.config)
        self.assertEqual(p['ranks'],2); self.assertEqual(p['threads'],1)
        self.assertFalse((self.root/'results').exists())
        self.assertEqual(before, (self.n.read_bytes(),self.s.read_bytes()))

    def test_scratch_defaults_to_shared_output_parent_and_can_be_overridden(self):
        plan=wf.plan_config(self.config)
        self.assertTrue(Path(plan['scratch_root']).is_relative_to(self.root/'.openmx-scratch'))
        self.config['runtime']['shared_scratch_root']=str(self.root/'shared-temporary')
        plan=wf.plan_config(self.config)
        self.assertTrue(Path(plan['scratch_root']).is_relative_to(self.root/'shared-temporary'))

    def test_duplicate_json_keys_rejected(self):
        with self.assertRaisesRegex(ValueError,'duplicate'):
            wf.loads_unique('{"schema":1,"schema":2}')

    def test_duplicate_keyword_case_insensitive_rejected(self):
        self.s.write_text(self.stext+'\nSCF.MaxIter 2\n')
        with self.assertRaisesRegex(ValueError,'duplicate'):
            wf.plan_config(self.config)

    def test_nonfinite_json_rejected(self):
        with self.assertRaisesRegex(ValueError,'nonfinite'):
            wf.loads_unique('{"schema":1,"timeout_seconds":NaN}')

    def test_symmetric_route_requires_explicit_basis(self):
        self.config['profiles']['soc']['environment']={'OPENMX_CPU_GSOC':'1'}
        with self.assertRaisesRegex(ValueError,'basis|BASIS'):
            wf.plan_config(self.config)

    def test_false_covariance_tolerance_not_accepted(self):
        self.config['symmetry']['covariance_tolerance']=False
        with self.assertRaisesRegex(ValueError,'covariance'):
            wf.plan_config(self.config)

    def test_nonfinite_physics_rejected(self):
        self.n.write_text(self.ntext.replace('scf.ElectronicTemperature 300','scf.ElectronicTemperature NaN'))
        self.s.write_text(self.stext.replace('scf.ElectronicTemperature 300','scf.ElectronicTemperature NaN'))
        with self.assertRaisesRegex(ValueError,'nonfinite'):
            wf.plan_config(self.config)

    def test_stale_sector_partition_rejected(self):
        self.config['profiles']['nsoc']['environment']={'OPENMX_CPU_NSOC_SECTOR_PARALLEL':'1','OPENMX_CPU_NSOC_SECTOR_RANKS':'48,8,8'}
        with self.assertRaisesRegex(ValueError,'sector|rank'):
            wf.plan_config(self.config)

    def test_memory_limit_preserved_in_plan(self):
        self.config['memory_limit_bytes_per_host']=8*1024**3
        self.assertEqual(wf.plan_config(self.config)['memory_limit_bytes_per_host'],8*1024**3)

    def test_sector_solver_flag_is_not_density_map_projection(self):
        basis=self.root/'scalar.nsym';basis.write_bytes(b'basis')
        self.config['profiles']['nsoc']['environment']={'OPENMX_CPU_NSOC_SYM':'1','OPENMX_CPU_NSOC_FULL_PROJECT':'1','OPENMX_CPU_NSOC_BASIS':str(basis),'OPENMX_CPU_REALSPACE_D3':'off'}
        p=wf.plan_config(self.config)
        self.assertFalse(p['nsoc_density_projection_explicit'])
        self.assertTrue(p['nsoc_sector_solver_all_steps'])
        self.config['profiles']['nsoc']['environment']['OPENMX_CPU_REALSPACE_D3']='scf'
        p=wf.plan_config(self.config)
        self.assertFalse(p['nsoc_density_projection_explicit'])
        self.assertTrue(p['nsoc_realspace_HS_projection_explicit'])
        self.assertEqual(p['nsoc_export_HS'],'PROJECTED_REALSPACE_H_S')

    def test_foreground_ssh_failure_returns_without_polling_forever(self):
        from unittest.mock import patch
        from types import SimpleNamespace
        plan=wf.plan_config(self.config);plan['build_manifest']='/dummy/manifest'
        with patch.object(wf,'plan_config',return_value=plan),patch.object(wf,'ssh',return_value=SimpleNamespace(returncode=255,stdout='',stderr='connection lost')),patch.object(wf.time,'sleep',side_effect=RuntimeError('infinite polling')):
            with self.assertRaisesRegex(ValueError,'SSH.*failed'):
                wf.submit(self.config,wait=True)

    def test_vna_shared_rejected_for_threaded_layout_in_either_stage(self):
        self.config['layout']['threads']=4
        self.config['layout']['cpu_ids_per_host']['node040']=list(range(8))
        for stage in ('nsoc','soc'):
            self.config['profiles'][stage]['environment']={'OPENMX_CPU_VNA_SHARED':'1'}
            with self.assertRaisesRegex(ValueError,'VNA_SHARED.*thread'):
                wf.plan_config(self.config)
            self.config['profiles'][stage]['environment']={'OPENMX_CPU_VNA_SHARED':'0'}
        self.assertEqual(wf.plan_config(self.config)['threads'],4)

    def test_export_only_realspace_projects_hs_without_density_write(self):
        self.config['profiles']['nsoc']['environment']={'OPENMX_CPU_REALSPACE_D3':'export'}
        p=wf.plan_config(self.config)
        self.assertFalse(p['nsoc_density_projection_explicit'])
        self.assertTrue(p['nsoc_realspace_HS_projection_explicit'])
        self.assertEqual(p['nsoc_export_HS'],'PROJECTED_REALSPACE_H_S')
        self.assertEqual(p['nsoc_realspace_HS_projection_scope'],'export_only')
        self.assertFalse(p['nsoc_sector_solver_all_steps'])

    def test_input_environment_symlinks_freeze_resolved_target(self):
        old=self.root/'basis-old';new=self.root/'basis-new';link=self.root/'basis-link'
        old.write_bytes(b'old');new.write_bytes(b'new');link.symlink_to(old)
        self.config['profiles']['nsoc']['environment']={'OPENMX_CPU_NSOC_SYM':'1','OPENMX_CPU_NSOC_BASIS':str(link),'OPENMX_CPU_AUDIT_INPUT_FILE':str(link)}
        p=wf.plan_config(self.config)
        link.unlink();link.symlink_to(new)
        e=p['profiles']['nsoc']['environment']
        self.assertEqual(e['OPENMX_CPU_NSOC_BASIS'],str(old))
        self.assertEqual(e['OPENMX_CPU_AUDIT_INPUT_FILE'],str(old))
        self.assertEqual(Path(e['OPENMX_CPU_NSOC_BASIS']).read_bytes(),b'old')
        wf.verify_dependencies(p)

    def test_sector_parallel_native_layout_scope_enforced(self):
        self.config['runtime']['transport']='shm:ofi'
        self.config['layout']={'hosts':['node040','node041','node042','node043'],'ranks_per_host':16,'threads':1,'cpu_ids_per_host':{h:list(range(16))for h in ['node040','node041','node042','node043']}}
        self.config['profiles']['nsoc']['environment']={'OPENMX_CPU_NSOC_SECTOR_PARALLEL':'1','OPENMX_CPU_NSOC_SECTOR_RANKS':'48,8,8'}
        with self.assertRaisesRegex(ValueError,'sector|SECTOR|layout'):
            wf.plan_config(self.config)
        self.config['layout']['threads']=4
        self.config['layout']['cpu_ids_per_host']={h:list(range(64))for h in self.config['layout']['hosts']}
        self.assertEqual(wf.plan_config(self.config)['threads'],4)
        self.config['layout']['ranks_per_host']=64;self.config['layout']['threads']=1
        self.config['profiles']['nsoc']['environment']['OPENMX_CPU_NSOC_SECTOR_RANKS']='192,32,32'
        self.assertEqual(wf.plan_config(self.config)['ranks'],256)

    def test_numa_grid_rejected_for_threaded_rank_affinity_in_either_stage(self):
        self.config['layout']['threads']=4
        self.config['layout']['cpu_ids_per_host']['node040']=list(range(8))
        for stage in ('nsoc','soc'):
            self.config['profiles'][stage]['environment']={'OPENMX_CPU_GRID_NUMA':'1'}
            with self.assertRaisesRegex(ValueError,'GRID_NUMA.*affinity|GRID_NUMA.*thread'):
                wf.plan_config(self.config)
            self.config['profiles'][stage]['environment']={'OPENMX_CPU_GRID_NUMA':'0'}
        self.assertEqual(wf.plan_config(self.config)['threads'],4)

    def test_inherited_gomp_affinity_is_removed(self):
        self.assertNotIn('GOMP_CPU_AFFINITY',wf.clean_environment({'GOMP_CPU_AFFINITY':'0-63'},{}))

    def test_sector_parallel_requires_native_cpu_id_range_and_three_partitions(self):
        hosts=['node040','node041','node042','node043']
        self.config['runtime']['transport']='shm:ofi'
        self.config['layout']={'hosts':hosts,'ranks_per_host':16,'threads':4,'cpu_ids_per_host':{h:list(range(64,128))for h in hosts}}
        env={'OPENMX_CPU_NSOC_SECTOR_PARALLEL':'1','OPENMX_CPU_NSOC_SECTOR_RANKS':'48,8,8'}
        self.config['profiles']['nsoc']['environment']=env
        with self.assertRaisesRegex(ValueError,'CPU|cpu'):
            wf.plan_config(self.config)
        self.config['layout']['cpu_ids_per_host']={h:list(range(64))for h in hosts}
        env['OPENMX_CPU_NSOC_SECTOR_RANKS']='32,32'
        with self.assertRaisesRegex(ValueError,'three|3|partition'):
            wf.plan_config(self.config)
        env['OPENMX_CPU_NSOC_SECTOR_RANKS']='48,8,8'
        self.assertEqual(wf.plan_config(self.config)['ranks'],64)

    def test_elpa_thread_override_requires_sector64mpi_four_threads(self):
        env={'OPENMX_CPU_NSOC_ELPA_OMP_THREADS':'1'}
        self.config['profiles']['nsoc']['environment']=env
        with self.assertRaisesRegex(ValueError,'ELPA|sector'):
            wf.plan_config(self.config)
        hosts=['node040','node041','node042','node043']
        self.config['runtime']['transport']='shm:ofi'
        self.config['layout']={'hosts':hosts,'ranks_per_host':16,'threads':4,'cpu_ids_per_host':{h:list(range(64))for h in hosts}}
        env.update(OPENMX_CPU_NSOC_SECTOR_PARALLEL='1',OPENMX_CPU_NSOC_SECTOR_RANKS='48,8,8')
        self.assertEqual(wf.plan_config(self.config)['threads'],4)
        self.config['layout'].update(ranks_per_host=64,threads=1)
        env['OPENMX_CPU_NSOC_SECTOR_RANKS']='224,16,16'
        with self.assertRaisesRegex(ValueError,'ELPA|sector'):
            wf.plan_config(self.config)

    def test_effective_default_charge_matches(self):
        self.s.write_text(self.stext+'\nscf.system.charge 0.0\n')
        self.assertEqual(wf.plan_config(self.config)['ranks'],2)

    def test_different_effective_charge_rejected(self):
        self.s.write_text(self.stext+'\nscf.system.charge 1\n')
        with self.assertRaisesRegex(ValueError,'charge'):
            wf.plan_config(self.config)

    def test_structure_and_basis_must_match(self):
        self.s.write_text(self.stext.replace('H5.0-s1','H5.0-s2'))
        with self.assertRaisesRegex(ValueError,'structure|basis'):
            wf.plan_config(self.config)

    def test_output_symlink_into_source_rejected(self):
        (self.root/'link').symlink_to(self.n.parent, target_is_directory=True)
        self.config['output']=str(self.root/'link'/'result')
        with self.assertRaisesRegex(ValueError,'outside'):
            wf.plan_config(self.config)

    def test_duplicate_hosts_and_overlapping_binding_rejected(self):
        self.config['layout']['hosts']=['node040','node040']
        with self.assertRaisesRegex(ValueError,'hosts'):
            wf.plan_config(self.config)
        self.config['layout']['hosts']=['node040']; self.config['layout']['cpu_ids_per_host']['node040']=[0,0]
        with self.assertRaisesRegex(ValueError,'distinct'):
            wf.plan_config(self.config)

    def test_thread_layout_requires_explicit_cpu_ownership(self):
        self.config['layout']['threads']=2
        with self.assertRaisesRegex(ValueError,'CPU|cpu'):
            wf.plan_config(self.config)
        self.config['layout']['cpu_ids_per_host']['node040']=[0,1,2,3]
        self.assertEqual(wf.plan_config(self.config)['total_cores'],4)

    def test_conflicting_runtime_profiles_rejected(self):
        self.config['profiles']['soc']['environment']={'OPENMX_CPU_GSOC':'1','OPENMX_CPU_C3':'3'}
        with self.assertRaisesRegex(ValueError,'conflict'):
            wf.plan_config(self.config)

    def test_projection_never_enabled_implicitly(self):
        self.config['profiles']['soc']['environment']={'OPENMX_CPU_GSOC':'1','OPENMX_CPU_GSOC_PROJECT':'1'}
        with self.assertRaisesRegex(ValueError,'project'):
            wf.plan_config(self.config)

    def test_stale_optimization_and_thread_environment_removed(self):
        e = wf.clean_environment({'HOME':'/home/u','PATH':'/bin','OPENMX_CPU_GSOC':'1','I_MPI_PIN':'0','OMP_NUM_THREADS':'128','MKL_DOMAIN_NUM_THREADS':'MKL_ALL=99','PYTHONPATH':'bad'}, {'OMP_NUM_THREADS':'2'})
        self.assertEqual(e['HOME'],'/home/u'); self.assertEqual(e['OMP_NUM_THREADS'],'2')
        for name in ('OPENMX_CPU_GSOC','I_MPI_PIN','MKL_DOMAIN_NUM_THREADS','PYTHONPATH'):
            self.assertNotIn(name,e)

    def test_source_checkpoint_rank_count_and_files_required(self):
        rst=self.root/'scalar_rst'; rst.mkdir()
        (rst/'scalar.crst_check').write_text('2 1 0\n')
        for rank in range(2):
            (rst/f'scalar.crst0_{rank}_0').write_bytes(b'x')
            (rst/f'scalar.adrst{rank}').write_bytes(b'x')
        self.assertEqual(wf.check_checkpoint(rst,'scalar',2),rst)
        with self.assertRaisesRegex(ValueError,'rank'):
            wf.check_checkpoint(rst,'scalar',4)
        (rst/'scalar.crst0_1_0').unlink()
        with self.assertRaisesRegex(ValueError,'restart|checkpoint'):
            wf.check_checkpoint(rst,'scalar',2)

    def test_native_completion_does_not_imply_convergence(self):
        m = wf.parse_native('The calculation was normally finished.\n', 'Utot. -12.5\nElapsed.Time. 1.25\n')
        self.assertTrue(m['normal_finish']); self.assertFalse(m['native_converged'])
        m = wf.parse_native('The calculation was normally finished.\nThe SCF was achieved at MD= 1\nMD= 1 SCF=32\n', 'Utot. -12.5\nElapsed.Time. 1.25\n')
        self.assertTrue(m['native_converged']); self.assertEqual(m['SCF_steps'],32)
        self.assertEqual(m['energy_Ha'],-12.5); self.assertEqual(m['native_wall_seconds'],1.25)

    def test_soc_must_read_own_scalar_checkpoint(self):
        self.s.write_text(self.stext.replace('scf.restart.filename scalar','scf.restart.filename somebodyelse'))
        with self.assertRaisesRegex(ValueError,'checkpoint|restart'):
            wf.plan_config(self.config)


class MemoryContracts(unittest.TestCase):
    def test_cluster_peak_is_simultaneous_not_sum_of_node_peaks(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            for h,vals in [('a',[10,2]),('b',[2,10])]:
                n=d/'nodes'/h;n.mkdir(parents=True)
                (n/'final.json').write_text(json.dumps({'slice_peak_bytes':10}))
                samples=[{'sample_id':i,'time_unix':i+0.001,'slice_current_bytes':v} for i,v in enumerate(vals)]
                (n/'samples.jsonl').write_text('\n'.join(map(json.dumps,samples)))
            m=wf.memory_summary(d,['a','b'],.25)
            self.assertEqual(m['synchronized_cluster_peak_bytes'],12)
            self.assertEqual(m['sum_independent_node_peaks_bytes'],20)

    def test_missing_common_sample_id_does_not_fabricate_peak(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            for h,i in [('a',0),('b',1)]:
                n=d/'nodes'/h;n.mkdir(parents=True)
                (n/'final.json').write_text(json.dumps({'slice_peak_bytes':10}))
                (n/'samples.jsonl').write_text(json.dumps({'sample_id':i,'time_unix':i*.1,'slice_current_bytes':10})+'\n')
            m=wf.memory_summary(d,['a','b'],.25)
            self.assertIsNone(m['synchronized_cluster_peak_bytes'])

class CleanupContracts(unittest.TestCase):
    def test_partial_monitor_start_is_registered_before_memory_limit_failure(self):
        from unittest.mock import patch
        from types import SimpleNamespace
        import subprocess
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            p={'output':str(d/'run'),'hosts':['node040'],'timeout_seconds':10,'environment_python':'/bin/python',
               'monitor_interval_seconds':.25,'memory_limit_bytes_per_host':1024}
            with patch.object(wf,'ssh',side_effect=[SimpleNamespace(stdout='',stderr=''),subprocess.CalledProcessError(1,['systemctl'])]):
                with self.assertRaises(subprocess.CalledProcessError):
                    wf.start_monitors(p,'nsoc',d)
            record=json.loads((d/'monitor_units.json').read_text())
            self.assertIn('node040',record)

    def test_stop_monitors_stops_only_this_task_slice(self):
        from unittest.mock import patch
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as t:
            d=Path(t); n=d/'nodes'/'node040';n.mkdir(parents=True)
            (n/'final.json').write_text('{}')
            p={'output':str(d.parent/'run')}
            slice_name,unit=wf._units(p,'nsoc')
            u={'node040':{'unit':unit,'slice':slice_name}}
            with patch.object(wf,'ssh',return_value=SimpleNamespace(stdout='ActiveState=inactive\n',returncode=0)) as remote:
                wf.stop_monitors(p,'nsoc',d,u)
            stopped=[call.args[1] for call in remote.call_args_list if 'stop' in call.args[1]]
            self.assertEqual(stopped,[['systemctl','--user','stop',slice_name]])
            self.assertTrue((d/'MONITOR_STOP').exists())

class PreparationContracts(unittest.TestCase):
    def test_common_metadata_command_never_enables_projection_by_default(self):
        from types import SimpleNamespace
        args=SimpleNamespace(environment_prefix=sys.prefix,group='C3',symprec_ladder=[1e-6,3e-6],geometry_tol=5e-5,project_nsoc=False)
        cmd=wf.preparation_command(args,Path('/case/nsoc.dat'),Path('/case/soc.dat'),Path('/new/meta'))
        self.assertIn(str(wf.REPO/'cpu/geometry/prepare.py'),cmd)
        self.assertEqual(cmd[cmd.index('--group')+1],'C3')
        self.assertNotIn('--project-nsoc',cmd)
        args.project_nsoc=True
        self.assertIn('--project-nsoc',wf.preparation_command(args,Path('/n'),Path('/s'),Path('/o')))

class ExportPolicyContracts(unittest.TestCase):
    def profile(self,**env): return {'name':'explicit','environment':env}
    def test_legacy_native_export_requires_explicit_symmetrize_policy(self):
        env={'OPENMX_CPU_C3':'3','OPENMX_CPU_D3':'1','OPENMX_CPU_C3_BASIS':'/x/basis','OPENMX_CPU_SOC_NATIVE_D3':'1','OPENMX_CPU_SOC_NATIVE_MAP':'/x/map'}
        with self.assertRaisesRegex(ValueError,'export|symmetrize'):
            wf.validate_profile('soc',self.profile(**env),{})
        self.assertEqual(wf.validate_profile('soc',self.profile(**env),{'soc_export_policy':'symmetrize'}),env)
    def test_native_group_selector_also_requires_export_policy(self):
        env={'OPENMX_CPU_C3':'3','OPENMX_CPU_D3':'0','OPENMX_CPU_C3_BASIS':'/x/basis','OPENMX_CPU_SOC_NATIVE_GROUP':'C3','OPENMX_CPU_SOC_NATIVE_MAP':'/x/map'}
        with self.assertRaisesRegex(ValueError,'export|symmetrize'):
            wf.validate_profile('soc',self.profile(**env),{'soc_export_policy':'raw'})
    def test_general_soc_does_not_claim_symmetrized_export(self):
        with self.assertRaisesRegex(ValueError,'generic|GSOC|raw'):
            wf.validate_profile('soc',self.profile(OPENMX_CPU_GSOC='1',OPENMX_CPU_GSOC_BASIS='/x/basis'),{'soc_export_policy':'symmetrize'})
    def test_legacy_c3_cannot_enable_d3_time_reversal(self):
        env={'OPENMX_CPU_C3':'3','OPENMX_CPU_D3':'0','OPENMX_CPU_C3_BASIS':'/x/basis','OPENMX_CPU_SOC_NATIVE_GROUP':'C3','OPENMX_CPU_SOC_NATIVE_MAP':'/x/map','OPENMX_CPU_SOC_TR':'1'}
        with self.assertRaisesRegex(ValueError,'D3|C3|native'):
            wf.validate_profile('soc',self.profile(**env),{'soc_export_policy':'symmetrize'})
    def test_d3_time_reversal_requires_native_map(self):
        env={'OPENMX_CPU_C3':'3','OPENMX_CPU_D3':'1','OPENMX_CPU_C3_BASIS':'/x/basis','OPENMX_CPU_SOC_TR':'1','OPENMX_CPU_SOC_REAL':'1'}
        with self.assertRaisesRegex(ValueError,'D3|native|map'):
            wf.validate_profile('soc',self.profile(**env),{})

class ForegroundLifecycleContracts(unittest.TestCase):
    def test_async_submit_rejected_before_output_mutation(self):
        from unittest.mock import patch
        with patch.object(wf,'plan_config') as planned:
            with self.assertRaisesRegex(ValueError,'wait|foreground'):
                wf.submit({},wait=False)
            planned.assert_not_called()
    def test_frozen_inputs_revalidated_before_every_stage(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'input.dat';p.write_text('a')
            plan={'dependencies_sha256':{str(p):wf.sha256(p)}}
            wf.verify_dependencies(plan)
            p.write_text('changed')
            with self.assertRaisesRegex(ValueError,'changed'):
                wf.verify_dependencies(plan)

if __name__=='__main__': unittest.main()
