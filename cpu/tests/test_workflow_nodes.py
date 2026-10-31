import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

P=Path(__file__).parents[1]/'tools/workflow_node.py'
S=importlib.util.spec_from_file_location('workflow_node',P)
node=None
if P.exists():
    node=importlib.util.module_from_spec(S); S.loader.exec_module(node)

class NodeContracts(unittest.TestCase):
    def setUp(self): self.assertIsNotNone(node,'node monitor/rank implementation missing')
    def test_local_rank_uses_intel_local_rank_first(self):
        self.assertEqual(node.rank_identity({'PMI_RANK':'10','MPI_LOCALRANKID':'2','PMI_SIZE':'16'}),(10,2))
        with self.assertRaisesRegex(ValueError,'MPI'):
            node.rank_identity({})
    def test_affinity_is_per_rank_thread_chunk(self):
        plan={'threads':2,'ranks_per_host':2,'layout':{'cpu_ids_per_host':{'node040':[4,5,8,9]}}}
        self.assertEqual(node.rank_cpus(plan,'node040',1),[8,9])
        with self.assertRaisesRegex(ValueError,'rank'):
            node.rank_cpus(plan,'node040',2)
    def test_physical_core_topology_rejects_duplicate_threads(self):
        with tempfile.TemporaryDirectory() as t:
            root=Path(t)
            for cpu in [0,64]:
                p=root/('cpu'+str(cpu))/'topology';p.mkdir(parents=True)
                (p/'physical_package_id').write_text('0')
                (p/'core_id').write_text('0')
            with self.assertRaisesRegex(ValueError,'physical'):
                node.cpu_topology([0,64],root)

    def test_cgroup_path_only_below_user_slice(self):
        with self.assertRaisesRegex(ValueError,'slice'):
            node.task_slice_path('0::/system.slice/foo.service\n')
        got=node.task_slice_path('0::/user.slice/user-1021.slice/user@1021.service/omxwfx.slice/a.service\n')
        self.assertTrue(str(got).endswith('omxwfx.slice'))


class DescriptorPreservation(unittest.TestCase):
    def test_pmi_descriptor_survives_wrapper_process(self):
        self.assertIsNotNone(node)
        import os, subprocess, sys
        read,write=os.pipe()
        self.addCleanup(os.close,read);self.addCleanup(os.close,write)
        os.set_inheritable(write,True)
        command=[sys.executable,'-c','import os,sys;os.write(int(sys.argv[1]),b"pmi")',str(write)]
        code=node.run_scope_command(command)
        self.assertEqual(code,0)
        self.assertEqual(os.read(read,3),b'pmi')

class ThreadBindingContracts(unittest.TestCase):
    def test_thread_places_follow_explicit_rank_cpus_and_remove_conflicts(self):
        env=node.binding_environment({'KMP_AFFINITY':'disabled','GOMP_CPU_AFFINITY':'0-63','OMP_PROC_BIND':'false','OMP_PLACES':'cores'},[4,6,8,10])
        self.assertNotIn('KMP_AFFINITY',env)
        self.assertNotIn('GOMP_CPU_AFFINITY',env)
        self.assertEqual(env['OMP_PROC_BIND'],'close')
        self.assertEqual(env['OMP_PLACES'],'{4},{6},{8},{10}')

if __name__=='__main__': unittest.main()
