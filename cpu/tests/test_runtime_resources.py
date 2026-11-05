"""A built CPU executable carries its automatic-symmetry preparation resources."""
from pathlib import Path
import importlib.util
import tempfile
import unittest

CPU=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_resource_builder',CPU/'build.py')
builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)


class RuntimeResources(unittest.TestCase):
    def test_build_plan_carries_python_and_helper_paths(self):
        with tempfile.TemporaryDirectory()as d:
            output=Path(d)/'output'
            import json
            selection=json.loads((CPU/'overrides/nsoc/UNIFIED_SELECTION.json').read_text())
            plan=builder.source_plan(selection,profile='unified')
            commands,_=builder.command_plan(plan,output,Path('/opt/intel/oneapi'),'CORE-AVX2')
            args=next(x['command']for x in commands if x['label']=='Input_std.o')
            self.assertTrue(any('OPENMX_SYMMETRY_BOOTSTRAP='in x for x in args))
            self.assertTrue(any('OPENMX_SYMMETRY_PYTHON='in x for x in args))

    def test_resource_bundle_is_self_contained_without_native_source_or_personal_tree(self):
        self.assertTrue(hasattr(builder,'copy_runtime_resources'),'missing relocatable helper bundle')
        with tempfile.TemporaryDirectory()as d:
            root=builder.copy_runtime_resources(Path(d))
            self.assertTrue((root/'cpu/geometry/bootstrap.py').is_file())
            self.assertTrue((root/'cpu/geometry/metadata/symmetry_detect.py').is_file())
            self.assertTrue((root/'cpu/overrides/nsoc/metadata/prepare_paired_d3.py').is_file())
            self.assertTrue((root/'cpu/tools/workflow.py').is_file())
            self.assertTrue((root/'cpu/tools/examples/profile_templates.json').is_file())
            self.assertFalse((root/'cpu/vendor').exists())


if __name__=='__main__':unittest.main()
