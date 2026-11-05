"""Calculations and their lifecycle supervision have no assistant time budget."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


contract = module('runtime_contract_fixture', ROOT/'tests/test_workflow_contract.py')
wf = contract.wf
node = module('runtime_node', ROOT/'tools/workflow_node.py')


class UnlimitedRuntime(unittest.TestCase):
    def fixture(self):
        fixture = contract.WorkflowContract('test_plan_is_read_only_and_preserves_inputs')
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        return fixture

    def test_plan_has_no_default_daily_limit(self):
        fixture = self.fixture()
        self.assertIsNone(wf.plan_config(fixture.config)['timeout_seconds'])

    def test_legacy_finite_limits_are_rejected_not_silently_dropped(self):
        fixture = self.fixture()
        for value in (1, 86400, 0, -1, True, '86400'):
            with self.subTest(value=value):
                fixture.config['timeout_seconds'] = value
                with self.assertRaisesRegex(ValueError, 'timeout_seconds.*finite|finite.*timeout_seconds'):
                    wf.plan_config(fixture.config)

    def test_explicit_unlimited_config(self):
        fixture = self.fixture()
        for value in (None, 'infinity'):
            with self.subTest(value=value):
                fixture.config['timeout_seconds'] = value
                self.assertIsNone(wf.plan_config(fixture.config)['timeout_seconds'])

    def test_ssh_live_remote_work_has_no_process_deadline(self):
        with patch.object(wf.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
            wf.ssh('node040', ['long-controller'])
        self.assertIsNone(run.call_args.kwargs.get('timeout'))
        self.assertIn('ConnectTimeout=8', run.call_args.args[0])

    def test_controller_service_and_foreground_wait_are_unlimited(self):
        fixture = self.fixture()
        plan = wf.plan_config(fixture.config)
        plan['build_manifest'] = '/dummy/manifest'
        with patch.object(wf, 'plan_config', return_value=plan), \
             patch.object(wf, 'ssh', return_value=SimpleNamespace(returncode=0, stdout='', stderr='')) as ssh, \
             patch.object(wf, 'inspect', return_value={'state':'complete'}):
            wf.submit(fixture.config, wait=True)
        command = ssh.call_args.args[1]
        self.assertIn('--property=RuntimeMaxSec=infinity', command)
        self.assertIn('--property=TimeoutStartSec=infinity', command)
        self.assertIn('--property=TimeoutStopSec=infinity', command)
        self.assertIsNone(ssh.call_args.kwargs.get('timeout'))

    def test_slow_monitor_startup_is_not_a_calculation_time_limit(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            plan = {'output':str(directory), 'hosts':['node040'], 'timeout_seconds':None,
                    'environment_python':'/bin/python', 'monitor_interval_seconds':.25}
            calls = []
            def remote(host, command, **kwargs):
                calls.append(command)
                return SimpleNamespace(returncode=0, stdout='ActiveState=active\nSubState=running\n', stderr='')
            sleeps = 0
            def progress(_):
                nonlocal sleeps
                sleeps += 1
                if sleeps == 3:
                    ready = directory/'nodes/node040/ready.json'
                    ready.parent.mkdir(parents=True)
                    ready.write_text('{}')
            with patch.object(wf, 'ssh', side_effect=remote), \
                 patch.object(wf.time, 'sleep', side_effect=progress), \
                 patch.object(wf.time, 'monotonic', side_effect=[0, 1e9]):
                units = wf.start_monitors(plan, 'nsoc', directory)
            self.assertIn('node040', units)
            self.assertIn('--property=RuntimeMaxSec=infinity', calls[0])
            self.assertIn('--property=TimeoutStartSec=infinity', calls[0])
            self.assertIn('--property=TimeoutStopSec=infinity', calls[0])

    def test_finite_saved_plan_cannot_bypass_unlimited_config_validation(self):
        with patch.object(wf, 'verify_dependencies') as verify:
            with self.assertRaisesRegex(ValueError, 'timeout_seconds.*finite'):
                wf.run_stage({'timeout_seconds':86400}, 'nsoc')
        verify.assert_not_called()

    def test_native_mpi_run_has_no_process_deadline(self):
        fixture = self.fixture()
        plan = wf.plan_config(fixture.config)
        Path(plan['output']).mkdir()
        with patch.object(wf, 'verify_dependencies'), \
             patch.object(wf, 'start_monitors', return_value={}), \
             patch.object(wf, '_library_environment', return_value={}), \
             patch.object(wf, '_verify_rank_records', return_value=[]), \
             patch.object(wf.subprocess, 'run', return_value=SimpleNamespace(returncode=1)) as run:
            wf.run_stage(plan, 'nsoc')
        self.assertIsNone(run.call_args.kwargs['timeout'])

    def test_slow_anchor_startup_is_not_cancelled(self):
        with tempfile.TemporaryDirectory() as temp:
            plan = {'output':temp, 'hosts':['node040'], 'environment_python':'/bin/python'}
            process = Mock()
            process.poll.return_value = None
            sleeps = 0
            def progress(_):
                nonlocal sleeps
                sleeps += 1
                if sleeps == 3:
                    (Path(temp)/'anchors/node040.ready.json').write_text('{}')
            with patch.object(wf.subprocess, 'Popen', return_value=process), \
                 patch.object(wf.time, 'sleep', side_effect=progress), \
                 patch.object(wf.time, 'monotonic', side_effect=[0, 1e9]):
                anchors = wf.start_anchors(plan)
            self.assertEqual(set(anchors), {'node040'})
            wf.stop_anchors(plan, anchors)

    def test_slow_monitor_final_accounting_does_not_stop_it_early(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            plan = {'output':temp}
            unit = {'node040':{'unit':'own-monitor.service', 'slice':'own.slice'}}
            sleeps = 0
            def progress(_):
                nonlocal sleeps
                sleeps += 1
                if sleeps == 3:
                    final = directory/'nodes/node040/final.json'
                    final.parent.mkdir(parents=True)
                    final.write_text('{}')
            def remote(host, command, **kwargs):
                if 'stop' in command:
                    self.assertTrue((directory/'nodes/node040/final.json').exists())
                return SimpleNamespace(returncode=0, stdout='ActiveState=active\n', stderr='')
            with patch.object(wf, 'ssh', side_effect=remote), \
                 patch.object(wf.time, 'sleep', side_effect=progress), \
                 patch.object(wf.time, 'monotonic', side_effect=[0, 1e9]):
                wf.stop_monitors(plan, 'nsoc', directory, unit)

    def test_failed_monitor_cleans_own_slice_and_reports_missing_accounting(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            units = {'node040':{'unit':'own-monitor.service', 'slice':'own.slice'}}
            with patch.object(wf, 'ssh', return_value=SimpleNamespace(returncode=0, stdout='ActiveState=failed\n', stderr='')) as remote:
                with self.assertRaisesRegex(RuntimeError, 'monitor.*(failed|ended)'):
                    wf.stop_monitors({'output':temp}, 'nsoc', directory, units)
            self.assertIn(['systemctl','--user','stop','own.slice'], [call.args[1] for call in remote.call_args_list])

    def test_failed_monitor_still_stops_startup_promptly(self):
        with tempfile.TemporaryDirectory() as temp:
            directory = Path(temp)
            plan = {'output':str(directory), 'hosts':['node040'], 'timeout_seconds':None,
                    'environment_python':'/bin/python', 'monitor_interval_seconds':.25}
            answers = [SimpleNamespace(returncode=0, stdout='', stderr=''),
                       SimpleNamespace(returncode=0, stdout='ActiveState=failed\n', stderr='')]
            with patch.object(wf, 'ssh', side_effect=answers):
                with self.assertRaisesRegex(RuntimeError, 'monitor.*(failed|ended)'):
                    wf.start_monitors(plan, 'nsoc', directory)

    def test_rank_scope_is_explicitly_unlimited(self):
        with tempfile.TemporaryDirectory() as temp:
            plan_path = Path(temp)/'plan.json'
            plan_path.write_text(json.dumps({'output':temp, 'environment_python':'/bin/python'}))
            with patch.object(node, 'rank_identity', return_value=(0, 0)), \
                 patch.object(node, 'run_scope_command', return_value=0) as command:
                node.rank(SimpleNamespace(plan=plan_path, stage='nsoc'))
            self.assertIn('--property=RuntimeMaxSec=infinity', command.call_args.args[0])
            self.assertIn('--property=TimeoutStopSec=infinity', command.call_args.args[0])

    def test_anchor_cleanup_waits_without_timer_kills(self):
        with tempfile.TemporaryDirectory() as temp:
            (Path(temp)/'anchors').mkdir()
            process = Mock()
            wf.stop_anchors({'output':temp}, {'node040':(process, Mock())})
            process.wait.assert_called_once_with(timeout=None)
            process.terminate.assert_not_called()
            process.kill.assert_not_called()


if __name__ == '__main__':
    unittest.main()
