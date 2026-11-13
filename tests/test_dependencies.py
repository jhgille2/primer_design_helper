import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from primer_helper import core, dependencies as deps

class Dependencies(unittest.TestCase):
    def test_missing_tool(self):
        with patch('primer_helper.dependencies.shutil.which',return_value=None):
            s=deps.probe('mafft','missing')
        self.assertFalse(s.usable);self.assertIsNone(s.path)
    def test_probe_version_and_primer3_boulder_io(self):
        outputs=[SimpleNamespace(returncode=0,stdout='',stderr='v7.526'),SimpleNamespace(returncode=0,stdout='PRIMER_ERROR=No sequence\n=',stderr='')]
        with patch('primer_helper.dependencies.shutil.which',side_effect=lambda x:x),patch('primer_helper.dependencies.subprocess.run',side_effect=outputs) as run:
            self.assertTrue(deps.probe('mafft','mafft').usable)
            self.assertTrue(deps.probe('primer3_core','primer3_core').usable)
            self.assertEqual(run.call_args.kwargs['input'],'=\n')
    def test_package_plan_deduplicates_blast_suite(self):
        statuses=[deps.ToolStatus(n,n,None,False,'') for n in deps.TOOL_NAMES]
        self.assertEqual(deps.package_plan(statuses),['mafft','primer3','blast'])
    def test_windows_is_not_silently_claimed_supported(self):
        with patch('primer_helper.dependencies.platform.system',return_value='Windows'):
            with self.assertRaises(core.InputError):deps.install_missing([], 'conda', Path('/tmp/nonexistent'))
    def test_install_is_isolated_and_uses_conda_forge_bioconda(self):
        statuses=[deps.ToolStatus('blastn','blastn',None,False,'')]
        with tempfile.TemporaryDirectory() as tmp:
            prefix=Path(tmp)/'isolated'
            with patch('primer_helper.dependencies.platform.system',return_value='Linux'),patch('primer_helper.dependencies.shutil.which',return_value='/usr/bin/conda'),patch('primer_helper.dependencies.subprocess.run',return_value=SimpleNamespace(returncode=0,stdout='',stderr='')) as run,patch('primer_helper.dependencies.probe',return_value=deps.ToolStatus('blastn','x','/installed/blastn',True,'ok')):
                result=deps.install_missing(statuses,'conda',prefix)
            cmd=run.call_args.args[0]
            self.assertEqual(cmd[1],'create');self.assertIn(str(prefix),cmd)
            self.assertIn('conda-forge',cmd);self.assertIn('bioconda',cmd);self.assertEqual(cmd[-1],'blast')
            self.assertEqual(result,{'blastn':'/installed/blastn'})
    def test_refuses_existing_non_conda_folder(self):
        statuses=[deps.ToolStatus('mafft','mafft',None,False,'')]
        with tempfile.TemporaryDirectory() as tmp,patch('primer_helper.dependencies.platform.system',return_value='Linux'),patch('primer_helper.dependencies.shutil.which',return_value='/conda'):
            with self.assertRaises(core.InputError):deps.install_missing(statuses,'conda',Path(tmp))

if __name__=='__main__':unittest.main()
