import gzip
import hashlib
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from primer_helper import core, soybean

class Response(io.BytesIO):
    def __init__(self,data):
        super().__init__(data); self.headers={'Content-Length':str(len(data))}

class LibraryTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.ref=soybean.BY_KEY['w05_v1']
        self.gz=gzip.compress(b'>chr1\nACGTACGTACGT\n')
        self.md5=hashlib.md5(self.gz).hexdigest()
        self.manifest=f'{self.md5}  ./{self.ref.filename}\n'.encode()
    def tearDown(self):self.temp.cleanup()
    def open_url(self,request,timeout=None):
        return Response(self.manifest if request.full_url.endswith('.md5') else self.gz)
    def install_mock(self):
        with patch('primer_helper.soybean.urllib.request.urlopen',side_effect=self.open_url):
            return soybean.download_and_prepare(self.ref,self.root)
    def fake_makeblastdb(self,args,**kwargs):
        if '-out' in args:
            prefix=args[args.index('-out')+1]
            for ext in ('.nhr','.nin','.nsq'):Path(prefix+ext).write_bytes(b'index')
            return ''
        if args[0]=='blastn':
            return 'chr1\t99.0\t4\t1\t4\t10\t7\t12\t1e-9\t42\n'
        if args[0]=='blastdbcmd':return 'ACGTACGT\n'
        raise AssertionError(args)
    def install_database(self):
        self.install_mock()
        with patch('primer_helper.soybean._run',side_effect=self.fake_makeblastdb):
            return soybean.build_database(self.ref,self.root)
    def test_checksum_formats_and_missing(self):
        self.assertEqual(soybean.parse_checksum(self.manifest.decode(),self.ref.filename),self.md5)
        with self.assertRaises(core.InputError):soybean.parse_checksum('not a checksum','missing')
    def test_download_extract_checksum_and_metadata(self):
        path=self.install_mock();self.assertEqual(path.read_bytes(),b'>chr1\nACGTACGTACGT\n')
        info=json.loads((path.parent/'reference.json').read_text())
        self.assertEqual(info['source_md5'],self.md5);self.assertEqual(info['reference']['assembly'],'W05.gnm1')
    def test_bad_checksum_preserves_existing_file(self):
        folder=soybean.reference_dir(self.root,self.ref);folder.mkdir(parents=True)
        old=folder/self.ref.filename;old.write_bytes(b'existing')
        def bad(req,timeout=None):return Response(self.manifest if req.full_url.endswith('.md5') else b'corrupt')
        with patch('primer_helper.soybean.urllib.request.urlopen',side_effect=bad):
            with self.assertRaises(core.InputError):soybean.download_and_prepare(self.ref,self.root)
        self.assertEqual(old.read_bytes(),b'existing');self.assertFalse((folder/(self.ref.filename+'.part')).exists())
    def test_build_and_find_minus_hit_and_extract_flanking_sequence(self):
        db=self.install_database()
        with patch('primer_helper.soybean._run',side_effect=self.fake_makeblastdb) as run:
            hits=soybean.find_hits(self.root,[self.ref],'A'*30)
            self.assertEqual(hits[0].strand,'minus');self.assertEqual(hits[0].subject_length,12)
            name,seq=soybean.extract_hit(self.root,hits[0],2)
            cmd=run.call_args.args[0]
            self.assertEqual(cmd[cmd.index('-range')+1],'5-12')
            self.assertEqual(cmd[cmd.index('-strand')+1],'minus')
            self.assertIn('w05_v1|chr1:5-12|strand=minus', name)
            self.assertEqual(seq,'ACGTACGT')
    def test_extract_converts_ambiguity_without_choosing_a_base(self):
        self.install_database();hit=soybean.Hit(self.ref.key,'chr1',100,4,1,4,1,4,12,0.,40.)
        with patch('primer_helper.soybean._run',return_value='ACRG'):
            name,seq=soybean.extract_hit(self.root,hit,0)
        self.assertEqual(seq,'ACNG');self.assertIn('ambiguous_to_N=1',name)
    def test_unready_or_stale_index_not_used(self):
        with self.assertRaises(core.InputError):soybean.database_prefix(self.root,self.ref)
        self.install_database();(soybean.reference_dir(self.root,self.ref)/'blast'/'genome.nin').unlink()
        with self.assertRaises(core.InputError):soybean.database_prefix(self.root,self.ref)
    @unittest.skipIf(os.name == "nt", "Unix executable fixture uses a shebang")
    def test_subprocess_pipeline_with_executables_in_path_with_spaces(self):
        import sys
        self.install_mock()
        bindir=self.root/'mock tools';bindir.mkdir()
        script=(f'#!{sys.executable}\n'
                'import pathlib,sys\n'
                'name=pathlib.Path(sys.argv[0]).name\n'
                'if name=="makeblastdb":\n'
                ' p=sys.argv[sys.argv.index("-out")+1]\n'
                ' for ext in (".nhr",".nin",".nsq"): pathlib.Path(p+ext).write_text("db")\n'
                'elif name=="blastn":\n'
                ' sys.stdin.read()\n'
                ' print("chr1\\t100\\t4\\t1\\t4\\t10\\t7\\t12\\t1e-9\\t42")\n'
                'else:\n'
                ' print("ACGTACGT")\n')
        for name in ('makeblastdb','blastn','blastdbcmd'):
            path=bindir/name;path.write_text(script);path.chmod(0o755)
        soybean.build_database(self.ref,self.root,str(bindir/'makeblastdb'))
        hits=soybean.find_hits(self.root,[self.ref],'A'*30,str(bindir/'blastn'))
        name,sequence=soybean.extract_hit(self.root,hits[0],2,str(bindir/'blastdbcmd'))
        self.assertIn('chr1:5-12',name);self.assertEqual(sequence,'ACGTACGT')
    def test_query_validation(self):
        self.assertEqual(soybean.normalize_query('>x\n'+('A'*15)+'--'+('N'*15)), 'A'*15+'N'*15)
        with self.assertRaises(core.InputError):soybean.normalize_query('>a\n'+'A'*30+'\n>b\n'+'C'*30)

if __name__=='__main__':unittest.main()
