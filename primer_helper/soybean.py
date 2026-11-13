"""SoyBase/LIS whole-genome reference library, local setup, search, and extraction.

No genomes are bundled. All coordinates in BLAST hits and extraction requests
are 1-based inclusive; gap/strand details follow NCBI BLAST's tabular output.
"""
from __future__ import annotations
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable
from .core import InputError

ROOT = 'https://data.legumeinfo.org/Glycine'
# Curated exact paths from SoyBase/LIS collection listings. Wm82 v2 and v6
# are versions of the same accession, not independent diversity samples.
@dataclass(frozen=True)
class Reference:
    key: str
    species: str
    accession: str
    assembly: str
    collection: str
    filename: str
    approx_mb: int
    citation: str
    @property
    def url(self):
        return f'{ROOT}/{self.species}/genomes/{self.collection}/{self.filename}'
    @property
    def checksum_url(self):
        return f'{ROOT}/{self.species}/genomes/{self.collection}/CHECKSUM.{self.collection}.md5'
    @property
    def collection_url(self):
        return f'{ROOT}/{self.species}/genomes/{self.collection}/'

REFERENCES = (
    Reference('wm82_v6','max','Williams 82','Wm82.gnm6','Wm82.gnm6.S97D','glyma.Wm82.gnm6.S97D.genome_main.fna.gz',255,'https://data.legumeinfo.org/Glycine/max/genomes/Wm82.gnm6.S97D/'),
    Reference('wm82_v2','max','Williams 82','Wm82.gnm2','Wm82.gnm2.DTC4','glyma.Wm82.gnm2.DTC4.genome_main.fna.gz',287,'https://data.legumeinfo.org/Glycine/max/genomes/Wm82.gnm2.DTC4/'),
    Reference('zh13_v1','max','Zhonghuang 13','Zh13.gnm1','Zh13.gnm1.N6C8','glyma.Zh13.gnm1.N6C8.genome_main.fna.gz',294,'https://data.legumeinfo.org/Glycine/max/genomes/Zh13.gnm1.N6C8/'),
    Reference('w05_v1','soja','W05','W05.gnm1','W05.gnm1.SVL1','glyso.W05.gnm1.SVL1.genome_main.fna.gz',296,'https://data.legumeinfo.org/Glycine/soja/genomes/W05.gnm1.SVL1/'),
    Reference('pi483463_v1','soja','PI 483463','PI483463.gnm1','PI483463.gnm1.YJWS','glyso.PI483463.gnm1.YJWS.genome_main.fna.gz',292,'https://data.legumeinfo.org/Glycine/soja/genomes/PI483463.gnm1.YJWS/'),
    Reference('soyw01_v1','soja','PI 562565 (SoyW01)','PI_562565.gnm1','PI_562565.gnm1.Z2VG','glyso.PI_562565.gnm1.Z2VG.genome_main.fna.gz',261,'https://data.legumeinfo.org/Glycine/soja/genomes/PI_562565.gnm1.Z2VG/'),
)
BY_KEY = {item.key: item for item in REFERENCES}
DEFAULT_ROOT = Path(__file__).resolve().parent.parent / 'data' / 'soybean'

# POWR1 / Glyma.20G085100 in Wm82.a2.v1 (the Wm82.gnm2 assembly).
# Coordinates are 1-based inclusive, with 500 bases added on each side.
POWR1_START = 31774770 - 500
POWR1_END = 31779804 + 500

def _log(callback: Callable[[str], None] | None, msg: str):
    if callback: callback(msg)

def reference_dir(root: Path, ref: Reference) -> Path:
    return Path(root) / ref.key

def powr1_query(root: Path):
    """Extract the published Wm82.gnm2 locus from the installed genome FASTA.

    A chromosome FASTA header must contain a Gm20 token. Refuse missing or
    ambiguous identifiers rather than silently extracting from the wrong contig.
    """
    ref=BY_KEY['wm82_v2'];folder=reference_dir(root,ref)
    try:
        info=json.loads((folder/'reference.json').read_text(encoding='utf-8'))
        if info.get('reference',{}).get('collection') != ref.collection:
            raise InputError('The installed Williams 82 assembly metadata does not match Wm82.gnm2.DTC4.')
        path=folder / ref.filename[:-3]
        pieces=[];position=0;matches=0;collect=False;header=''
        with path.open(encoding='ascii') as stream:
            for line in stream:
                if line.startswith('>'):
                    name=line[1:].split()[0]
                    collect=bool(re.search(r'(?:^|[._|])Gm20(?:$|[._|])',name,re.IGNORECASE))
                    if collect:matches+=1;header=name;position=0
                    continue
                if not collect:continue
                bases=''.join(line.split()).upper();next_position=position+len(bases)
                if position < POWR1_END and next_position >= POWR1_START:
                    pieces.append(bases[max(0,POWR1_START-1-position):min(len(bases),POWR1_END-position)])
                position=next_position
    except (OSError,ValueError,UnicodeError,IndexError) as exc:
        raise InputError('Download Williams 82 Wm82.gnm2 in Soybean library before loading POWR1; its FASTA and metadata must be available.') from exc
    if matches != 1:
        raise InputError(f'Expected exactly one Wm82.gnm2 Gm20 chromosome in the FASTA; found {matches}. Check the assembly and FASTA headers.')
    seq=''.join(pieces)
    if len(seq)!=POWR1_END-POWR1_START+1 or not set(seq)<=set('ACGTNRYKMSWBDHV'):
        raise InputError('POWR1 extraction failed: the chromosome is too short or contains invalid bases.')
    ambiguous=sum(c not in 'ACGTN' for c in seq)
    seq=''.join(c if c in 'ACGTN' else 'N' for c in seq)
    return (f'Wm82.gnm2|{header}:{POWR1_START}-{POWR1_END}|strand=plus|ambiguous_to_N={ambiguous}',seq)

def parse_checksum(text: str, filename: str) -> str:
    for line in text.splitlines():
        m = re.match(r'^([0-9a-fA-F]{32})\s+[* ]?(?:\./)?(.+?)\s*$',line)
        if m and m.group(2) == filename:
            return m.group(1).lower()
    raise InputError(f'Official checksum file has no entry for {filename}.')

def _read_checksum(url: str, filename: str) -> str:
    try:
        req = urllib.request.Request(url,headers={'User-Agent':'PrimerDesignHelper/1.0'})
        with urllib.request.urlopen(req,timeout=60) as response:
            text = response.read(1024 * 1024).decode('utf-8')
    except (OSError,UnicodeError,urllib.error.URLError) as exc:
        raise InputError(f'Could not obtain the published checksum from {url}: {exc}') from exc
    return parse_checksum(text,filename)

def verify_md5(path: Path, expected: str) -> bool:
    digest=hashlib.md5()  # Published integrity checksum, not an authenticity signature.
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(4*1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()==expected

def download_and_prepare(ref: Reference, root: Path, progress=None) -> Path:
    """Stream verified .gz and uncompressed FASTA to disk, leaving no partial output.

    The published MD5 is fetched from the same provider; protects against
    accidental corruption, but is not an independent authenticity guarantee.
    """
    folder=reference_dir(root,ref); folder.mkdir(parents=True,exist_ok=True)
    gz_path=folder/ref.filename
    fasta_path=folder/(ref.filename[:-3])
    metadata=folder/'reference.json'
    _log(progress, f'{ref.accession} {ref.assembly}: checking published checksum')
    expected=_read_checksum(ref.checksum_url,ref.filename)
    if gz_path.exists() and verify_md5(gz_path,expected):
        _log(progress,f'{ref.accession}: existing compressed file verified')
    else:
        # Do not replace a possibly valid older file until the new one checks out.
        part=folder/(ref.filename+'.part')
        try:
            request=urllib.request.Request(ref.url,headers={'User-Agent':'PrimerDesignHelper/1.0'})
            with urllib.request.urlopen(request,timeout=90) as response, part.open('wb') as out:
                expected_bytes=response.headers.get('Content-Length')
                downloaded=0; next_notice=0
                while True:
                    chunk=response.read(1024*1024)
                    if not chunk:break
                    out.write(chunk);downloaded+=len(chunk)
                    if downloaded>=next_notice:
                        _log(progress,f'{ref.accession}: downloaded {downloaded // (1024*1024)} MiB')
                        next_notice=downloaded+25*1024*1024
                if expected_bytes and downloaded!=int(expected_bytes):
                    raise InputError('Download ended before its advertised content length.')
            if not verify_md5(part,expected):raise InputError(f'{ref.accession}: checksum mismatch; downloaded file discarded.')
            os.replace(part,gz_path)
        except (OSError,urllib.error.URLError,ValueError) as exc:
            raise InputError(f'{ref.accession}: download failed: {exc}') from exc
        finally:
            part.unlink(missing_ok=True)
    # Recreate FASTA if metadata does not establish its exact verified source.
    old={}
    try:old=json.loads(metadata.read_text(encoding='utf-8'))
    except (OSError,ValueError):pass
    if not fasta_path.is_file() or old.get('source_md5') != expected:
        temp=folder/(fasta_path.name+'.part')
        try:
            _log(progress,f'{ref.accession}: decompressing and validating FASTA (requires additional disk space)')
            with gzip.open(gz_path,'rb') as inp,temp.open('wb') as out:
                first=inp.read(1)
                if first != b'>':raise InputError('Downloaded file is not a FASTA file.')
                out.write(first);shutil.copyfileobj(inp,out,1024*1024)
            if temp.stat().st_size<10:raise InputError('FASTA file is empty or incomplete.')
            os.replace(temp,fasta_path)
        except (OSError,EOFError,gzip.BadGzipFile) as exc:
            raise InputError(f'{ref.accession}: decompression failed: {exc}') from exc
        finally:
            temp.unlink(missing_ok=True)
    info={'reference':asdict(ref),'source_url':ref.url,'checksum_url':ref.checksum_url,'source_md5':expected,'fasta_filename':fasta_path.name}
    tempmeta=folder/'reference.json.part';tempmeta.write_text(json.dumps(info,indent=2),encoding='utf-8');os.replace(tempmeta,metadata)
    _log(progress,f'{ref.accession}: downloaded and verified; FASTA ready')
    return fasta_path

def _run(args, input_text='', timeout=1200):
    if not shutil.which(args[0]):raise InputError(f'Cannot find {args[0]!r}; install NCBI BLAST+ or enter its full path.')
    try:
        proc=subprocess.run(args,input=input_text,capture_output=True,text=True,timeout=timeout,check=False)
    except (OSError,subprocess.TimeoutExpired) as exc:raise InputError(f'Could not run {args[0]}: {exc}') from exc
    if proc.returncode:raise InputError(f'{args[0]} failed: {proc.stderr.strip()[:3000]}')
    return proc.stdout

def build_database(ref: Reference, root: Path, executable='makeblastdb', progress=None) -> Path:
    folder=reference_dir(root,ref);meta=folder/'reference.json'
    if not meta.exists():raise InputError(f'{ref.accession}: download and verify the genome first.')
    info=json.loads(meta.read_text(encoding='utf-8'));fasta=folder/info['fasta_filename']
    if not fasta.is_file():raise InputError(f'{ref.accession}: local genome FASTA is missing.')
    target=folder/'blast';marker=target/'ready.json';prefix=target/'genome'
    if marker.exists():
        try:
            old=json.loads(marker.read_text(encoding='utf-8'))
            if old.get('source_md5')==info['source_md5'] and all(Path(str(prefix)+suffix).is_file() for suffix in ('.nhr','.nin','.nsq')):
                _log(progress,f'{ref.accession}: existing BLAST database is ready');return prefix
        except (OSError,ValueError):pass
    staging=Path(tempfile.mkdtemp(prefix='.blast-build-',dir=folder))
    try:
        temp_prefix=staging/'genome'
        _log(progress,f'{ref.accession}: building nucleotide BLAST database')
        # v4 avoids v5 LMDB's unusually large virtual-address reservation.
        _run([executable,'-in',str(fasta),'-dbtype','nucl','-parse_seqids','-blastdb_version','4','-out',str(temp_prefix)],timeout=3600)
        if not all(Path(str(temp_prefix)+suffix).is_file() for suffix in ('.nhr','.nin','.nsq')):
            raise InputError('makeblastdb did not produce the expected database files.')
        (staging/'ready.json').write_text(json.dumps({'source_md5':info['source_md5']}),encoding='utf-8')
        if target.exists():
            raise InputError(f'An incomplete BLAST directory already exists at {target}; move it aside, then retry.')
        os.replace(staging,target)
    finally:
        if staging.exists():shutil.rmtree(staging)
    _log(progress,f'{ref.accession}: database ready')
    return prefix

def database_prefix(root: Path, ref: Reference) -> Path:
    folder=reference_dir(root,ref);prefix=folder/'blast'/'genome';marker=folder/'blast'/'ready.json'
    try:
        m=json.loads(marker.read_text(encoding='utf-8'))
        info=json.loads((folder/'reference.json').read_text(encoding='utf-8'))
    except (OSError,ValueError):raise InputError(f'{ref.accession}: no ready BLAST database; download and build it first.')
    if m.get('source_md5')!=info.get('source_md5') or not all(Path(str(prefix)+s).exists() for s in ('.nhr','.nin','.nsq')):
        raise InputError(f'{ref.accession}: BLAST database is incomplete or stale.')
    return prefix

@dataclass(frozen=True)
class Hit:
    reference_key: str
    subject: str
    identity: float
    length: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    subject_length: int
    evalue: float
    bitscore: float
    @property
    def strand(self):return 'plus' if self.subject_start<=self.subject_end else 'minus'

def normalize_query(text: str) -> str:
    raw=text.strip()
    if raw.startswith('>'):
        lines=raw.splitlines()
        if any(line.startswith('>') for line in lines[1:]):raise InputError('Provide a single query sequence.')
        raw=''.join(lines[1:])
    seq=''.join(raw.split()).upper().replace('-','')
    if len(seq)<30 or len(seq)>50000 or not set(seq)<=set('ACGTN'):
        raise InputError('Query must contain 30-50,000 A/C/G/T/N bases (FASTA or plain sequence).')
    return seq

def find_hits(root: Path, refs: list[Reference], query: str, executable='blastn', max_hits=30, progress=None):
    seq=normalize_query(query)
    if not refs:raise InputError('Select at least one installed reference.')
    if not 1<=max_hits<=1000:raise InputError('Maximum hits must be between 1 and 1000.')
    hits=[]
    for ref in refs:
        db=database_prefix(root,ref);_log(progress,f'Searching {ref.accession} {ref.assembly}...')
        output=_run([executable,'-task','megablast','-query','-','-db',str(db),'-max_target_seqs',str(max_hits),'-max_hsps','3',
                     '-outfmt','6 sseqid pident length qstart qend sstart send slen evalue bitscore'],input_text=f'>query\n{seq}\n',timeout=3600)
        for line in output.splitlines():
            fields=line.split('\t')
            if len(fields)!=10:raise InputError(f'Unexpected BLAST output from {ref.assembly}.')
            try:
                hits.append(Hit(ref.key,fields[0],float(fields[1]),int(fields[2]),int(fields[3]),int(fields[4]),int(fields[5]),int(fields[6]),int(fields[7]),float(fields[8]),float(fields[9])))
            except ValueError as exc:raise InputError('Malformed numeric field in BLAST output.') from exc
    return hits

def extract_hit(root: Path, hit: Hit, flank: int, executable='blastdbcmd'):
    if hit.reference_key not in BY_KEY:raise InputError('The hit refers to an unknown assembly.')
    if flank<0 or flank>100000:raise InputError('Flank must be between 0 and 100,000 bp.')
    ref=BY_KEY[hit.reference_key];db=database_prefix(root,ref)
    start=max(1,min(hit.subject_start,hit.subject_end)-flank)
    end=min(hit.subject_length,max(hit.subject_start,hit.subject_end)+flank)
    if start>end:raise InputError('Invalid hit coordinates.')
    result=_run([executable,'-db',str(db),'-entry',hit.subject,'-range',f'{start}-{end}','-strand',hit.strand,'-outfmt','%s'])
    seq=''.join(result.split()).upper()
    if not seq or not set(seq)<=set('ACGTNRYKMSWBDHV') or len(seq)!=(end-start+1):
        raise InputError('Extracted sequence is empty, contains invalid bases, or has an unexpected length.')
    # Existing sequence workflow accepts ACGTN; represent IUPAC ambiguity as N,
    # never silently as an exact base. Tell users in the interface.
    ambiguous=sum(base not in 'ACGTN' for base in seq)
    seq=''.join(base if base in 'ACGTN' else 'N' for base in seq)
    return (f'{ref.key}|{hit.subject}:{start}-{end}|strand={hit.strand}|ambiguous_to_N={ambiguous}',seq)

def main():
    import argparse
    parser=argparse.ArgumentParser(description='Download and prepare curated SoyBase/LIS soybean assemblies.')
    parser.add_argument('--list',action='store_true',help='List available assemblies')
    parser.add_argument('--install',nargs='+',choices=list(BY_KEY),metavar='KEY',help='Download and index these assemblies')
    parser.add_argument('--directory',type=Path,default=DEFAULT_ROOT,help='Library directory (default: repository data/soybean)')
    parser.add_argument('--makeblastdb',default='makeblastdb',help='Path to NCBI makeblastdb')
    args=parser.parse_args()
    if args.list:
        for ref in REFERENCES:print(f'{ref.key:15} G. {ref.species:5} {ref.accession:22} {ref.assembly:18} ~{ref.approx_mb} MB compressed')
    if args.install:
        if not shutil.which(args.makeblastdb):parser.error('makeblastdb is required and was not found.')
        for key in args.install:
            ref=BY_KEY[key]
            download_and_prepare(ref,args.directory,print)
            build_database(ref,args.directory,args.makeblastdb,print)
    elif not args.list:parser.print_help()

if __name__=='__main__':
    main()
