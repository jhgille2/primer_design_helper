"""Detect external tools and optionally install into a dedicated Conda prefix."""
from __future__ import annotations
import os
import platform
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from .core import InputError

TOOL_NAMES=('mafft','primer3_core','blastn','makeblastdb','blastdbcmd')
PACKAGES={'mafft':'mafft','primer3_core':'primer3','blastn':'blast','makeblastdb':'blast','blastdbcmd':'blast'}
OFFICIAL_LINKS={
    'MAFFT':'https://mafft.cbrc.jp/alignment/software/',
    'Primer3':'https://primer3.org/releases.html',
    'BLAST+':'https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html',
    'Conda (Miniforge)':'https://github.com/conda-forge/miniforge',
}
DEFAULT_PREFIX=Path(__file__).resolve().parent.parent/'data'/'tools'/'bioinformatics'

@dataclass(frozen=True)
class ToolStatus:
    name: str
    requested: str
    path: str | None
    usable: bool
    detail: str

def managed_paths(prefix: Path) -> dict[str,str]:
    # Bioconda's supported platforms are Linux and macOS; keep the Windows
    # lookup for manually populated prefixes, without assuming package support.
    folder=Path(prefix)/('Scripts' if os.name=='nt' else 'bin')
    return {name:str(folder/(name+('.exe' if os.name=='nt' else ''))) for name in TOOL_NAMES}

def probe(name: str, requested: str, timeout: int=12) -> ToolStatus:
    requested=requested.strip()
    if name not in TOOL_NAMES:raise InputError('Unknown software dependency.')
    path=shutil.which(requested) if requested else None
    if path is None:return ToolStatus(name,requested,None,False,'Not found. Enter a path or install it.')
    if name=='primer3_core':
        # primer3_core has no portable --version flag. Check its Boulder-IO
        # interface with an empty record rather than pretending version output.
        args=[path];input_text='=\n'
    else:
        args=[path,'--version' if name=='mafft' else '-version'];input_text=''
    try:
        result=subprocess.run(args,input=input_text,text=True,capture_output=True,timeout=timeout,check=False)
    except (OSError,subprocess.TimeoutExpired) as exc:return ToolStatus(name,requested,path,False,f'Could not run: {exc}')
    output=(result.stdout+'\n'+result.stderr).strip()
    # Primer3 may reject an empty input record while still revealing PRIMER_ERROR.
    if name=='primer3_core':
        usable=('PRIMER_ERROR=' in output or 'PRIMER_LEFT_' in output or 'PRIMER_PAIR_' in output) and result.returncode in (0,1)
        detail='Boulder-IO response detected' if usable else (output[:220] or 'No recognizable Primer3 response.')
    else:
        usable=result.returncode==0
        detail=next((line.strip() for line in output.splitlines() if line.strip()),'No version output')[:220]
    return ToolStatus(name,requested,path,usable,detail)

def package_plan(statuses: list[ToolStatus]) -> list[str]:
    return list(dict.fromkeys(PACKAGES[status.name] for status in statuses if not status.usable))

def install_missing(statuses: list[ToolStatus], conda: str, prefix: Path, progress=None):
    """Explicitly install missing tools via conda/mamba into the app prefix.

    This requires an existing conda-family executable on Linux or macOS. No
    elevated privileges, shell execution, or modification of base is used.
    """
    if platform.system() not in ('Linux','Darwin'):
        raise InputError('Automatic Bioconda installation requires Linux or macOS. On Windows, use the official download links and browse to executables, or run the entire app in WSL.')
    manager=shutil.which(conda.strip()) if conda.strip() else None
    if not manager:raise InputError('Conda or Mamba was not found. Install Miniforge, or enter the full path to an existing conda/mamba executable.')
    packages=package_plan(statuses)
    if not packages:return {}
    prefix=Path(prefix).expanduser().resolve()
    prefix.parent.mkdir(parents=True,exist_ok=True)
    exists=(prefix/'conda-meta'/'history').is_file()
    if prefix.exists() and not exists:raise InputError(f'The target {prefix} exists but is not a valid Conda environment; choose a different location or move it aside.')
    cmd=[manager,'install' if exists else 'create','--yes','--prefix',str(prefix),'--override-channels','--channel','conda-forge','--channel','bioconda','--strict-channel-priority',*packages]
    if progress:progress('Installing '+', '.join(packages)+' into the app environment. This can take several minutes...')
    try:result=subprocess.run(cmd,text=True,capture_output=True,timeout=3600,check=False)
    except (OSError,subprocess.TimeoutExpired) as exc:raise InputError(f'Conda installation could not finish: {exc}') from exc
    if result.returncode:raise InputError('Conda installation failed. '+(result.stderr or result.stdout)[-3000:])
    managed=managed_paths(prefix); found={}
    for status in statuses:
        if status.usable:continue
        checked=probe(status.name,managed[status.name])
        if checked.usable:found[status.name]=checked.path
    missing=[s.name for s in statuses if not s.usable and s.name not in found]
    if missing:raise InputError('Installation completed, but these executables could not be run: '+', '.join(missing)+'. Check Conda output and platform support.')
    if progress:progress('Installed and checked: '+', '.join(found))
    return found
