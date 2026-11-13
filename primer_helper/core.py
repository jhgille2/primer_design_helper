"""Sequence and reaction calculations for the local primer design helper.

Coordinates in the Python API are zero-based, half-open unless stated otherwise.
"""
from __future__ import annotations
import csv
import io
import math
import re
import shutil
import subprocess
from dataclasses import dataclass
from itertools import zip_longest

class InputError(ValueError):
    pass

def parse_fasta(text: str, *, aligned: bool = False) -> list[tuple[str, str]]:
    records = []
    name = None
    chunks = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name is not None:
                records.append((name, ''.join(chunks)))
            name, chunks = line[1:].strip(), []
            if not name:
                raise InputError('Every FASTA record must have a name.')
        else:
            if name is None:
                raise InputError('FASTA must start with a >header.')
            chunks.append(''.join(line.split()).upper())
    if name is not None:
        records.append((name, ''.join(chunks)))
    if len(records) < 2:
        raise InputError('Provide at least two FASTA records.')
    if len({n for n, _ in records}) != len(records):
        raise InputError('FASTA names must be unique.')
    allowed = set('ACGTN-') if aligned else set('ACGTN')
    for name, seq in records:
        if not seq or not set(seq) <= allowed:
            raise InputError(f'{name}: use only A, C, G, T, N' + (' and gaps (-).' if aligned else '.'))
    if aligned and len({len(seq) for _, seq in records}) != 1:
        raise InputError('Aligned FASTA sequences must all have the same length.')
    return records

def fasta(records):
    return ''.join(f'>{name}\n{seq}\n' for name, seq in records)

def run_tool(command: list[str], *, input_text: str = '', timeout: int = 120) -> str:
    if not shutil.which(command[0]):
        raise InputError(f'Cannot find {command[0]!r}. Install it or provide its full executable path.')
    try:
        proc = subprocess.run(command, input=input_text, text=True, capture_output=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        raise InputError(f'{command[0]} exceeded {timeout} seconds.') from exc
    except OSError as exc:
        raise InputError(f'Could not start {command[0]}: {exc}') from exc
    if proc.returncode:
        raise InputError(f'{command[0]} failed: {proc.stderr.strip()[:2000]}')
    return proc.stdout

def mafft_align(text: str, executable: str = 'mafft') -> str:
    records = parse_fasta(text)
    result = run_tool([executable, '--auto', '-'], input_text=fasta(records), timeout=300)
    aligned = parse_fasta(result, aligned=True)
    # MAFFT normally keeps sequence order, but explicitly restore it by unique names.
    by_name = dict(aligned)
    if set(by_name) != {name for name, _ in records}:
        raise InputError('MAFFT output names do not match the input names.')
    return fasta([(name, by_name[name]) for name, _ in records])

@dataclass(frozen=True)
class Variant:
    start: int  # alignment start, inclusive
    end: int    # alignment end, exclusive
    ref_start: int  # number of non-gap reference bases before start
    ref: str
    alt: str
    kind: str
    @property
    def label(self):
        return f'alignment {self.start + 1}–{self.end} | ref offset {self.ref_start} (0-based) | {self.kind} | {self.ref or "−"} → {self.alt or "−"}'

def variants(ref: str, alt: str) -> list[Variant]:
    if len(ref) != len(alt):
        raise InputError('Aligned sequences must have the same length.')
    out = []
    refpos = 0
    i = 0
    while i < len(ref):
        a,b = ref[i],alt[i]
        if a == b or a == 'N' or b == 'N':
            refpos += (a != '-')
            i += 1
            continue
        start, start_ref = i, refpos
        gap = a == '-' or b == '-'
        while i < len(ref):
            a,b = ref[i],alt[i]
            if a == b or a == 'N' or b == 'N' or ((a == '-' or b == '-') != gap):
                break
            refpos += (a != '-')
            i += 1
        r, q = ref[start:i].replace('-', ''), alt[start:i].replace('-', '')
        kind = ('insertion in comparison' if not r else 'deletion in comparison' if not q else 'complex indel') if gap else ('SNP' if len(r)==1 else 'substitution block')
        out.append(Variant(start,i,start_ref,r,q,kind))
    return out

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def primer3_input(ref_aligned: str, variant: Variant, minimum: int, maximum: int, count: int) -> str:
    ref = ref_aligned.replace('-', '')
    if minimum < 40 or maximum < minimum or maximum > 10000 or not 1 <= count <= 50:
        raise InputError('Use a product range of 40–10,000 bp and 1–50 pairs.')
    if len(ref) < minimum:
        raise InputError('Reference sequence is shorter than the minimum product size.')
    # A reference-gap insertion has no reference base: target one flanking base, if available.
    start = min(variant.ref_start, len(ref)-1)
    target_len = max(1, len(variant.ref))
    return (f'SEQUENCE_ID=selected_marker\nSEQUENCE_TEMPLATE={ref}\nSEQUENCE_TARGET={start},{target_len}\n'
            f'PRIMER_TASK=generic\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_RIGHT_PRIMER=1\n'
            f'PRIMER_NUM_RETURN={count}\nPRIMER_PRODUCT_SIZE_RANGE={minimum}-{maximum}\n'
            'PRIMER_MIN_SIZE=18\nPRIMER_OPT_SIZE=20\nPRIMER_MAX_SIZE=25\n'
            'PRIMER_MIN_TM=57.0\nPRIMER_OPT_TM=60.0\nPRIMER_MAX_TM=63.0\n'
            'PRIMER_EXPLAIN_FLAG=1\n=\n')

def parse_primer3(text: str) -> dict:
    out = {}
    for line in text.splitlines():
        if '=' in line and line != '=':
            k,v = line.split('=', 1)
            out[k] = v
    if out.get('PRIMER_ERROR'):
        raise InputError('Primer3: ' + out['PRIMER_ERROR'])
    return out

def binding_check(records: list[tuple[str,str]], start: int, end: int, expected: str) -> list[dict]:
    # [start,end) is on ungapped reference plus strand. Include alignment gaps
    # between the first and last reference bases, so indels are not hidden.
    ref = records[0][1]
    cols = [i for i,base in enumerate(ref) if base != '-']
    if start < 0 or end > len(cols) or start >= end:
        raise InputError('Primer position is outside the alignment.')
    left,right = cols[start], cols[end-1] + 1
    result = []
    for name, seq in records:
        observed = seq[left:right].replace('-', '')
        result.append({'name':name, 'match':observed == expected, 'observed':observed})
    return result

def design(records, v: Variant, minimum: int, maximum: int, count: int, executable='primer3_core', fixed_left: str | None = None) -> list[dict]:
    data = primer3_input(records[0][1], v, minimum, maximum, count)
    if fixed_left is not None:
        # Primer3 requires SEQUENCE_TARGET to be flanked, so it cannot contain
        # the SNP that intentionally lies in this allele-specific primer.
        data = re.sub(r'^SEQUENCE_TARGET=.*\n', '', data, flags=re.MULTILINE)
        primer_start = v.ref_start - len(fixed_left) + 1
        data = data.replace('SEQUENCE_ID=selected_marker\n', f'SEQUENCE_ID=selected_marker\nSEQUENCE_PRIMER={fixed_left}\nSEQUENCE_FORCE_LEFT_START={primer_start}\nSEQUENCE_FORCE_LEFT_END={v.ref_start}\n')
    parsed = parse_primer3(run_tool([executable], input_text=data))
    n = int(parsed.get('PRIMER_PAIR_NUM_RETURNED','0'))
    if n == 0:
        details = '; '.join(parsed.get(k,'') for k in ('PRIMER_PAIR_EXPLAIN','PRIMER_LEFT_EXPLAIN','PRIMER_RIGHT_EXPLAIN') if parsed.get(k))
        raise InputError('Primer3 returned no pairs. ' + details)
    results = []
    for i in range(n):
        lp,ll = map(int,parsed[f'PRIMER_LEFT_{i}'].split(','))
        # Primer3 right coordinate is the rightmost plus-strand coordinate.
        rp,rl = map(int,parsed[f'PRIMER_RIGHT_{i}'].split(','))
        left = parsed[f'PRIMER_LEFT_{i}_SEQUENCE']
        right = parsed[f'PRIMER_RIGHT_{i}_SEQUENCE']
        if fixed_left is not None and (left != fixed_left or lp + ll - 1 != v.ref_start):
            continue
        lcheck = binding_check(records,lp,lp+ll,left)
        rcheck = binding_check(records,rp-rl+1,rp+1,revcomp(right))
        # Product size for each aligned sequence between mapped terminal binding bases.
        refcols=[j for j,b in enumerate(records[0][1]) if b!='-']
        a,b=refcols[lp],refcols[rp]+1
        if not (a <= v.start and v.end <= b):
            continue
        if fixed_left is None and (v.start < refcols[lp+ll-1]+1 or v.end > refcols[rp-rl+1]):
            continue
        sizes={name:len(seq[a:b].replace('-','')) for name,seq in records}
        results.append({'rank':i+1,'forward':left,'reverse':right,'forward_tm':parsed.get(f'PRIMER_LEFT_{i}_TM','?'),
                        'reverse_tm':parsed.get(f'PRIMER_RIGHT_{i}_TM','?'),'sizes':sizes,
                        'forward_match':lcheck,'reverse_match':rcheck,
                        'ref_left_1based':lp+1,'ref_right_1based':rp+1})
    if not results:
        raise InputError('Primer3 returned candidates, but none span the complete alignment variant with the requested primer positions.')
    return results

def kasp_candidates(records, v: Variant, primer_length: int, minimum: int, maximum: int, count: int, executable='primer3_core') -> list[dict]:
    """Propose forward allele-specific binding segments and Primer3 common reverse primers.

    The SNP is at the final 3' base. No universal tails, mismatch engineering,
    proprietary assay chemistry or performance claims are supplied.
    """
    if v.kind != 'SNP' or len(v.ref) != 1 or len(v.alt) != 1 or v.ref not in 'ACGT' or v.alt not in 'ACGT':
        raise InputError('This initial KASP path supports a single A/C/G/T SNP only; small indels require specialized assay design.')
    if not 18 <= primer_length <= 25:
        raise InputError('Allele-specific binding length must be 18–25 bases.')
    ref=records[0][1].replace('-', '')
    start=v.ref_start-primer_length+1
    if start < 0:
        raise InputError('There is not enough upstream reference sequence for the allele-specific primer.')
    f1=ref[start:v.ref_start+1]
    f2=f1[:-1]+v.alt
    # Primer sites are assessed on each alignment, including internal gaps.
    checks=binding_check(records,start,v.ref_start+1,f1)
    second=checks[1]['observed']
    if second != f2:
        raise InputError('The comparison allele-specific binding site contains additional differences or an indel. Consider the opposite orientation or a specialized assay design.')
    for item in design(records,v,minimum,maximum,count,executable,fixed_left=f1):
        item['allele_1']=f1;item['allele_2']=f2
        item['allele_1_binding']=checks
        yield item


def blast_primers(forward: str, reverse: str, db: str, executable='blastn', max_hits=1000) -> str:
    for seq in (forward,reverse):
        if not seq or not set(seq.upper()) <= set('ACGT'):
            raise InputError('Primer sequences must contain only A, C, G and T.')
    if not db.strip():
        raise InputError('Provide a local BLAST nucleotide database prefix.')
    return run_tool([executable,'-task','blastn-short','-db',db,'-query','-',
                     '-max_target_seqs',str(max_hits),'-dust','no',
                     '-outfmt','6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq'],
                    input_text=f'>forward\n{forward.upper()}\n>reverse\n{reverse.upper()}\n', timeout=300)

def reaction_mix(total: float, master_factor: float, fstock: float, rstock: float, final: float, template: float, wells: int, overage: float):
    if min(total,master_factor,fstock,rstock,final,wells) <= 0 or template < 0 or overage < 0 or master_factor <= 1:
        raise InputError('Enter positive volumes and concentrations, nonnegative overage, and a master mix stronger than 1×.')
    n = math.ceil(wells*(1+overage/100))
    master = total/master_factor
    fw = final*total/fstock
    rv = final*total/rstock
    water = total-master-fw-rv-template
    if water < -1e-9:
        raise InputError('The component volumes exceed the total reaction volume.')
    return n, [('2× or specified master mix',master),('Forward primer',fw),('Reverse primer',rv),('Water',max(0.0,water)),('Template added separately',template)]

def oligo_resuspension(nmol: float, target_um: float, stock_um: float, work_um: float, work_volume_ul: float):
    if min(nmol,target_um,stock_um,work_um,work_volume_ul)<=0 or work_um>stock_um:
        raise InputError('Enter positive concentrations and volumes; working concentration must not exceed stock.')
    volume=nmol*1000/target_um
    stock=work_um*work_volume_ul/stock_um
    return volume,stock,work_volume_ul-stock

def plate_csv(samples: str) -> str:
    ids=[line.strip() for line in samples.splitlines() if line.strip()]
    if not ids:
        raise InputError('Enter one sample or control name per line.')
    if len(ids)>96:
        raise InputError('At most 96 wells are supported per plate.')
    if len(set(ids)) != len(ids):
        raise InputError('Well labels must be unique; add replicate suffixes to repeated IDs.')
    result=io.StringIO(newline=''); w=csv.writer(result); w.writerow(['well','sample_or_control'])
    for i,name in enumerate(ids):
        w.writerow([f'{chr(65+i%8)}{i//8+1}',name])
    return result.getvalue()
