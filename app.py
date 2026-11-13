#!/usr/bin/env python3
"""Primer Design Helper: local GUI for a soybean marker design workflow."""
from __future__ import annotations
import csv
import io
import json
import os
import random
import threading
import tkinter as tk
from pathlib import Path
from tkinter import ttk, filedialog, messagebox
from primer_helper import core

class Application(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Primer Design Helper')
        self.geometry('1180x800'); self.minsize(880,640)
        self.option_add('*Font', ('Arial', 10))
        style=ttk.Style(self)
        try: style.theme_use('clam')
        except tk.TclError: pass
        style.configure('Title.TLabel',font=('Arial',19,'bold'),foreground='#173d3b')
        style.configure('Sub.TLabel',foreground='#425d5b')
        style.configure('Action.TButton',font=('Arial',10,'bold'))
        self.fields={}; self.records=[]; self.variant_items=[]; self.designs=[]
        self.var_status=tk.StringVar(value='Ready. Your sequences stay on this computer.')
        top=ttk.Frame(self,padding=(16,14,16,6)); top.pack(fill='x')
        ttk.Label(top,text='Primer Design Helper',style='Title.TLabel').pack(anchor='w')
        ttk.Label(top,text='From aligned parental sequences to assay design, plate planning and validation',style='Sub.TLabel').pack(anchor='w',pady=(4,0))
        bar=ttk.Frame(self,padding=(16,0));bar.pack(fill='x')
        for label,callback in [('Open project',self.open_project),('Save project',self.save_project),('Load example',self.example)]:
            ttk.Button(bar,text=label,command=callback).pack(side='left',padx=(0,6))
        ttk.Label(bar,text='Local workspace • save a project to keep your work',style='Sub.TLabel').pack(side='right')
        self.nb=ttk.Notebook(self);self.nb.pack(fill='both',expand=True,padx=16,pady=10)
        self.make_sequences();self.make_design();self.make_specificity();self.make_bench();self.make_plate();self.make_validation();self.make_kasp()
        ttk.Label(self,textvariable=self.var_status,relief='sunken',anchor='w',padding=6).pack(fill='x',side='bottom')
        self.example()
    def tab(self,name):
        frame=ttk.Frame(self.nb,padding=12);self.nb.add(frame,text=name);return frame
    def textfield(self,parent,key,height=8,wrap='none'):
        frame=ttk.Frame(parent); frame.pack(fill='both',expand=True,pady=(3,8))
        text=tk.Text(frame,height=height,wrap=wrap,undo=True,font=('Consolas',10),borderwidth=1,relief='solid')
        sy=ttk.Scrollbar(frame,orient='vertical',command=text.yview);text.configure(yscrollcommand=sy.set)
        sy.pack(side='right',fill='y')
        if wrap == 'none':
            sx=ttk.Scrollbar(frame,orient='horizontal',command=text.xview);text.configure(xscrollcommand=sx.set);sx.pack(side='bottom',fill='x')
        text.pack(fill='both',expand=True)
        self.fields[key]=text;return text
    def value(self,key): return self.fields[key].get('1.0','end-1c')
    def put(self,key,text):
        w=self.fields[key];w.delete('1.0','end');w.insert('1.0',str(text))
    def header(self,parent,title,sub):
        ttk.Label(parent,text=title,font=('Arial',14,'bold')).pack(anchor='w')
        ttk.Label(parent,text=sub,wraplength=1060,justify='left').pack(anchor='w',pady=(2,8))
    def make_sequences(self):
        f=self.tab('1  Sequences & variants')
        self.header(f,'Compare the parental region','Paste genomic FASTA with two or more sequences. The first is the reference and the second is the comparison. Supply orthologous loci in the same orientation. The example is already aligned.')
        controls=ttk.Frame(f);controls.pack(fill='x')
        ttk.Button(controls,text='Open FASTA',command=self.open_fasta).pack(side='left',padx=(0,6))
        ttk.Button(controls,text='Align or realign with MAFFT',command=self.align).pack(side='left',padx=6)
        ttk.Button(controls,text='Find variants in aligned FASTA',command=self.find_variants).pack(side='left',padx=6)
        ttk.Label(controls,text='MAFFT executable').pack(side='left',padx=(16,4));self.mafft=tk.StringVar(value='mafft');ttk.Entry(controls,textvariable=self.mafft,width=16).pack(side='left')
        ttk.Label(f,text='FASTA input or aligned FASTA (MAFFT realignment removes existing gaps before aligning)' ).pack(anchor='w',pady=(8,0));self.textfield(f,'sequences',height=9)
        ttk.Label(f,text='Select a variant for primer design (coordinates refer to the alignment and first record)').pack(anchor='w')
        self.variant_list=tk.Listbox(f,height=8,exportselection=False,font=('Consolas',10));self.variant_list.pack(fill='both',expand=True,pady=5)
        ttk.Label(f,text='Already aligned? Select Find variants directly. For realignment, existing - gaps are removed first and may be placed differently; N (unknown bases) is retained. Inspect alignments and verify orthology.',wraplength=1050).pack(anchor='w')
    def make_design(self):
        f=self.tab('2  Primer3 design')
        self.header(f,'Design shared flanking primers','Select a variant in the first tab, then run Primer3 locally. The app checks each returned binding site across every aligned input sequence; it does not establish genome-wide specificity.')
        row=ttk.Frame(f);row.pack(fill='x')
        self.p3path=tk.StringVar(value='primer3_core');self.minlen=tk.StringVar(value='100');self.maxlen=tk.StringVar(value='400');self.npairs=tk.StringVar(value='5')
        for label,var,width in [('Primer3 executable',self.p3path,18),('Minimum product',self.minlen,7),('Maximum product',self.maxlen,7),('Pairs',self.npairs,5)]:
            ttk.Label(row,text=label).pack(side='left',padx=(4,3));ttk.Entry(row,textvariable=var,width=width).pack(side='left',padx=(0,10))
        ttk.Button(row,text='Run Primer3',command=self.run_design,style='Action.TButton').pack(side='left')
        ttk.Label(f,text='Candidates and per-parent binding checks').pack(anchor='w',pady=(10,0));self.textfield(f,'design_output',height=18,wrap='word')
        row=ttk.Frame(f);row.pack(fill='x');ttk.Label(row,text='Candidate number').pack(side='left');self.candidate=tk.StringVar(value='1');ttk.Entry(row,textvariable=self.candidate,width=5).pack(side='left',padx=5)
        ttk.Button(row,text='Use pair for specificity tab',command=self.copy_pair).pack(side='left',padx=5)
        ttk.Label(f,text='For a SNP, a common pair alone does not genotype the variant. For KASP, use the separate KASP planning notes on the validation tab.',wraplength=1050).pack(anchor='w',pady=7)
    def make_kasp(self):
        f=self.tab('7  KASP candidates')
        self.header(f,'Propose allele-specific SNP primers','Select a SNP on the first tab. This experimental path proposes two forward target-binding segments with the SNP at the 3′ end and uses Primer3 for a common reverse primer. Evaluate every candidate before ordering.')
        row=ttk.Frame(f);row.pack(fill='x')
        self.kasp_length=tk.StringVar(value='20')
        ttk.Label(row,text='Allele-specific binding length').pack(side='left');ttk.Entry(row,textvariable=self.kasp_length,width=6).pack(side='left',padx=6)
        ttk.Button(row,text='Propose KASP candidates',command=self.run_kasp,style='Action.TButton').pack(side='left',padx=10)
        ttk.Label(f,text='Target-specific segments and common primer').pack(anchor='w',pady=(8,0));self.textfield(f,'kasp_output',height=12,wrap='word')
        ttk.Label(f,text='Assay mix: 72× KASP assay mix with 2× KASP-TF master mix, using the vendor exact-volume example for 96-well plates.',wraplength=1060).pack(anchor='w')
        row2=ttk.Frame(f);row2.pack(fill='x',pady=5)
        self.kasp_wells=tk.StringVar(value='96');self.kasp_overage=tk.StringVar(value='10')
        ttk.Label(row2,text='Wells').pack(side='left');ttk.Entry(row2,textvariable=self.kasp_wells,width=6).pack(side='left',padx=5)
        ttk.Label(row2,text='Overage %').pack(side='left');ttk.Entry(row2,textvariable=self.kasp_overage,width=6).pack(side='left',padx=5)
        ttk.Button(row2,text='Calculate KASP mix',command=self.kasp_mix).pack(side='left',padx=8)
        self.textfield(f,'kasp_mix_output',height=6,wrap='word')
        ttk.Label(f,text='Not an order-ready assay. Universal 5′ tails, dye mapping, allele-specific Tm/dimer review, common-primer specificity, assay-mix formulation and wet-lab validation are required. This path does not design small-indel KASP assays.',wraplength=1060).pack(anchor='w')
    def run_kasp(self):
        def start():
            records=core.parse_fasta(self.value('sequences'),aligned=True)
            if records!=self.records:raise core.InputError('Sequences changed. Find variants again.')
            selected=self.variant_list.curselection()
            if not selected:raise core.InputError('Select a SNP on the first tab.')
            v=self.variant_items[selected[0]];length=int(self.kasp_length.get());minimum=int(self.minlen.get());maximum=int(self.maxlen.get());count=int(self.npairs.get());exe=self.p3path.get()
            def done(items):
                try:current=core.parse_fasta(self.value('sequences'),aligned=True)
                except core.InputError:current=[]
                if current!=records:self.var_status.set('KASP results discarded because sequences changed.');return
                out=['EXPERIMENTAL SNP CANDIDATES — sequence-specific segments only; not ready to order.','First record allele '+v.ref+'; second record allele '+v.alt+'.','']
                for item in items:
                    out += [f"PAIR {item['rank']}",f"  Allele {v.ref} target-binding segment: 5′-{item['allele_1']}-3′",f"  Allele {v.alt} target-binding segment: 5′-{item['allele_2']}-3′",f"  Common reverse primer: 5′-{item['reverse']}-3′; Primer3 Tm {item['reverse_tm']} °C",'  Common-primer binding: '+', '.join(x['name']+(' exact' if x['match'] else ' MISMATCH/INDEL') for x in item['reverse_match']),'  Expected product lengths: '+', '.join(f'{k} {v} bp' for k,v in item['sizes'].items()),'']
                out.append('The final base is the discriminating SNP. Check both allele-specific primer Tms, cross-dimers, local genome hits, non-focal variation, correct KASP tails and dye assignment. Tails alter synthesized primer length. Do not use without assay validation.')
                self.put('kasp_output','\n'.join(out));self.var_status.set('KASP candidate segments generated. Genome and bench validation required.')
            self.background(lambda:list(core.kasp_candidates(records,v,length,minimum,maximum,count,exe)),done)
        self.guard(start)
    def kasp_mix(self):
        def run():
            import math
            wells=int(self.kasp_wells.get());over=float(self.kasp_overage.get())
            if wells<=0 or over<0:raise core.InputError('Use positive wells and nonnegative overage.')
            n=math.ceil(wells*(1+over/100))
            self.put('kasp_mix_output',f'For {wells} wells prepare {n} equivalents:\n2× KASP-TF master mix: {5*n:.2f} µL\n72× assay mix: {0.14*n:.2f} µL\nMix and dispense 5.14 µL per well; add 4.86 µL DNA or NTC water separately for 10 µL total.\nFollow the assay sheet and instrument-specific requirements.')
        self.guard(run)
    def make_specificity(self):
        f=self.tab('3  Specificity')
        self.header(f,'Inspect potential genomic binding sites','Run BLAST against a local nucleotide database, then inspect both primers at all relevant loci. A single hit is not an amplicon; possible products depend on orientation and spacing.')
        row=ttk.Frame(f);row.pack(fill='x')
        self.blastpath=tk.StringVar(value='blastn');self.dbpath=tk.StringVar()
        ttk.Label(row,text='blastn').pack(side='left');ttk.Entry(row,textvariable=self.blastpath,width=18).pack(side='left',padx=5)
        ttk.Label(row,text='Database prefix').pack(side='left',padx=(12,2));ttk.Entry(row,textvariable=self.dbpath,width=55).pack(side='left',padx=4,fill='x',expand=True)
        ttk.Label(f,text='Primer sequences, both written 5′ to 3′; one line per primer').pack(anchor='w',pady=(8,0));self.textfield(f,'primers',height=2)
        row=ttk.Frame(f);row.pack(fill='x');ttk.Button(row,text='Search local database',command=self.run_blast,style='Action.TButton').pack(side='left')
        ttk.Button(row,text='Export hit table',command=lambda:self.export_text('blast_output','blast_hits.tsv')).pack(side='left',padx=7)
        ttk.Label(f,text='BLAST hits (tab separated)').pack(anchor='w',pady=(8,0));self.textfield(f,'blast_output',height=16)
        ttk.Label(f,text='Results depend on database completeness and search settings. Short-primer BLAST can miss sites; examine 3′ mismatches and paired products, and verify experimentally. Consider both parental assemblies and paralogues.',wraplength=1050).pack(anchor='w')
    def make_bench(self):
        f=self.tab('4  Mix calculator');self.header(f,'Primer stocks and reaction master mix','Use the actual oligo yield and reagent concentrations. These calculations do not choose polymerase conditions for you.')
        grid=ttk.Frame(f);grid.pack(fill='x')
        self.calcvars={}
        values=[('oligo_nmol','Delivered oligo (nmol)','9'),('stock_um','Resuspension target (µM)','100'),('work_um','Working concentration (µM)','10'),('work_vol','Working volume (µL)','100'),('total','Reaction volume (µL)','25'),('strength','Master mix strength (×)','2'),('fstock','Forward primer stock (µM)','10'),('rstock','Reverse primer stock (µM)','10'),('final','Final each primer (µM)','0.2'),('template','DNA per well (µL)','2'),('wells','Number of wells','96'),('overage','Overage (%)','10')]
        for i,(key,label,val) in enumerate(values):
            row=i//3;col=(i%3)*2
            ttk.Label(grid,text=label).grid(row=row,column=col,sticky='w',padx=(4,6),pady=5)
            self.calcvars[key]=tk.StringVar(value=val);ttk.Entry(grid,textvariable=self.calcvars[key],width=10).grid(row=row,column=col+1,sticky='w',padx=(0,16))
        ttk.Button(f,text='Calculate',command=self.calculate,style='Action.TButton').pack(anchor='w',pady=8)
        self.textfield(f,'mix_output',height=12,wrap='word')
        ttk.Label(f,text='This calculator is for a conventional two-primer reaction. KASP assay mixes have different chemistry and ratios; use the KASP instructions in the validation tab.',wraplength=1050).pack(anchor='w')
    def make_plate(self):
        f=self.tab('5  Plate planner');self.header(f,'Prepare a 96-well sample map','Enter one unique well label per line, in loading order. Wells fill down columns: A1, B1 … H1, A2. Explicitly include parents, true heterozygote when available, extraction blank and at least two NTCs.')
        ttk.Label(f,text='Sample and control labels').pack(anchor='w');self.textfield(f,'samples',height=10)
        row=ttk.Frame(f);row.pack(fill='x');ttk.Button(row,text='Preview plate',command=self.make_plate_map).pack(side='left');ttk.Button(row,text='Export CSV',command=self.export_plate).pack(side='left',padx=6)
        self.textfield(f,'plate_output',height=10)
    def make_validation(self):
        f=self.tab('6  Validate & record');self.header(f,'Validation, KASP and experimental notes','Document what was actually tested. Keep uncertain wells as no calls; no software result proves biological specificity.')
        self.textfield(f,'notes',height=18,wrap='word')
        ttk.Button(f,text='Insert assay guidance',command=self.insert_guidance).pack(anchor='w')
    def open_fasta(self):
        path=filedialog.askopenfilename(filetypes=[('FASTA','*.fa *.fasta *.fna *.fas'),('All files','*')])
        if path:
            try:self.put('sequences',Path(path).read_text(encoding='utf-8'));self.invalidate_sequences();self.var_status.set('FASTA loaded; align or find variants.')
            except OSError as e:messagebox.showerror('Open FASTA',str(e))
    def invalidate_sequences(self):
        self.records=[];self.variant_items=[];self.designs=[];self.variant_list.delete(0,'end')
    def guard(self,fn):
        try: fn()
        except (core.InputError,ValueError,IndexError,KeyError,OSError) as e:messagebox.showerror('Check input',str(e));self.var_status.set('Check input and try again.')
    def background(self,func,done):
        self.var_status.set('Running external program…')
        def worker():
            try:
                result=func();self.after(0,lambda:done(result))
            except Exception as exc:
                message=str(exc);self.after(0,lambda:messagebox.showerror('Program error',message));self.after(0,lambda:self.var_status.set('External program failed.'))
        threading.Thread(target=worker,daemon=True).start()
    def align(self):
        inp=self.value('sequences');exe=self.mafft.get()
        def done(result):
            # Do not overwrite input if the user changed it while MAFFT was running.
            if self.value('sequences') != inp:
                self.var_status.set('Alignment discarded because input changed while MAFFT was running.');return
            self.put('sequences',result);self.invalidate_sequences();self.find_variants();self.var_status.set('MAFFT realignment completed; existing gaps were reset. Inspect orthology and orientation.')
        self.background(lambda:core.mafft_align(inp,exe),done)
    def find_variants(self):
        def run():
            records=core.parse_fasta(self.value('sequences'),aligned=True);v=core.variants(records[0][1],records[1][1])
            self.records=records;self.variant_items=v;self.designs=[];self.variant_list.delete(0,'end')
            for item in v:self.variant_list.insert('end',item.label)
            self.var_status.set(f'{len(v)} variant groups found comparing {records[0][0]} with {records[1][0]}. Select one.')
        self.guard(run)
    def run_design(self):
        def start():
            current=core.parse_fasta(self.value('sequences'),aligned=True)
            if current != self.records:raise core.InputError('Sequences changed. Find variants again before designing.')
            sel=self.variant_list.curselection()
            if not sel:raise core.InputError('Select a variant on the first tab.')
            selected=self.variant_items[sel[0]];records=self.records[:]
            minimum,maximum,count=int(self.minlen.get()),int(self.maxlen.get()),int(self.npairs.get());exe=self.p3path.get()
            def done(result):
                if self.value('sequences') != core.fasta(records):
                    # Compare parsed data rather than formatting to support wrapped FASTA.
                    try:unchanged=core.parse_fasta(self.value('sequences'),aligned=True)==records
                    except core.InputError:unchanged=False
                    if not unchanged:self.var_status.set('Primer results discarded because sequences changed.');return
                self.designs=result;out=[]
                for d in result:
                    out.append(f"PAIR {d['rank']}\n  Forward 5′–3′: {d['forward']}   Tm {d['forward_tm']} °C\n  Reverse 5′–3′: {d['reverse']}   Tm {d['reverse_tm']} °C\n  Reference positions: {d['ref_left_1based']} to {d['ref_right_1based']}\n  Predicted lengths: "+', '.join(f'{k} {v} bp' for k,v in d['sizes'].items()))
                    for typ,key in [('Forward','forward_match'),('Reverse','reverse_match')]:
                        out.append('  '+typ+' binding: '+', '.join(x['name']+(' exact' if x['match'] else ' MISMATCH/INDEL') for x in d[key]))
                    out.append('')
                self.put('design_output','\n'.join(out));self.var_status.set(f'{len(result)} pairs returned. Screen genome-wide and validate in the lab.')
            self.background(lambda:core.design(records,selected,minimum,maximum,count,exe),done)
        self.guard(start)
    def copy_pair(self):
        def run():
            i=int(self.candidate.get())-1
            if i<0 or i>=len(self.designs):raise core.InputError('Choose an available candidate number.')
            d=self.designs[i];self.put('primers',d['forward']+'\n'+d['reverse']);self.nb.select(2)
        self.guard(run)
    def run_blast(self):
        def start():
            seq=self.value('primers').splitlines()
            if len(seq)!=2:raise core.InputError('Enter exactly two primer sequences, one per line.')
            f,r=[s.strip().upper() for s in seq];db=self.dbpath.get();exe=self.blastpath.get()
            def done(output):
                header='query\tsubject\tidentity%\tlength\tmismatches\tgaps\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tquery_seq\tsubject_seq\n'
                self.put('blast_output',header+output);self.var_status.set('BLAST complete. Review paired sites and database coverage; not a specificity guarantee.')
            self.background(lambda:core.blast_primers(f,r,db,exe),done)
        self.guard(start)
    def calculate(self):
        def run():
            v={k:float(x.get()) for k,x in self.calcvars.items()};nmol=v['oligo_nmol']; stock=v['stock_um']
            vol,stockvol,water=core.oligo_resuspension(nmol,stock,stock,v['work_um'],v['work_vol'])
            if not v['wells'].is_integer():raise core.InputError('Well count must be an integer.')
            n,parts=core.reaction_mix(v['total'],v['strength'],v['fstock'],v['rstock'],v['final'],v['template'],int(v['wells']),v['overage'])
            text=f"OLIGO\nAdd {vol:.2f} µL buffer to {nmol:g} nmol for {stock:g} µM stock.\nFor {v['work_vol']:g} µL working solution: {stockvol:.2f} µL stock + {water:.2f} µL water.\n\nMASTER MIX — {int(v['wells'])} wells; prepare {n} equivalents\n"
            for name,x in parts:text+=f'{name:32} {x:8.3f} µL per well   {x*n:9.2f} µL for {n} equivalents'+(' (add separately)' if name.startswith('Template') else '')+'\n'
            text+='\nDispense shared components together, then add DNA separately; use water in no-template wells.'
            self.put('mix_output',text)
        self.guard(run)
    def make_plate_map(self):
        def run():
            data=core.plate_csv(self.value('samples'));rows=list(csv.DictReader(io.StringIO(data)));mapping={r['well']:r['sample_or_control'] for r in rows}
            lines=['       '+''.join(f'{i:<17}' for i in range(1,13))]
            for r in 'ABCDEFGH': lines.append(f'{r:<7}'+''.join(f'{mapping.get(r+str(c),"—")[:15]:<17}' for c in range(1,13)))
            self.put('plate_output','\n'.join(lines));self.var_status.set(f'{len(rows)} wells mapped. Confirm controls before export.')
        self.guard(run)
    def export_plate(self):
        def run():
            data=core.plate_csv(self.value('samples'));path=filedialog.asksaveasfilename(defaultextension='.csv',initialfile='plate_map.csv',filetypes=[('CSV','*.csv')])
            if path:Path(path).write_text(data,encoding='utf-8',newline='');self.var_status.set('Plate map saved.')
        self.guard(run)
    def export_text(self,key,name):
        path=filedialog.asksaveasfilename(initialfile=name)
        if path:Path(path).write_text(self.value(key),encoding='utf-8');self.var_status.set('File saved.')
    def save_project(self):
        path=filedialog.asksaveasfilename(defaultextension='.json',initialfile='primer_project.json',filetypes=[('Project JSON','*.json')])
        if not path:return
        data={'format':'primer-design-helper','version':1,'fields':{k:self.value(k) for k in self.fields},'settings':{'mafft':self.mafft.get(),'primer3':self.p3path.get(),'blastn':self.blastpath.get(),'db':self.dbpath.get(),'min':self.minlen.get(),'max':self.maxlen.get(),'pairs':self.npairs.get(),'calc':{k:v.get() for k,v in self.calcvars.items()}}}
        try:Path(path).write_text(json.dumps(data,indent=2),encoding='utf-8');self.var_status.set('Project saved locally.')
        except OSError as e:messagebox.showerror('Save project',str(e))
    def open_project(self):
        path=filedialog.askopenfilename(filetypes=[('Project JSON','*.json')]);
        if not path:return
        def run():
            data=json.loads(Path(path).read_text(encoding='utf-8'))
            if data.get('format')!='primer-design-helper' or data.get('version')!=1:raise core.InputError('Unrecognized project format.')
            self.invalidate_sequences()
            for k,v in data.get('fields',{}).items():
                if k in self.fields and isinstance(v,str):self.put(k,v)
            settings=data.get('settings',{})
            for key,var in [('mafft',self.mafft),('primer3',self.p3path),('blastn',self.blastpath),('db',self.dbpath),('min',self.minlen),('max',self.maxlen),('pairs',self.npairs)]:
                if key in settings:var.set(settings[key])
            for k,v in settings.get('calc',{}).items():
                if k in self.calcvars:self.calcvars[k].set(v)
            self.var_status.set('Project loaded. Find variants again to reconstruct selections; rerun Primer3 to select a candidate.')
        self.guard(run)
    def example(self):
        rng=random.Random(21);ref=''.join(rng.choice('ACGT') for _ in range(450));alt=ref[:180]+('A' if ref[180]!='A' else 'G')+ref[181:]
        alt=alt[:260]+'ACGTTGCAACGTACGTTGCA'+alt[260:];aligned_ref=ref[:260]+'-'*20+ref[260:]
        self.put('sequences',f'>hypothetical_Gmax\n{aligned_ref}\n>hypothetical_Gsoja\n{alt}\n')
        self.put('samples','\n'.join([f'S{i:02}' for i in range(1,91)]+['Gmax_parent','Gsoja_parent','verified_heterozygote','extraction_blank','NTC_1','NTC_2']))
        self.put('notes','Example data are synthetic and are not orderable or biologically validated.\n\nRecord assembly/version, coordinate, strand, parent identities, variant evidence, homology checks, primer sequences and conditions, controls, instrument, trace/gel locations, and call decisions.\n')
        self.invalidate_sequences();self.find_variants();self.calculate();self.make_plate_map();self.var_status.set('Synthetic example loaded. Replace it with your own sequences.')
    def insert_guidance(self):
        text='''\n\nKASP AND HRM CHECKLIST\n• KASP: use two allele-specific primers and one common primer, compatible universal tails, assay mix, master mix and assay-specific dye mapping. Use vendor-specific design and instructions; common Primer3 pairs are not a complete KASP assay. Screen binding sites against both parents and paralogues.\n• For KASP use controls for both homozygotes and a sequence-verified heterozygote where possible, plus extraction blank and at least two NTCs. Interpret FAM/HEX clusters against controls. Do not force ambiguous clusters.\n• Indel gel: compare predicted parental fragment sizes with ladder resolution and validate heterozygotes. Additional polymorphisms can change size.\n• HRM: use appropriate dye, instrument, calibration and analysis; compare replicates of verified parental homozygotes and heterozygotes. Other variants can alter melt profiles; confirm representative groups by sequencing.\n• Document sequence-confirmed controls, repeats, concordance, missingness, allele dropout and off-target evidence before production use.\n\nSources: https://www.biosearchtech.com/how-does-kasp-work ; https://primer3.org/manual.html ; https://www.ncbi.nlm.nih.gov/tools/primer-blast/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC6198901/\n'''
        self.fields['notes'].insert('end',text)

if __name__=='__main__':
    Application().mainloop()
