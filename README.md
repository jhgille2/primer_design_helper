# Primer Design Helper

A local desktop application for developing PCR genotyping assays from candidate genomic regions. It is oriented toward a *Glycine max* × *Glycine soja* project, but accepts other DNA sequences. The app keeps sequence data on your computer and never submits it to external services.

> **Research prototype.** Review the underlying alignments, genomic specificity and experimental controls. This software does not establish orthology or prove that an assay is locus-specific. Candidate primers and genotype calls require independent and experimental validation before use.

## What it does

1. Accepts two or more genomic FASTA sequences, calls locally installed MAFFT for alignment or accepts aligned FASTA, and finds differences between the first and second records. Additional records participate in primer-binding checks. The supplied example is synthetic.
2. Runs local `primer3_core` to propose flanking primer pairs around a selected SNP or indel. It reports calculated primer Tms from Primer3, product sizes within the supplied alignments, and whether each primer binding site matches each provided sequence.
3. For a single SNP, proposes two allele-specific forward **binding segments** with the polymorphism at the 3′ end and a common reverse primer selected by Primer3. This is an experimental, plus-strand path. It does not design a complete order-ready KASP assay, automatically calculate the second primer's Tm, select universal tails or design small-indel KASP assays. Different allele-specific binding-site variation is rejected; consider the opposite orientation or specialist design.
4. Runs `blastn -task blastn-short` against a user-supplied local nucleotide database and exports individual primer hits. This is an aid to manual analysis: it does **not** pair hits, predict all off-target amplicons, prove the intended hit is orthologous or establish genome-wide specificity. BLAST can miss short primer matches and reports only the configured maximum of 1,000 target sequences.
5. Calculates stock dilutions and conventional two-primer reactions, and provides a specific example KASP mix calculation. Produces a 96-well, column-major sample/control CSV and saves a local JSON project with notes and displayed results.
6. Provides prompts for parental, heterozygote and negative controls, gel/HRM validation and KASP cluster review. It does **not** control a thermal cycler, read fluorescence/melt curves, make genotype calls, evaluate qPCR standard curves or automatically assess lab specificity.

## Install and run

You need **Python 3.10 or later with Tkinter**. The GUI and core calculations use the standard library; no Python packages are required. A graphical desktop is necessary, including when launching from a remote server. Check your installation:

```sh
python -m tkinter
```

This should open a small test window. On some Linux installations Tkinter is packaged separately; use the matching system Python/Tk package. On Windows, `py -m tkinter` may be appropriate. See [Python's Tkinter documentation](https://docs.python.org/3/library/tkinter.html).

Clone this repository and launch:

```sh
git clone https://github.com/jhgille2/primer_design_helper.git
cd primer_design_helper
python app.py
```

The GitHub account must have access to the private repository and Git must be authenticated. On Windows, substitute `py app.py` if needed. Alternatively, download the repository source using GitHub's interface and launch `app.py` from the extracted folder.

Install the following external programs for the corresponding actions. They must be available on your `PATH`, or enter the full path to the executable in the app:

- **MAFFT** for alignment: [official MAFFT downloads and instructions](https://mafft.cbrc.jp/alignment/software/). You can instead paste a pre-aligned FASTA.
- **Primer3** (`primer3_core`) for primer design: [Primer3](https://primer3.org/) and [manual](https://primer3.org/manual.html). Ensure the thermodynamic parameter files can be located by your installed build; consult its installation instructions.
- **NCBI BLAST+** (`blastn` and `makeblastdb`) for local binding-site searches: [BLAST+ documentation](https://www.ncbi.nlm.nih.gov/books/NBK279690/). A nucleotide FASTA alone is not a BLAST database. Create a database once, for example:

  ```sh
  makeblastdb -in soybean_genome.fasta -dbtype nucl -out soybean_db
  ```

  Enter the **database prefix** `soybean_db` or its full path in the app, not the name of a particular `.nhr`, `.nin`, or `.nsq` file. Search each parent/reference assembly that is relevant; record assembly names and versions. The app will not download genomes.

External executables and databases are independent installations and are not included in this repository. An unavailable executable gives an explanatory error. Long-running tools execute in a background thread.

## A first assay

1. Identify the actual parental lines, locus, variant and intended readout. Obtain genomic sequences for orthologous regions with substantial flanking sequence. Confirm assembly/version, strand and coordinates. A public reference is a proxy unless it represents your actual parent.
2. On **Sequences & variants**, open or paste FASTA. Put the design/reference parent first and comparison parent second. Align with MAFFT or provide equal-length aligned sequences. Inspect orientation, gaps, N bases, repeats, gene structure, genomic neighborhood and orthology independently. Press **Find variants** and select one. Alignment positions shown are one-based; the displayed reference offset is zero-based, and for insertions marks a boundary. Consecutive differences are grouped and an N is not called as a variant.
3. On **Primer3 design**, set the product size range appropriate for gel/fragment sizing, sequencing or HRM, then run Primer3. Review all sequence comparisons; a binding-site mismatch, especially at a 3′ end, can cause allele dropout. Product lengths are derived from the supplied alignment and are not a whole-genome prediction.
4. On **Specificity**, enter a local BLAST database prefix and query a candidate pair. Review both primers, orientation, product distance, 3′ matches and paralogues across suitable assemblies using a complementary in-silico PCR or Primer-BLAST workflow. An empty hit table is not proof of specificity.
5. For a SNP KASP candidate, select a *single SNP* and use **KASP candidates**. Review both target-specific segments, the common primer, Tms and potential dimers. Obtain the correct tails, dye assignments, primer ratios and chemistry from the assay vendor. The default conventional-PCR mix is not a KASP mix.
6. Use **Mix calculator** and **Plate planner** for working solutions, master-mix overage and well mapping. Include both parents, a sequence-confirmed heterozygote when available, extraction blank and at least two NTCs. A mixture of parental DNA is not a substitute for validating a real heterozygote.
7. Test the assay experimentally. Sequence representative products, compare with a trusted independent genotype, inspect gel bands or appropriate HRM/FAM/HEX profiles, examine repeated controls and samples, and retain unclear results as no calls. Standard melt curves are not automatically suitable for HRM genotyping. Keep your records in **Validate & record**, and save your project.

The built-in SNP example is synthetic. **Do not order its primers or use it as a biological result.**

## Project and data behavior

- A project JSON file contains the sequences, screen text, results, local executable paths, database path, plate labels, notes and calculator settings; treat it as sensitive if your material is unpublished.
- Loading a project restores displayed results but not executable selection state. Click **Find variants** and rerun Primer3 before choosing a candidate to send to the specificity tab.
- Changes are not autosaved. Use **Save project** to retain work. Long-running results are discarded if the sequence input changes before completion.
- The application is local; reference and software downloads are separate actions initiated by you.

## Assumptions and limitations

- Input sequences must already refer to orthologous genomic regions in the same orientation, and aligned FASTA must be equal length and contain A/C/G/T/N/gaps. Only the first two records are used for variant discovery. Unknown bases are ignored, not interpreted as evidence of equality.
- For a reference-gap insertion, Primer3 targets an adjacent reference base and the app checks the returned products span the full insertion alignment. For shared primers it also filters candidates overlapping the variant. Primer binding comparisons use the supplied alignment only.
- KASP candidates are limited to a single SNP and one direction; the special primer is forced on the reference allele. Different Tm for the alternate allele, proprietary tail and primer ratios, additional paralogue/allele effects and actual cluster separation still need to be evaluated independently. Multiple variants within a binding region can prevent this method from producing a valid pair.
- BLAST results are individual local alignments. They may be partial or truncated, are limited by database contents and search sensitivity, and do not automatically account for amplification orientation, spacing or 3′ mismatches. Review with an appropriate independent specificity workflow.
- No claim is made that a length difference can be resolved on a given agarose gel, or that any allele will have a particular HRM melting temperature or fluorescence cluster. Validate against controls.

## Testing

```sh
python -m unittest discover -s tests -v
python -m compileall -q app.py primer_helper tests
```

The automated tests cover sequence parsing, variant positions, primer coordinates, KASP allele selection, plate order and reagent calculations. Tests use controlled Primer3 output rather than a real Primer3 installation. Actual MAFFT, Primer3, BLAST and graphical-desktop runs should be tested on the user's target machine with the installed versions and a small known reference before use.

## References

- [Primer3 manual](https://primer3.org/manual.html)
- [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
- [NCBI BLAST+ manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [LGC: how KASP works](https://www.biosearchtech.com/how-does-kasp-work)
- [LGC: KASP assay FAQ and mix examples](https://www.biosearchtech.com/support/faqs/kasp-genotyping-assays)
- [Plant breeding HRM marker design study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6198901/)
