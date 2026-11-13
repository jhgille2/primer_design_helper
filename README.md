# Primer Design Helper

A local desktop application for developing PCR genotyping assays from candidate genomic regions. It is oriented toward a *Glycine max* × *Glycine soja* project, but accepts other DNA sequences. The app keeps sequence data on your computer and never submits it to external services.

> **Research prototype.** Review the underlying alignments, genomic specificity and experimental controls. This software does not establish orthology or prove that an assay is locus-specific. Candidate primers and genotype calls require independent and experimental validation before use.

## What it does

1. Downloads and indexes selected public soybean reference genomes, helps review candidate homologous regions, and accepts two or more genomic FASTA sequences, calls locally installed MAFFT for alignment or accepts aligned FASTA, and finds differences between the first and second records. Additional records participate in primer-binding checks. The supplied example is synthetic.
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

  Enter the **database prefix** `soybean_db` or its full path in the app, not the name of a particular `.nhr`, `.nin`, or `.nsq` file. Search each parent/reference assembly that is relevant; record assembly names and versions. The Soybean library can also download and index curated genomes as described below.

External executables and databases are independent installations and are not included in this repository. An unavailable executable gives an explanatory error. Long-running tools execute in a background thread.

## Detect and install software in the app

Open **Software setup** (or the **Software** tab). The app checks `mafft`, `primer3_core`, `blastn`, `makeblastdb`, and `blastdbcmd` on startup; click **Check software** to check again. Select a row and use **Browse for selected executable** if a program is installed but not on your PATH. BLAST+ supplies the last three executables. The check runs each program; Primer3 is checked using a small Boulder-IO request. A successful response does not validate a specific analysis.

On **Linux and macOS**, if Conda or Mamba is already installed, select **Install missing tools**. The app installs missing Bioconda/conda-forge packages into `data/tools/bioinformatics`, a dedicated environment under the repository. It does not change the base environment. It then checks the resulting executables and fills their app paths; subsequent launches discover that managed environment. Internet access, free disk space, package availability for your architecture, and a working Conda/Mamba installation are required. If Conda is not on PATH, enter or browse to its executable; the app links to the official Miniforge instructions. The managed tool directory is excluded from Git. If the target directory already exists but is not a valid Conda environment, the app refuses to overwrite it.

**Native Windows:** Bioconda does not publish Windows packages. Automatic Bioconda installation is therefore unavailable in a native Windows session. The same tab links to the official MAFFT, Primer3 and NCBI BLAST+ downloads; follow each project's Windows instructions and browse to the executables after installation. Alternatively, run the *entire application*, Python, and dependencies inside Linux/WSL with a working graphical display. An executable in WSL is not automatically usable by native Windows Python. The app does not install WSL or Conda for you. Some tools/platform combinations may require manual installation; availability varies.

## Download a soybean reference library

The **Soybean library** tab downloads selected public *genome assemblies* from the shared [SoyBase / Legume Information System datastore](https://www.legumeinfo.org/download/glycine/). The code repository does not contain the genomes. The catalog contains:

| Species | Accession | Assembly | Catalog key | Compressed genome, approx. |
| --- | --- | --- | --- | ---: |
| *G. max* | Williams 82 | Wm82.gnm6.S97D | `wm82_v6` | 255 MB |
| *G. max* | Williams 82 | Wm82.gnm2.DTC4 | `wm82_v2` | 287 MB |
| *G. max* | Zhonghuang 13 | Zh13.gnm1.N6C8 | `zh13_v1` | 294 MB |
| *G. soja* | W05 | W05.gnm1.SVL1 | `w05_v1` | 296 MB |
| *G. soja* | PI 483463 | PI483463.gnm1.YJWS | `pi483463_v1` | 292 MB |
| *G. soja* | PI 562565 (SoyW01) | PI_562565.gnm1.Z2VG | `soyw01_v1` | 261 MB |

These are full genome assemblies, not the 1,110-accession variant pangenome or all assemblies available at SoyBase. Wm82 v2 and v6 are two assemblies of the same accession, not two biological accessions. W05 is a reference for that accession; it is not automatically the sequence of a different *G. soja* parent. Source URLs, assembly names and published checksums are recorded locally. Review the source collection pages for publication/citation and use conditions.

### Set up from the desktop app

1. Install **NCBI BLAST+**, including `makeblastdb`, `blastn` and `blastdbcmd`. You also need MAFFT to align extracted regions. Confirm they work at a terminal or enter full executable paths in the app. Primer3 is separate.
2. Open **Soybean library**. The default library folder is `data/soybean/` in the cloned repository; choose another folder if preferred. By default Williams 82 v6 and W05 are highlighted. Use Ctrl/Shift (Command on some macOS configurations) to select the assemblies you want.
3. Click **Download and build selected**. This performs a large download for each assembly, retrieves the source's MD5 checksum, verifies the compressed data, decompresses it, and creates a local nucleotide BLAST database. Watch the status at the bottom of the app. Downloads run one at a time and can take a while. Plan for several GB free per assembly for the compressed file, uncompressed genome and BLAST index, plus temporary working space.
4. Wait until the assembly status shows **ready**. If the network, checksum or BLAST build fails, read the message and try again. Completed compressed downloads are reused after checksum verification; interrupted partial downloads are removed. If an old incomplete `blast` directory is reported, move that directory aside and rebuild. Do not delete a reference that another process is using.

The `data/soybean/` default is in `.gitignore`, so the genomes and indexes will not be pushed to GitHub. If you choose another folder within the repository, add that path to `.gitignore` yourself. The files are not created until you explicitly download them. A source checksum verifies accidental transfer errors, not independent authenticity.

### Start from the command line instead

From the repository directory:

```sh
python -m primer_helper.soybean --list
python -m primer_helper.soybean --install wm82_v6 w05_v1
```

Select any catalog keys shown by `--list`. For another storage location or a non-default executable, add `--directory /path/to/soybean_data` and/or `--makeblastdb /path/to/makeblastdb`. On Windows, use `py` if that is how Python is installed. Restart the app or use **Refresh status** after a command-line installation.

### Search and extract candidate loci

1. In **Soybean library**, select the **ready** assemblies you want to search. Paste one genomic query sequence, 30–50,000 bases, with or without a FASTA header. Query gaps are removed. For a target variant, a unique genomic region covering the locus with useful flanks is generally more informative than a short coding fragment; high-complexity 300–1,000 bp or longer can be a starting point, depending on the region.
2. Click **Search locus in selected**. The app uses nucleotide `megablast` against each selected whole-genome database and displays individual local HSPs, their query and subject spans, percent identity, strand, and bit score. The search is limited to 30 subject sequences per database and up to three HSPs per subject; it is a candidate discovery tool, not an exhaustive paralogue screen or a whole-genome alignment.
3. Inspect hits carefully using position, coverage, annotation and genomic context in SoyBase. Similarity or the top score does **not** establish orthology, and multiple rows can be different HSPs from the same locus. Select two or more hits; their displayed order determines the order in the next tab, with the first as reference. Choose extra flanking bases and click **Use selected hits for alignment**. The app extracts each displayed HSP's subject span plus flanks and orients it to the query. An HSP might cover only part of a gene or query, so select and review accordingly. Subject positions are 1-based inclusive, with flanks clipped to contig boundaries. IUPAC ambiguity codes other than `N` are converted to `N`; their count is recorded in each FASTA header.
4. In **Sequences**, align the extracted regions with MAFFT, inspect the alignment, then **Find variants**. Verify both parental genotypes and binding sites using the actual parents; the downloaded genomes are reference samples.
5. For primer specificity, paste or copy a candidate pair into **Specificity**, return to the library tab, select the ready assemblies and click **Search primers in selected**. Results are individual `blastn-short` binding hits for each assembly. Review pairing, orientation, product sizes and 3′ mismatches with an appropriate in-silico PCR method; a lack of hits is not proof of uniqueness. This is not an automated genome-wide primer-pair specificity assessment.

No sequence or query is uploaded by the search: after initial downloads, genomic searches and extraction run locally. The app does not yet search by gene ID, download gene annotations, look up orthology groups, or automatically transfer coordinates across assemblies. Obtain a query region and its metadata from SoyBase or your own data, and review candidate hits before alignment.

### POWR1 example using actual reference assemblies

1. In **Soybean library**, select **Williams 82 Wm82.gnm2** and **W05 W05.gnm1**, then click **Download and build selected**. Wm82.gnm2 is required here because the published coordinates below belong to that assembly, not Wm82.gnm6. The downloads are large.
2. Click **Load POWR1 query**. The app reads the installed Wm82.gnm2 chromosome 20 FASTA and extracts the published POWR1 gene, Glyma.20G085100, plus 500 bp on either side: Gm20:31,774,270–31,780,304, 1-based inclusive, on the assembly's plus strand. The FASTA header records the assembly, chromosome, coordinates, strand and ambiguity conversions. The query comes from the genome; it is not a synthetic sequence. This action does not itself load W05.
3. Click **Search locus in selected** with both assemblies selected. Check candidate hits in both references, including their query coverage, location, orientation and neighboring genes; W05 orthology is not established by a top BLAST score. Multiple rows can be separate partial alignments of the same locus. Select the Wm82 and candidate W05 rows and choose **Use selected hits for alignment**. Check that extracted regions cover the variant and adequate flanks; partial hits and repetitive sequence may complicate this. Align with MAFFT, inspect the result and identify the indel. The selected FASTA headers document the extracted intervals. Save the project to retain the sequences.
4. A published study describes a 321 bp difference comprising a 304 bp transposable element and a 17 bp target-site duplication; another comparison reports a 318 bp W05 deletion. Variant calls and representations differ, including one analysis that did not call a W05 deletion. Do not assume the same gap size or allele in every dataset: check this actual alignment and its genomic context. Design and validate against the actual parents and relevant paralogues. See [Goettel et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9160092/) and [Marsh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC9033719/); the [Wm82.a2 gene coordinates](https://pmc.ncbi.nlm.nih.gov/articles/PMC9372550/) are specific to that assembly.

The existing startup example remains synthetic; **Load POWR1 query** prepares a real genomic query after reference download. No pre-extracted W05 sequence or verified W05 allele is bundled. A successful local extraction/search is required to produce the comparison.


## A first assay

1. Identify the actual parental lines, locus, variant and intended readout. Obtain genomic sequences for orthologous regions with substantial flanking sequence. Confirm assembly/version, strand and coordinates. A public reference is a proxy unless it represents your actual parent.
2. On **Sequences & variants**, open or paste FASTA. Put the design/reference parent first and comparison parent second. Align with MAFFT or provide equal-length aligned sequences. The built-in example is already aligned: use **Find variants** directly, or choose **Align or realign with MAFFT** to start over. For a fresh alignment the app strips existing `-` gap characters, retains `N` for unknown bases, and MAFFT may place new gaps differently. A gap is an alignment placeholder; it is not itself an unknown nucleotide. Inspect orientation, gaps, N bases, repeats, gene structure, genomic neighborhood and orthology independently. Press **Find variants** and select one. Alignment positions shown are one-based; the displayed reference offset is zero-based, and for insertions marks a boundary. Consecutive differences are grouped and an N is not called as a variant.
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

The automated tests cover local library setup with mocked downloads and BLAST, checksum handling, strand-aware sequence extraction, sequence parsing, variant positions, primer coordinates, KASP allele selection, plate order and reagent calculations. Tests use controlled tool output rather than real external installations. Actual MAFFT, Primer3, BLAST and graphical-desktop runs should be tested on the user's target machine with the installed versions and a small known reference before use.

## References

- [Primer3 manual](https://primer3.org/manual.html)
- [NCBI Primer-BLAST](https://www.ncbi.nlm.nih.gov/tools/primer-blast/)
- [NCBI BLAST+ manual](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
- [LGC: how KASP works](https://www.biosearchtech.com/how-does-kasp-work)
- [LGC: KASP assay FAQ and mix examples](https://www.biosearchtech.com/support/faqs/kasp-genotyping-assays)
- [Plant breeding HRM marker design study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6198901/)
