# Exorcise

Exome-guided reannotation of nucleotide sequences.

Exorcise is described in Lam, S., Thomas, J.C. & Jackson, S.P. Genome-aware
annotation of CRISPR guides validates targets in variant cell lines and enhances
discovery in screens. *Genome Med* **16**, 139 (2024).
<https://doi.org/10.1186/s13073-024-01414-4>.

You can try Exorcise on the web at <https://sjlab.cruk.cam.ac.uk/app/ddrcs/> >
Exorcise. It runs on x86_64 and Apple silicon.

## Introduction

Sequences are usually annotated against one particular genome assembly, and it
is not always right to carry those annotations over to another use. A library of
sgRNAs for CRISPR-Cas9 might have been designed against GRCh37, with each guide
described as targeting some gene or as non-targeting. Searching for those same
guides in GRCh38 will turn up discrepancies, because the assembly has moved on;
and a cell line that does not faithfully reflect GRCh37 will respond to the
guides in ways its annotations do not predict.

Exorcise reannotates sequences according to whether they are present in a genome
and exome that you supply. If that genome and exome reflect the material you are
actually working with, you can be confident about which sequences are present in
it and which annotations hold.

It works by aligning the sequences to the genome with BLAT, then using set
operations between the alignment coordinates and the exon coordinates to transfer
annotations onto the perfect alignments. Each sequence is considered
independently, so a guide that turns out to cut two genes is reported against
both.

Exorcise is written in R and developed by Dr Simon Lam.

## The commands

Exorcise itself reannotates sequences. Four other commands cover the rest of a
screen analysis, and all five are in the Docker image.

| Command             | Does                                             | Language |
|---------------------|--------------------------------------------------|----------|
| `exorcise`          | Reannotate sequences against a genome and exome  | R        |
| `exorcise-trace`    | Plot nucleotide frequency by sequencing cycle    | R        |
| `exorcise-count`    | Count reads and map them to a guide library      | Python   |
| `exorcise-analyse`  | Run MAGeCK, DrugZ or Chronos over a screen       | Python   |
| `exorcise-database` | Build or update a screen results database        | Python   |

A typical run goes `exorcise-trace` to find the guide region in your reads,
`exorcise` to reannotate the library, `exorcise-count` to count against it,
`exorcise-analyse` to score the screen, then `exorcise-database` to collect the
results. `example/examples.txt` works through each step.

The Python commands were previously the separate
[crispr_tools](https://github.com/SimonLammmm/crispr_tools) and
[crispr_screen_viewer](https://github.com/johncthomas/crispr_screen_viewer)
projects. They now live in `py/` as one package.

### Deprecated names

The commands were renamed in 3.0.0. Every old name still works and behaves
identically; each prints a warning naming its replacement.

| Old name                        | Use instead          |
|---------------------------------|----------------------|
| `ntByCycle`                     | `exorcise-trace`     |
| `count_reads`                   | `exorcise-count`     |
| `crispr_pipeline`               | `exorcise-analyse`   |
| `crispr-screen-viewer database` | `exorcise-database`  |

`crispr-screen-viewer` only ever accepted `database` here. Its other
subcommands, and the Dash web viewer, are not part of this project.

## Installation

The Docker image is the recommended way to install everything at once. Conda and
source installations are also supported.

### Prebuilt Docker image (recommended)

Pick a version and architecture from
[Docker Hub](https://hub.docker.com/r/simonlammmm/exorcise/tags), then pull it:

```bash
docker pull simonlammmm/exorcise:<tag>
```

Bind mount a working directory at `/data` and name the command you want:

```bash
docker run --rm -v .:/data/ simonlammmm/exorcise:<tag> exorcise [arguments]
```

Run the image with no command, or with `--help`, to list what it accepts.

With Singularity:

```bash
singularity pull docker://simonlammmm/exorcise:<tag>
singularity run -B .:/data/ exorcise_<tag>.sif exorcise [arguments]
```

Ask your cluster administrator which host directories are mounted by default and
whether you may bind others with `-B`.

### Build the Docker image yourself

The image builds from your working tree, so this is also how you test a change.

```bash
git clone https://github.com/SimonLammmm/exorcise
cd exorcise
docker build -f docker/Dockerfile -t simonlammmm/exorcise .
```

The build context is the repository root, which is why `-f` is needed and why the
trailing `.` matters. To target an architecture other than the host's:

```bash
docker buildx build --platform linux/amd64 -f docker/Dockerfile \
  -t simonlammmm/exorcise:amd64 .
```

### conda

There are two environments, one per language. Install whichever you need.

```bash
git clone https://github.com/SimonLammmm/exorcise
cd exorcise

# For exorcise and exorcise-trace
conda env create -f env/exorcise.yaml
conda activate exorcise
chmod u+x bin/*
export PATH="$PATH:$(pwd)/bin"
exorcise --help

# For exorcise-count, exorcise-analyse and exorcise-database
conda env create -f env/exorcise-py.yaml
conda activate exorcise-py
pip install ./py
exorcise-count --help
```

They are kept separate so that a dependency conflict on one side cannot block the
other. If you want both available at once, activate the Python environment and
add the R environment's `bin/` to your `PATH`, or just use the image.

### From source

1. Clone the repository.
2. For the R commands, install R, the Bioconductor packages, BLAT and twoBitToFa.
   The pinned versions are in `env/exorcise.yaml`.
3. For the Python commands, `pip install ./py`. MAGeCK must be on your `PATH` for
   `exorcise-analyse` to run MAGeCK analyses; Chronos analyses additionally need
   `pip install ./py[chronos]`.
4. Make `bin/*` executable and add `bin/` to your `PATH`.

Exorcise finds its own R modules relative to `bin/`, so it does not matter where
you put the clone. If you need to relocate `bin/` away from `R/`, set
`EXORCISE_HOME` to the repository root. If `blat` and `twoBitToFa` are not on
your `PATH`, point `EXORCISE_BLAT` and `EXORCISE_TWOBITTOFA` at them.

## Syntax

Worked examples are in `example/examples.txt`.

| Long flag (short)       | Required | Value      | Description                                                                            |
|-------------------------|----------|------------|----------------------------------------------------------------------------------------|
| `--infile` (`-i`)       | Yes      | File       | Input file containing the sequences to exorcise.                                       |
| `--outdir` (`-o`)       | Yes      | Directory  | Output directory. Created if it does not exist.                                        |
| `--seq` (`-g`)          | Yes      | Number     | 1-based column number in `--infile` holding the sequences.                              |
| `--genome` (`-v`)       | Yes      | File       | Genome assembly in 2bit format.                                                        |
| `--exome` (`-w`)        | Yes      | File       | Exome annotation in UCSC Table Browser format.                                          |
| `--pam` (`-z`)          | No       | Nucleotide | PAM appended to each sequence before alignment. `[ACGTN]`.                              |
| `--mode` (`-q`)         | No       | String     | CRISPR chemistry. See below. Default `ko`.                                             |
| `--harm` (`-n`)         | No       | Number     | 1-based column number in `--infile` holding the existing annotations.                    |
| `--control` (`-c`)      | No       | Patterns   | Comma-separated regular expressions matching control guides in the `--harm` column.     |
| `--control_type` (`-d`) | No       | Names      | Comma-separated control names, one per `--control` pattern.                             |
| `--expression` (`-x`)   | No       | File       | Two-column file of gene symbol and expression value. Requires `--harm`.                 |
| `--exprcutoff` (`-k`)   | No       | Number     | Expression value below which a gene counts as not expressed. Default `10`.              |
| `--ref`                 | No       | Flag       | Print the citation and exit.                                                            |
| `--version`             | No       | Flag       | Print the version and exit.                                                             |
| `--help` (`-h`)         | No       | Flag       | Print the syntax and exit.                                                              |

Files may be plain or gzipped.

## CRISPR chemistry

Set the chemistry with `--mode`. It decides which features a guide is considered
to act on.

| `--mode`   | Behaviour                                                                                      |
|------------|------------------------------------------------------------------------------------------------|
| `ko`       | CRISPR knockout, the default. Annotates each guide with the features overlapping its cut site. |
| `a`        | CRISPR activation. Annotates each guide with features within 500 bases of its target site.     |
| `i`        | CRISPR interference. Identical to `a`; the distinction is for your records.                    |
| an integer | As `a` and `i`, but with the distance you give instead of 500 bases.                            |
| `cbe`      | Cytosine base editor. C to T within the window [2, 8].                                          |
| `abe`      | Adenine base editor. A to G within the window [4, 9].                                           |
| `beXYnnmm` | Custom base editor. See below.                                                                  |

For a custom base editor, `X` is the original base and `Y` the edited base, both
in `[ACGT]`; `nn` and `mm` are the first and last positions of the editing
window, two digits each, counting from the PAM-distal end of the guide. So
`beCT0208` is equivalent to `cbe`, and `beAG0409` to `abe`.

Base editor modes need an exome that carries CDS coordinates and exon frames,
which means the `name`, `cdsStart`, `cdsEnd` and `exonFrames` columns as well as
the usual ones. The exome bundled in `data/` has only five columns, so download a
fuller refGene table from the UCSC Table Browser for base editing.

## Controls

Guides that fail to align inside any exon are named `exo_Non-targeting_1`,
`exo_Non-targeting_2` and so on. If you know that some of them are declared
controls, pass `--harm` along with `--control` and `--control_type` and they will
be named after the control they belong to instead. For example:

```bash
--harm 2 --control "NO_CURRENT,^NonTargeting" --control_type "Unmapped,NonTargeting"
```

names every unannotated guide whose column 2 value matches `/NO_CURRENT/` as
`Unmapped_1`, `Unmapped_2` and so on. Patterns are tried in the order given, and
each guide takes the first that matches.

## Expression masking

`--expression` takes a two-column file of gene symbol and expression value.
Guides whose `--harm` annotation names a gene scoring below `--exprcutoff` are
dropped before alignment, which is a fast way to restrict a screen to the genes
a particular cell line actually expresses. It requires `--harm`, since that is
the column the symbols are matched against.

## Outputs

Exorcise writes `exorcise.tsv` plus the intermediate files each stage produced.
The intermediates are safe to delete once you are happy with the result.

Each stage is checkpointed: if its output is already present and non-empty, the
stage is skipped. That makes it cheap to rerun with a different `--mode` against
the same genome, since the BLAT alignment is reused. It also means that if you
change `--infile` you must either use a fresh `--outdir` or delete the
intermediates, or the stale alignment will be picked up. Exorcise warns when it
notices alignments that are not in the current input.

| File                                        | Description                                                       |
|---------------------------------------------|-------------------------------------------------------------------|
| `exorcise.1-seq.fa`                         | Sequences plus PAM, as sent to BLAT.                              |
| `exorcise.2-<genome>_BLAT.psl`              | BLAT alignments against the genome.                               |
| `exorcise.3-<genome>_genomicRanges.tsv`     | Coordinates of the perfect alignments, and their cut sites.       |
| `exorcise.3-<genome>_genomicSeqSpecs.tsv`   | Those coordinates, as sent for sequence verification.             |
| `exorcise.3-<genome>_genomicSeqs.fa`        | Sequences read back out of the genome to verify each alignment.    |
| `exorcise.4-<exome>_exonHits.tsv`           | Alignments with the annotations inherited from the exome.          |
| `exorcise.5-<exome>_exonDist.tsv`           | Alignments with the distance to the nearest exon.                  |
| `exorcise.6-<exome>_baseEdited.tsv`         | Predicted base edits. Base editor modes only.                      |
| `exorcise.tsv`                              | The exorcised library.                                            |
| `logfile_exorcise_<timestamp>.log`          | Everything Exorcise reported during the run.                       |

In `exorcise.tsv`, Exorcise's own columns are prefixed `exo_` and come first.
Every column from your input file follows, unchanged. Columns inherited from the
exome are prefixed `inherit.`.

| Column       | Description                                                                 |
|--------------|-----------------------------------------------------------------------------|
| `exo_id`     | Unique identifier for this reannotation.                                    |
| `exo_seq`    | The sequence.                                                               |
| `exo_symbol` | Reannotation from the exome, or a non-targeting or control name.            |
| `exo_orig`   | The original annotation, if `--harm` was given.                             |
| `exo_target` | Genome coordinates of the alignment, including the PAM.                     |
| `exo_cut`    | Genome coordinates of the cut site.                                         |
| `exo_be`     | Symbol and predicted amino acid change. Base editor modes only.             |
| `...`        | Your input columns, then the `inherit.` columns from the exome.             |

One input sequence can produce several rows: a sequence aligning to two loci is
two reannotations, and in base editor mode each predicted protein change is its
own reannotation.

Base editor modes add these columns:

| Column                  | Description                                                           |
|-------------------------|-----------------------------------------------------------------------|
| `be_position_in_genome` | Genome coordinates of the edit.                                       |
| `be_nt_mutation`        | The edit in transcript coordinates, e.g. `C123T`.                     |
| `be_aa_mutation`        | The predicted protein change, e.g. `Q41*`.                            |
| `be_consequence`        | One of `synonymous`, `missense`, `stopgain` or `stoploss`.             |
| `transcript`            | The transcript the prediction was made against.                        |

## exorcise-trace

Plots the proportion of each nucleotide called at each sequencing cycle in a
FASTQ file, which is a quick way to find where the fixed region of an amplicon
ends and the variable guide region begins. Use the answer as `--slice` for
`exorcise-count`.

```bash
exorcise-trace --file reads.fastq.gz --nrows 10000 --out trace-output/
```

| Long flag (short) | Required | Value      | Description                                                      |
|-------------------|----------|------------|------------------------------------------------------------------|
| `--file` (`-f`)   | Yes      | Paths      | Comma-separated FASTQ files or directories to search recursively. |
| `--out` (`-o`)    | No       | Directory  | Output directory. Default the working directory.                  |
| `--nrows` (`-n`)  | No       | Number     | Read only the first N lines of each file. Four lines per read.     |
| `--start` (`-s`)  | No       | Number     | First cycle to plot. Default 1.                                   |
| `--end` (`-e`)    | No       | Number     | Last cycle to plot. Default the longest read.                     |

It writes a PDF, an interactive HTML plot and a TSV of the counts per file.

## exorcise-count

Counts the unique sequences in the sliced region of each read, then maps them
onto a guide library. Output is one dereplicated count file per sample, plus a
single counts table when `--library` is given.

```bash
exorcise-count --slice 0,20 --prefix counts/run \
  --library exorcise-output/exorcise.tsv \
  --seqhdr exo_seq --guidehdr exo_id --genehdr exo_symbol \
  reads/*.fastq.gz
```

| Long flag (short)     | Required | Description                                                                    |
|-----------------------|----------|--------------------------------------------------------------------------------|
| `files`               | Yes      | FASTQ/FASTA files, or directories of them.                                      |
| `--slice` (`-s`)      | Yes      | The guide region of each read, as `M,N` zero-based offsets, end exclusive.       |
| `--prefix` (`-p`)     | Yes      | Prefix for output files. May include a directory.                               |
| `--library`           | No       | Guide library. Also writes `<prefix>.counts.tsv`.                                |
| `--seqhdr` (`-g`)     | No       | Library column with guide sequences. Default `seq`.                             |
| `--guidehdr` (`-j`)   | No       | Library column with guide names. Default `guide`.                               |
| `--genehdr` (`-n`)    | No       | Library column with gene names. Default `gene`.                                 |
| `--allow-mismatch`    | No       | Also count reads one substitution from exactly one library guide.                |
| `--merge-samples`     | No       | Sum counts across files differing only by lane.                                  |
| `--fn-split`          | No       | Substring at which filenames are split to get the sample name. Default `_R1_`.    |
| `--suffix`            | No       | Suffix for per-sample count files. Default `.rawcount`.                          |
| `--file-type`         | No       | `a`, `q` or `infer`. Default `infer`.                                            |
| `--overwrite`         | No       | Recount samples whose output already exists.                                     |

An Exorcise output makes a good library here: pass `--seqhdr exo_seq
--guidehdr exo_id --genehdr exo_symbol` and your counts arrive already carrying
genome-aware annotations.

Guides shorter than 16 bases are not counted, because below that a sequence stops
being specific enough to identify one guide. If the library's guides vary in
length, or are shorter than the slice, a window is slid along the sliced region.

## exorcise-analyse

Runs MAGeCK, DrugZ or Chronos over a screen. What runs is decided by the Analyses
sheet of an `.xlsx` workbook; start from
`example/example-analysis-workbook-template.xlsx`. A `.json` configuration
carrying `sample_reps`, `analyses` and `control_groups` also works.

```bash
exorcise-analyse --counts counts/ --output-dir results/ my-screen.xlsx
```

| Long flag             | Required | Description                                                                 |
|-----------------------|----------|-----------------------------------------------------------------------------|
| `config_file`         | Yes      | An `.xlsx` workbook or `.json` configuration.                                |
| `--counts`            | Yes      | Counts file, or the directory holding the counts files the workbook names.    |
| `--output-dir`        | No       | Where the results directory is created. Default the working directory.         |
| `--file-prefix`       | No       | Prefix for generated files. Default `result`.                                 |
| `--skip-method`       | No       | Comma-separated methods to skip.                                             |
| `--run-groups`        | No       | Comma-separated control groups to include. All by default.                     |
| `--run-analyses`      | No       | Comma-separated analysis names to include. All by default.                     |
| `--analysis-version`  | No       | Results are filed under this, if set.                                        |
| `--overwrite`         | No       | Rerun analyses whose output already exists.                                   |
| `--dry-run`           | No       | Validate the configuration without running anything.                          |

Every analysis is checkpointed on its output file, so a rerun only does what is
missing. `--dry-run` is worth using first: it reports every problem it can find
in the workbook at once.

Chronos only applies to comparisons that span more than one timepoint with the
control at the earliest one, so it needs `Days grown` and `Trajectory` columns in
Sample details. Comparisons that do not qualify are skipped with a reason.

## exorcise-database

Collects the tables `exorcise-analyse` wrote into a queryable database, alongside
two gzipped CSVs of experiment and comparison metadata.

```bash
exorcise-database --out-dir db/ --results-dir results/ --counts-dir counts/ \
  --new-db my-screen.xlsx
```

| Long flag (short)        | Required | Description                                                        |
|--------------------------|----------|--------------------------------------------------------------------|
| `details_xlsx`           | Yes      | One or more analysis workbooks.                                     |
| `--out-dir` (`-o`)       | Yes      | Where the database and metadata are written.                        |
| `--results-dir` (`-r`)   | No       | Directory holding the analysis results. Default `results`.           |
| `--counts-dir` (`-c`)    | No       | Directory holding the counts files. Default `counts`.                |
| `--filename-prefix` (`-p`)| No      | Prefix the analysis files were written with. Default `result`.        |
| `--new-db` (`-n`)        | No       | Create a new database rather than adding to one.                     |
| `--update-existing` (`-u`)| No      | Replace experiments already present instead of skipping them.         |
| `--force-overwrite` (`-f`)| No      | With `--new-db`, replace existing files without asking.               |
| `--verbosity` (`-v`)     | No       | 0 warnings, 1 info, 2 debug. Default 1.                              |

It writes `database.db`, `experiments_metadata.csv.gz` and
`comparisons_metadata.csv.gz`. Nothing here touches the network.

`--new-db` over a database that already exists asks before replacing it. Under
`docker run` there is no terminal to answer with unless you pass `-it`, so in a
script or a scheduled job add `--force-overwrite`. Without either, it refuses and
leaves the database alone rather than proceeding. Omit `--new-db` altogether to
add to the existing database instead of replacing it.

## Repository layout

| Path          | Contents                                                                        |
|---------------|---------------------------------------------------------------------------------|
| `bin/`        | The executables, including shims for the deprecated names.                        |
| `R/`          | The Exorcise pipeline, sourced by `bin/exorcise.R`.                              |
| `py/`         | `exorcise_screens`: counting, analysis and the database builder.                   |
| `env/`        | The two conda environments, with pinned versions.                                 |
| `docker/`     | Dockerfile and entrypoint.                                                       |
| `data/`       | A bundled human exome annotation.                                                |
| `example/`    | Example inputs, a workbook template, and commands for every step.                  |
| `manuscript/` | Analysis scripts for the 2024 paper. Pinned; not maintained alongside the rest.    |

## About

Exorcise is developed and maintained by Dr Simon Lam, University of Cambridge,
<sl681@cam.ac.uk>.

`exorcise-count`, `exorcise-analyse` and `exorcise-database` derive from
crispr_tools and crispr_screen_viewer by John C. Thomas. `exorcise-analyse`
bundles [DrugZ](https://github.com/hart-lab/drugz) by Medina Colic and Traver
Hart, whose numerical routines are reproduced unchanged.

## Licence

This code is licensed under the Creative Commons Zero v1.0 Universal licence.

Exorcise uses [BLAT](https://kentinformatics.com/#BLAT), which is free for
academic, personal and non-profit use. `exorcise-analyse` calls
[MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) and, optionally,
[Chronos](https://github.com/broadinstitute/chronos); both carry their own
licences.

## Citation

If you use Exorcise in your work, please cite Lam, S., Thomas, J.C. & Jackson,
S.P. Genome-aware annotation of CRISPR guides validates targets in variant cell
lines and enhances discovery in screens. *Genome Med* **16**, 139 (2024).
<https://doi.org/10.1186/s13073-024-01414-4>.
