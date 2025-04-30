# Exorcise
Exome-guided reannotation of nucleotide sequences (Exorcise) is available at:
https://github.com/SimonLammmm/exorcise/

Exorcise is described in Lam, S., Thomas, J.C. & Jackson, S.P. Genome-aware annotation of CRISPR guides validates targets in variant cell lines and enhances discovery in screens. Genome Med 16, 139 (2024). https://doi.org/10.1186/s13073-024-01414-4.

You can try Exorcise on the web at https://sjlab.cruk.cam.ac.uk/app/ddrcs/ > Exorcise.

Exorcise is easily installed using Docker and runs on x86_64 and Apple silicon architectures.

This document describes the concept, installation, basic usage, and syntax of Exorcise.

## Introduction

Sequences are often prescribed with annotations based on a certain genome assembly. It is not always appropriate to accept these annotations for all use cases of those sequences. For instance, a library of sgRNA sequences for CRISPR-Cas9 may be designed based on GRCh37 and sgRNAs annotated as "targeting" a certain gene or "non-targeting". Sequence search of those sgRNAs in GRCh38 might reveal discrepancies because of updates to the assembly; further, cell lines which do not perfectly reflect GRCh37 will show differential reactivity to the CRISPR guides as compared their annotations. Exorcise reannotates sequences based on their presence or absence in a user-supplied genome and exome. If the supplied genome and exome reflect the subject under consideration, then the user can be confident as to the presence or absence of sequences in the subject and the validity of the annotations.

Exorcise aligns sequences to the genome and implements set mathematics operations between alignment targets and exon coordinates to transfer annotations to the perfect alignments. It is written entirely in R and relies on BLAT for perfect sequence alignment.

Exorcise calculates sequence-level reannotations whereby each sequence is considered in turn. If provided a grouping variable (for example, original annotations), Exorcise also calculates group-level reannotations where a consensus sequence-level reannotation is determined for each group and applied to all sequences in the group. We refer to this behaviour as harmonisation.

Exorcise is developed and maintained by Dr Simon Lam and is available at https://github.com/SimonLammmm/exorcise/.

## Installation

Docker installation of Exorcise is recommended. You can build your own Docker image for your operating system or pull a prebuilt one from Docker Hub to install Exorcise in one command.

You can also install Exorcise using conda or from source.

### Prebuilt Docker image (recommended)

Navigate to the Exorcise repo on [Docker Hub](https://hub.docker.com/r/simonlammmm/exorcise/tags) and choose the desired version and OS. Then, pull:
```
docker pull simonlammmm/exorcise:<tag>
```

For example, to pull version 1.5.3 for Apple silicon:
```
docker pull simonlammmm/exorcise:1.5.3_arm64
```

To pull Exorcise with Singularity:
```
singularity pull docker://simonlammmm/exorcise:<tag>
```

Run. At runtime, bind mount a directory to the `/data` location in the Docker container.
```
docker run --rm -v .:/data/ simonlammmm/exorcise exorcise [arguments]
```
For Singularity, speak to your cluster administrator to see what host directories are mounted in the image by default and if you are allowed to bind other directories into the container with `-B`. Run with the command:
```
singularity run -B .:/data/ exorcise_<tag>.sif exorcise [arguments]
```

The Docker and Singularity installations also include [crispr_tools](https://github.com/SimonLammmm/crispr_tools) and [crispr_screen_viewer](https://github.com/johncthomas/crispr_screen_viewer). Possible commands are:
* `docker run --rm -v .:/data/ simonlammmm/exorcise exorcise [arguments]`
* `docker run --rm -v .:/data/ simonlammmm/exorcise ntByCycle [arguments]`
* `docker run --rm -v .:/data/ simonlammmm/exorcise count_reads [arguments]`
* `docker run --rm -v .:/data/ simonlammmm/exorcise crispr_pipeline [arguments]`
* `docker run --rm -v .:/data/ simonlammmm/exorcise crispr-screen-viewer [arguments]`

### Build from Dockerfile

To build your own Docker image with Exorcise, follow the below steps. You might need to do this if your desired version or OS are not available on Docker Hub.

1. Clone this repo and `cd` into it.
```
git clone https://github.com/SimonLammmm/exorcise
```
2. Build from the Dockerfile.
```
docker build -t simonlammmm/exorcise docker/
```
3. Run, as above for the prebuilt Docker images.
```
docker run --rm -v .:/data/ simonlammmm/exorcise exorcise [arguments]
```

### conda

1. Clone this repo and `cd` into it.
```
git clone https://github.com/SimonLammmm/exorcise
```
2. Install dependencies using the conda environment yaml file.
```
conda env create -f env/exorcise.yaml
```
3. Add the executables in `bin/` to your `PATH` using the appropriate method for your system and shell.
4. Make the executables executable.
```
chmod u+x bin/*
```
5. Run by calling `exorcise` on the command line.
```
exorcise
```

### Install from source

1. Clone this repo and `cd` into it.
```
git clone https://github.com/SimonLammmm/exorcise
```
2. Install R, Bioconductor, and BLAT dependencies listed in `env/exorcise.yaml` using the appropriate method for your system.
3. Add the executables in `bin/` to your `PATH` using the appropriate method for your system and shell.
4. Make the executables executable.
```
chmod u+x bin/*
```
5. Run by calling `exorcise` on the command line.
```
exorcise
```

## Syntax

Example commands for the Docker image can be found in the `example/` folder.

Exorcise takes mandatory and optional arguments.

| Long flag (short flag) | Mandatory?                     | Value      | Description                                                                                                                      |
|------------------------|--------------------------------|------------|----------------------------------------------------------------------------------------------------------------------------------|
| --infile (-i)          | Yes                            | File       | Input file with sequences for Exorcise.                                                                                          |
| --outdir (-o)          | Yes                            | Directory  | Output directory.                                                                                                                |
| --seq (-g)             | Yes                            | Number     | Column number in the infile with sequences, 1-based integer.                                                                     |
| --pam (-z)             | No                             | Nucleotide | PAM sequence. [ACGTN] supported.                                                                                                 |
| --library (-l)         | Yes, if -v -w not specified    | File       | Existing Exorcise library, if using post-hoc mode (see Modes section).                                                           |
| --mode (-q)            | No                             | String     | CRISPR chemistry: ko (knockout, default), i (inhibition), a (activation), cbe (cytosine base editor), abe (adenine base editor). |
| --genome (-v)          | Yes, if -l not specified       | File       | Genome in 2bit format.                                                                                                           |
| --exome (-w)           | Yes, if -l not specified       | File       | Exome from UCSC Table Browser.                                                                                                   |
| --harm (-n)            | No                             | Number     | Column number in the infile with prior symbols, 1-based integer. Enables harmonisation (see Modes).                              |
| --priorities (-y)      | Yes, if -n specified and not 0 | File       | Feature priorities list in NCBI Datasets format. Enables harmonisation (see Modes).                                              |
| --control (-c)         | No                             | String     | Regex strings in the -n column to treat as controls, comma-separated list. Enables control reannotation (see Modes).             |
| --control_type (-d)    | No                             | String     | Control types, comma-separated list, same length as -c.                                                                          |
| --help (-h)            | No                             | Flag       | Show help and then exit.                                                                                                         |

## CRISPR chemistry

Use the `-q` flag to set the CRISPR chemistry.

* Use `-q ko` (default) for CRISPR knockout. Guides are annotated with features overlapping the Cas9 cut site.
* Use `-q a` for CRISPR activation and `-q i` for CRISPR interference. In both of these modes, guides are annotated with features within 500 bp of the guide target site. The distinction in the syntax is for documentation purposes only.
* Use `-q NUM` where `NUM` is an integer for CRISPR activation/interference chemistry considering features within `NUM` bp of the guide target site.
* Use `-q cbe` for cytosine base editor. This calculates C -> T base edits within the base editing window [2,8] relative to the guide target site.
* Use `-q abe` for adenine base editor. This calculates A -> G base edits within the base editing window [4,9] relative to the guide target site.
* Use `-q beXY0123` for custom base editor chemistry, where `X` is the original base, `Y` is the new base, `01` is the start of the base editing window, and `23` is the end of the base editing window. `X` and `Y` must be in [ACGT]; and `01` and `23` must be two digits each, 1-based indexing from the PAM-distal end of the guide.

## Modes

Exorcise can use previous Exorcise outputs (Exorcise libraries) to reannotate sequences without recalculation of the entire sequence vector each time. This is called post-hoc mode and is extremely fast. Post-hoc mode takes precedence over ad-hoc mode (see below).

When an Exorcise library is not specified, then it is calculated from the genome and exome. This is called ad-hoc mode. Ad-hoc mode can take a long time (from 20 minutes to over a week, depending on inputs) to run BLAT. The output can be used as an Exorcise library the next time that sequences need to be reannotated in the same genome/exome context. 

Both post-hoc and ad-hoc modes support harmonisation on the harm column if specified along with a feature priorities file. The harmonisation algorithm can take a long time if there are many groups in the harm column.

The following table shows Exorcise behaviours when different sets of inputs are received. Combined behaviour is allowed by specifying the necessary combinations of inputs. (Note that the inputs -i, -o, and -g are mandatory and so are omitted from the table.) 

| Inputs                   | Behaviours                                                                        |
|--------------------------|-----------------------------------------------------------------------------------|
| `-v`; `-w`               | Ad-hoc mode: aligns guides to the genome and uses set theory to find exonic hits. |
| `-l`                     | Post-hoc mode: uses an existing Exorcise library to reannotate input sequences.   |
| `-n`; `-y`               | Harmonisation: reannotates prior symbols.                                         |
| `-n`; `-y`; `-c`, (`-d`) | Control reannotation: reannotates controls to the optional argument `-d`.         |
| `-l`; omit `-n`, `-y`    | Use post-hoc mode and inherit existing harmonisations from the library, if any.   |
| `-l`; `-n 0`             | Use post-hoc mode and explicitly specify not to inherit existing harmonisations.  |

## Outputs

Exorcise outputs a single exorcise.tsv file as well as intermediate files for diagnostic purposes. If you are happy with the results, then the intermediate files can be removed.

Exorcise uses checkpointing. It checks for the presence of an intermediate file and uses it rather than regenerating it if it exists from a previous run.

|     File                                          |     Description                                                                   |
|---------------------------------------------------|-----------------------------------------------------------------------------------|
|     exorcise.1-seq.fa                             |     Sequences plus PAM to   be sent to BLAT.                                      |
|     exorcise.2-genome.2bit_BLAT.psl               |     BLAT genome alignment   results.                                              |
|     exorcise.3-genome.2bit_genomicRanges.tsv      |     Coordinates of BLAT   alignments and calculated cut sites.                    |
|     exorcise.3-genome.2bit_genomicSeqSpecs.tsv    |     Coordinates of BLAT   alignments to be sent for sequence validation.          |
|     exorcise.3-genome.2bit_genomicSeqs.fa         |     Validated sequences.                                                          |
|     exorcise.4-exome.gz_exonHits.tsv              |     BLAT alignments with   annotations inherited from the exome.                  |
|     exorcise.5-exome.gz_exonDist.tsv              |     BLAT alignments and   exome annotations with distance to the nearest exon.    |
|     exorcise.tsv                                  |     Main Exorcise output   file.                                                  |

The main Exorcise output file is exorcise.tsv. Exorcise outputs are in columns starting with exo_. Arbitrary columns included in the input file are appended after the Exorcise columns.

|     Column        |     Description                                                                    |
|-------------------|------------------------------------------------------------------------------------|
|     exo_id        |     Unique identifier for   this Exorcise reannotation.                            |
|     exo_id_harm   |     Unique identifier for   this Exorcise harmonisation.                           |
|     exo_seq       |     Sequence.                                                                      |
|     exo_symbol    |     Sequence-level   reannotation from exome.                                      |
|     exo_harm      |     Group-level annotation   (harmonisation) from exome.                           |
|     exo_orig      |     Original annotation.                                                           |
|     exo_target    |     Genome coordinates of   alignment between sequence plus PAM and the genome.    |
|     exo_cut       |     Genome coordinates of   the cut site.                                          |
|     ...           |     Arbitrary columns   included in the input file.                                |

## About
Exorcise is developed and maintained by Dr Simon Lam, University of Cambridge, mailto:sl681@cam.ac.uk.

## License
This code is licensed under the Creative Commons Zero v1.0 Universal license.

Ad-hoc Exorcise uses [BLAT](https://kentinformatics.com/#BLAT), which is available free for academic, personal, and non-profit uses.

## Citation
If you use Exorcise in your work, please cite Lam, S., Thomas, J.C. & Jackson, S.P. Genome-aware annotation of CRISPR guides validates targets in variant cell lines and enhances discovery in screens. Genome Med 16, 139 (2024). https://doi.org/10.1186/s13073-024-01414-4.





