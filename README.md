# Exorcise
Exome-guided reannotation of nucleotide sequences (Exorcise) is available at:
https://github.com/SimonLammmm/exorcise/

This document describes the basic usage and syntax of Exorcise.

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

To install Exorcise with Singularity:
```
singularity pull docker://simonlammmm/exorcise:<tag>
```

Run. At runtime, bind mount a directory to the `/data` location in the Docker container.
```
docker run --rm -v .:/data/ simonlammmm/exorcise exorcise [arguments]
```
For Singularity, speak to your cluster administrator to see what host directories are mounted in the image by default.
```
singularity run exorcise_<tag>.sif exorcise [arguments]
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

Exorcise takes mandatory and optional arguments.

|     Long flag (short flag)    |     Mandatory?                         |     Value         |     Description                                                                                                                                                                                               |
|-------------------------------|----------------------------------------|-------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|     --infile   (-i)           |     Yes                                |     File          |     File containing   sequences to be Exorcised. Arbitrary columns will be returned in the output   file.                                                                                                     |
|     --outdir   (-o)           |     Yes                                |     Directory     |     Working directory to   place intermediate files and the output file exorcise.tsv.                                                                                                                         |
|     --seq   (-g)              |     Yes                                |     Number        |     Column number in the infile   that corresponds to sequences to be Exorcised. 1-based integer, must not be   greater than the number of columns in the infile.                                             |
|     --pam   (-z)              |     No                                 |     Nucleotide    |     Sequence to be appended   to the 5' end of each sequence for BLAT purposes only. Must be in the   single-letter nucleotide alphabet (ACGTN).                                                              |
|     --library   (-l)          |     Yes, if -v   -w   not specified    |     File          |     Exorcised library to be   used for the re-annotation. If this is not specified, then genome   and exome   must be specified. Must be an Exorcise library. Enables post-hoc mode (see   Modes section).    |
|     --mode   (-q)             |     No                                 |     String        |     CRISPR chemistry. Can be any  of `ko` (knockout, default), `i` (inhibition), or `a` (activation), or a number (arbitrary).                                                                                                         |
|     --genome   (-v)           |     Yes, if -l   not specified         |     File          |     Genome in 2bit format. Enables   ad-hoc mode (see Modes section). Ignored if library   is specified.                                                                                                      |
|     --exome   (-w)            |     Yes, if -l   not specified         |     File          |     Exome from UCSC Table   Browser. Enables ad-hoc mode (see Modes section). Ignored if library   is specified.                                                                                              |
|     --harm   (-n)             |     No                                 |     Number        |     Column number in the infile   containing groups. 1-based integer, must not be greater than the number of   columns in the infile. Enables harmonisation   (see Modes).                                    |
|     --priorities   (-y)       |     Yes, if -n   specified and not 0   |     File          |     Feature priorities list in NCBI Datasets format. Ignored if harm is not specified. Enables   harmonisation (see Modes).                                                                                                         |
|     --control   (-c)          |     No                                 |     String        |     Comma-separated list of   strings to treat as control sequences. R regex allowed. Ignored if harm   is not specified. Enables control reannotation (see Modes).                                                       |
|     --control_type   (-d)     |     No                                 |     String        |     Comma-separated list of   strings of the same length as control   indicating control type. Ignored if control   is not specified. Ignored if harm is not specified.                                       |
|     --help   (-h)             |     No                                 |     Flag          |     Show help and then   exit.                                                                                                                                                                                |

## Modes

Exorcise can use previous Exorcise outputs (Exorcise libraries) to reannotate sequences without recalculation of the entire sequence vector each time. This is called post-hoc mode and is extremely fast. Post-hoc mode takes precedence over ad-hoc mode (see below).

When an Exorcise library is not specified, then it is calculated from the genome and exome. This is called ad-hoc mode. Ad-hoc mode can take a long time (from 20 minutes to over a week, depending on inputs) to run BLAT. The output can be used as an Exorcise library the next time that sequences need to be reannotated in the same genome/exome context. 

Both post-hoc and ad-hoc modes support harmonisation on the harm column if specified along with a feature priorities file. The harmonisation algorithm can take a long time if there are many groups in the harm column.

The following table shows Exorcise behaviours when different sets of inputs are received. Combined behaviour is allowed by specifying the necessary combinations of inputs. (Note that the inputs -i, -o, and -g are mandatory and so are omitted from the table.) 

|     Inputs                 |     Behaviours                                                                                |
|----------------------------|-----------------------------------------------------------------------------------------------|
|     -v, -w                 |     Ad-hoc mode: aligns   guides to the genome and uses set theory to find exonic hits.       |
|     -l                     |     Post-hoc mode: uses an   existing Exorcise library to reannotate input sequences.         |
|     -n, -y                 |     Harmonisation: also performs   group-level reannotation.                                  |
|     -n, -y,   -c   (-d)    |     Control reannotation: also   reannotates controls to the optional vector control_type.    |
|     -l, omit -n,   -y      |     Use post-hoc mode and inherit   existing harmonisations from the library, if any.         |
|     -l, -n 0               |     Use post-hoc mode and   explicitly specify not to inherit existing harmonisations.        |

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
Exorcise was developed by Dr Simon Lam, University of Cambridge, mailto:sl681@cam.ac.uk.

## License
This code is licensed under the Creative Commons Zero v1.0 Universal license.








