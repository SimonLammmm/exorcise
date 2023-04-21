# exorcise
EXOme-guided Reannotation of nuCleotIde SEquences
or exome-guided reannotation of CRISPR sequences, as it was originally designed

This software was previously codenamed crispieR (CRISPR Identification Enhancement in R) and filenames and references still reflect this old name.

## Description
exorcise accepts a vector of nucleotide sequences, aligns them to a supplied reference genome, and annotates them according to a supplied reference feature set based on co-ordinates.

## Background
Annotations to standard reference assemblies such as GRCh38 and feature sets such as RefSeq are not universally applicable to all use cases. Experimental cell lines are highly likely to differ from references in at least a small but detectable way and these discrepancies influence, for example, targeting efficiency of CRISPR-Cas9 or RNA interference guide RNAs. This can be overcome by design of guide RNAs to arbitrary genomes and feature sets that are specific to the target experimental cell line.

exorcise is a software suite that verifies the presence or absence of nucleotide sequences in a supplied reference genome and provides annotations according to a supplied reference feature set according to a consistent genomic co-ordinate system. This software can be used to inform whether pre-constructed libraries of nucleotides exist in a supplied genome/exome and can reannotate sequencing counts data to the correct targets for the supplied genome/exome.

## Algorithm
exorcise requires a vector of short nucleotide sequences as an input as well as a genome assembly in `.2bit` format and a UCSC-style feature set file. Examples of these files can be found in `tutorial/tutorial/data`. Nucleotide sequences are perfectly aligned to the genome using BLAT, genomic co-ordinates of the attention site (e.g., Cas9 cut site) are determined, and features occurring at the attention site are obtained from the feature set file. The output is a file containing the nucleotide sequences along with the genomic co-ordinates of the perfect alignments, genomic co-ordinates of the attention site, and features obtained at that site. Additionally, Ensembl, HGNC, and NCBI gene IDs are attached using `biomaRt`.

A vector of original features can be given as an optional input. If provided, then exorcise additionally computes a one-to-one mapping between original features and reannotated features. For each original feature, exorcise compiles a list of reannotated features across all nucleotide sequences pre-annotated with that original feature. One reannotated feature is selected according to mode; if multimodel, then a feature priority list is used (see example in `tutorial/tutorial/data/symbol_ids_table.csv.gz` and `tutorial/tutorial/data/MRK_List2.rpt.gz` for feature priority lists for RefSeq gene symbols in GRCh38 and GRCm39, respectively). All nucleotide sequences are reannotated with the one selected feature.

An optional specification of control (e.g., non-targeting) nucleotide sequences can be given as an input to exorcise. exorcise will accept these sequences into its algorithm as normal, but these sequences will additionally be analysed to produce a report detailing whether they mapped to any genomic region with features.

## Features
* CRISPRko sequence support, attention site at the Cas9 cut site between -4 and -3 nucleotides protospacer-adjacent motif (PAM)-distal from the PAM.
* Automatic attachment of Ensembl, HGNC, and NCBI gene IDs for human and mouse genomes.
* Support for gzipped file inputs (except for the `.2bit` genome file).
* Retains intermediate files to use as a checkpoint. These files will be accessed and won't be created again in the case that the analysis was stopped and needs to be run again.
* Reannotation of sequences using a pre-generated exorcise library.
* Reannotation of features using a pre-generated exorcise feature mapping file.
* Ad-hoc mode, calling the main exorcise algorithm where the input file contains additional data (e.g., counts data) in addition to the vector of nucleotide sequences.
* Verbose logging.

## Coming soon
* Support for CRISPRa and CRISPRi mode with attention sites based on experimentally determined target sites.
* Support for direct mode where the attention site matches the alignment coordinates.
* Support for arbitrary mode where the attention site is specified by the user.
* Option to disable joining Ensembl, HGNC, and NCBI annotations.
* Checking to see if intermediate files are out of date and need replacing.
* Option to remove intermediate files after a successful run.

## Usage
The main script, `bin/crispieR-lib.R`, performs the reannotation algorithm and is called by the peripheral scripts.

Options include:
* `-i` input file containing a vector of nucleotide sequences and optionally a vector of pre-annotated features. This file should be a plaintext `.tsv`, `.csv`, or `.txt` file, gzipped plaintext file, or `.xls`, or `.xlsx` file.
* `-k` input file specification: can be "tsv", "csv", "txt", "xls", or "xlsx". exorcise will infer the filetype if this is missing. (Optional)
* `-g` integer specifying the column number of the vector of nucleotide sequences in the input file.
* `-n` integer specifying the column number of the vector of pre-annotated features in the input file. (Optional)
* `-a`, either "Human" or "Mouse". This affects where to look to attach Ensembl, HGNC, and NCBI gene IDs.
* `-p` prefix to use for all output files.
* `-c` comma-separated list of control regexes. Any features matching these regexes are detected as control sequences.
* `-d` comma-separated list of control types. Control sequences identified by `-c` will be reannotated to these types.
* `-v` path to the genome in `.2bit` format.
* `-w` path to the reference feature set. Accepts gzipped files. UCSC format. See example in `tutorial/tutorial/data`
* `-y` path to the feature priority set. Accepts gzipped files. See example in `tutorial/tutorial/data`
* `-o` path where all output files will be written.

Call `bin/crispieR-lib.R --help` for detailed usage instructions.

Peripheral scripts `bin/crispieR-cts.R` and `bin/crispieR-up.R` use pre-generated exorcise outputs as sources of reannotation, specifying the exorcise file using `-l`. If `-l` is not specified, then these scripts run in ad-hoc mode and the above options are required.

## Output
The main output of exorcise is the master reannotation file: `6a-LIBRARY_master.tsv`, which contains mappings of the input nucleotide sequences with the alignment coordinates, attention coordinates, and features.

If an original feature vector is supplied, then an inferred mapping file, `6b-LIBRARY_inferred.tsv` is also created. This contains one-to-one mappings of pre-annotated features with exorcise reannotated features.

If any controls were specified, then a negative control report file `6c-LIBRARY_negctrl.tsv` is also created. This file contains only controls and is useful for seeing if any control sequences mapped to non-control regions, i.e. genomic regions with features.

## About
exorcise was developed by Dr Simon Lam, University of Cambridge, mailto:sl681@cam.ac.uk.

## License
This code is licensed under the Creative Commons Zero v1.0 Universal license.








