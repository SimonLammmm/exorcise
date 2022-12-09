#!/usr/bin/env Rscript
#
# author: Dr Simon Lam
# contact: sl681@cam.ac.uk
#
# Synopsis:
# Update a raw CRISPR counts file with a crispieR library
#
# Description:
#
# This program is to be run after generating a crispieR library. The user supplies
# a crispieR library, original raw counts file, and the column number containing
# sgRNA sequences. The sgRNA column will be used as the join key to the crispieR
# library.
#
# crispieR-cts outputs a corrected crispieR counts file which can be used for
# downstream processing as in a normal CRISPR screening pipeline.
#
# crispieR-cts relies on a crispieR library which has sgRNA and Symbol (updated
# gene annotation) columns. Library files of the expected format are created by
# crispieR-lib, which should be run before this program if you do not already have
# a crispieR library.

#### INIT ####

library(data.table)
library(readxl)
library(dplyr)
library(optparse)
library(logger)

#### VERSION HISTORY ####

# Version   Timestamp             Description
# 0.1       2022-12-07T16-48-10   Initial release

version <- 0.1

##### FUNCTIONS ####

# Main reannotation function
reannotateCts <- function(x) {
  
  log_info("Attaching library ", x$lib, "...")
  lib <- fread(x$lib)
  lib <- lib %>% transmute(sgRNA, gene = Symbol)
  
  log_info("Attaching counts file ", x$cts, "...")
  cts <- fread(x$cts)
  cts_sgRNA_col <- names(cts)[x$sgRNA]
  
  log_info("Using ", cts_sgRNA_col, " to join guide RNA sequences to new annotations...")
  guides <- cts[[cts_sgRNA_col]]
  if(any(grepl("[^ACGTacgt]", guides))) {
    log_error("Non-nucleotide (ACTG) character found in column ", x$sgRNA, ": \"", cts_sgRNA_col, "\" in file: \"", x$cts, "\". Exiting.")
    panic()
  }
  
  log_info("Reannotating...")
  cts <- cts %>% left_join(lib, by = setNames(c("sgRNA"), cts_sgRNA_col))
  cts <- cts %>% unique()
  cts <- cts %>% mutate(guide = paste0("ID_", 1:nrow(cts)))
  
  log_info("Writing new annotations to ", x$out)
  fwrite(cts, x$out, sep = "\t")
  
}

# Helper function to vectorise inputs into a rectangular table
vectoriseArgs <- function(option) {
  option_args <- unlist(opt[option])
  option_first <- option_args[1]
  option_length <- length(option_args)
  if (option_length != n_analyses) {
    if (option_length != 1) {
      log_warn(paste0("--", option, " received ", option_length, " inputs, but expected 1 or ", n_analyses, ". Only the first input, \"", option_first, "\" will be used."))
    }
    return(rep(option_first, n_analyses))
  } else {
    return(option_args)
  }
}

# Helper function to deal with file not found errors
checkIfFileExists <- function(file) {
  if(file.exists(file)) {
    return()
  } else {
    log_error("Input file not found: ", file, ". Quitting.")
    panic()
  }
}

# Stop if error
panic <- function() stop("There were errors while running the script. Check the logfile for details.", call. = F)

#### MAIN ####

if (interactive()) {
  opt <- list()
  opt$cts <- "~/bio/Projects/improving-guide-design/benslimane-paper/author/cts/BenslimaneHarrington2020.GSE150232_Nalm6_sgRNA_read_counts.txt"
  opt$lib <- "~/bio/Projects/DDRcs/lib/EKO/EKO_master.tsv"
  opt$out <- "~/bio/Sandbox/out.tsv"
  opt$sgRNA <- c("1", "2", "3")
}

if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-c", "--cts"), type = "character", default = NULL,
                help = "Path to original counts file(s) containing sgRNA sequences and authors' gene symbols.", metavar = "character"),
    make_option(opt_str = c("-l", "--lib"), type = "character", default = NULL,
                help = "Path to crispieR library file(s) that you would like to use for the reannotation.", metavar = "character"),
    make_option(opt_str = c("-o", "--out"), type = "character", default = NULL,
                help = "Path to save the new counts file(s) after re-annotation", metavar = "character"),
    make_option(opt_str = c("-g", "--sgRNA"), type = "character", default = NULL,
                help = "Column number in the counts file that specifies the sgRNAs", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

logfile <- paste0(dirname(opt$out), "/logfile_crispieR-cts_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
log_appender(appender_tee(logfile))
log_info("Welcome to crispieR-cts, version ", version, ".")

# Find out how many times we need to run the re-annotation
n_analyses <- max(sapply(opt, FUN = length))

# Vectorise inputs, warn if any inputs are ignored
opt2 <- lapply(names(opt), FUN = function(x) vectoriseArgs(x))
names(opt2) <- names(opt)

# Handle file not found errors
files <- unique(opt2$cts, opt2$lib)
for (file in files) {
  checkIfFileExists(file)
}

# Change input to --sgRNA to integer, error if fails
sgRNAcols <- as.integer(opt2$sgRNA)
if(any(is.na(sgRNAcols))) {
  log_error("Non-numeric input found in --sgRNA. Quitting.")
  panic()
} else {
  opt2$sgRNA <- sgRNAcols
}

# Reannotate
opt2 <- as_tibble(opt2)
for (i in 1:n_analyses) {
  commandArgs <- paste(paste0("--", names(opt2[i, ])), opt2[i, ], collapse = " ")
  log_info("Running crispieR-cts, ", i, " of ", n_analyses, "...")
  log_info("cwd: ", getwd())
  log_info("Command: crispieR-cts.R ", commandArgs)
  reannotateCts(opt2[i, ])
}

