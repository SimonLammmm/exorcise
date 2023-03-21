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
library(foreach)

#### VERSION HISTORY ####

# Version   Timestamp             Description
# 0.1       2022-12-07T16-48-10   Initial release
# 0.2       2022-12-16T11-48-01   Added crispieR-lib wrapper
# 0.21      2023-02-20T17-00-09   Fixed bugs
# 0.22      2023-02-27T12-32-43   Added globbed filename functionality
# 0.23      2023-03-13T14:36:48   Fixed to retain the "guide" column if already exists in input file
# 0.24      2023-03-13T14:38:33   Fixed to retain the metadata columns other than numerics

version <- 0.24

##### FUNCTIONS ####

# fread with hanging glob
fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

# Main reannotation function
reannotateCts <- function(opt2) {
  
  if(suppressWarnings(is.null(opt2$lib))) {
    log_info("No library supplied. Sending to crispieR-lib using sgRNAs as an ad-hoc library...")
    crispieRlibWrapper(opt2)
    opt2$lib <- paste0(dirname(opt2$out), "/crispieRLib/6a-crispieRLib_master.tsv")
  }
  
  log_info("Attaching library ", opt2$lib, "...")
  lib <- fread(opt2$lib)
  lib <- lib %>% transmute(sgRNA, gene = Symbol)
  
  log_info("Attaching counts file ", opt2$cts, "...")
  cts <- fread(opt2$cts)
  cts_sgRNA_col <- names(cts)[opt2$sgRNA]
  
  log_info("Using ", cts_sgRNA_col, " to join guide RNA sequences to new annotations...")
  guides <- cts[[cts_sgRNA_col]]
  if(any(grepl("[^ACGTacgt]", guides))) {
    log_error("Non-nucleotide (ACTG) character found in column ", opt2$sgRNA, ": \"", cts_sgRNA_col, "\" in file: \"", opt2$cts, "\". Exiting.")
    panic()
  }
  
  log_info("Reannotating...")
  cts <- cts %>% left_join(lib, by = setNames(c("sgRNA"), cts_sgRNA_col), suffix = c(".orig", ""))
  cts <- cts %>% unique()
  
  if("guide" %in% names(cts)) {
    cts$guide.orig <- cts$guide
  }
  
  cts <- cts %>% transmute(cts, guide = paste0("ID_", 1:nrow(cts)), gene)
  #cts <- cts %>% select(guide, gene, all_of(which(sapply(cts, FUN = function(x) class(x) == "numeric" | class(x) == "integer"))))
  
  log_info("Writing new annotations to ", opt2$out)
  fwrite(cts, opt2$out, sep = "\t")
  
}

# Call crispieR-lib.R if no library is supplied
crispieRlibWrapper <- function(opt2) {
  
  # Obtain sgRNAs
  sg <- foreach(i = 1:length(opt2$cts), .combine = "rbind", .final = function(x) unique(x)) %do% {
    fread(opt2$cts[i]) %>% select(all_of(c(opt2$sgRNA[i], opt2$symbol_col[i])))
  }
  
  if(any(grepl("[^ACGTacgt]", unlist(sg[,1])))) {
    log_error("Non-nucleotide (ACTG) character found in input. Exiting.")
    panic()
  }
  
  # Prepare crispieR-lib function call
  crispieRLibopt <- list()
  crispieRLibopt$infiles <- paste0(dirname(opt2$out), "/crispieRLib/0-crispieRLib_sgRNAs.tsv")
  crispieRLibopt$kind <- "txt"
  crispieRLibopt$grna_col <- "1"
  crispieRLibopt$project_name <- "crispieRLib"
  crispieRLibopt$outdir <- dirname(opt2$out)
  suppressWarnings({
    if(!is.null(opt2$symbol_col             )) crispieRLibopt$symbol_col              <- "2"                          else crispieRLibopt$symbol_col              <- NULL
    if(!is.null(opt2$mode                   )) crispieRLibopt$mode                    <- opt2$mode                    else crispieRLibopt$mode                    <- "KO"
    if(!is.null(opt2$PAM                    )) crispieRLibopt$PAM                     <- opt2$PAM                     else crispieRLibopt$PAM                     <- "NGG"
    if(!is.null(opt2$species                )) crispieRLibopt$species                 <- opt2$species                 else crispieRLibopt$species                 <- "Human"
    if(!is.null(opt2$file_genome            )) crispieRLibopt$file_genome             <- opt2$file_genome             else crispieRLibopt$file_genome             <- "data/hg38.2020-09-22.2bit"
    if(!is.null(opt2$file_exons             )) crispieRLibopt$file_exons              <- opt2$file_exons              else crispieRLibopt$file_exons              <- "data/hg38.refseq.exons.tsv"
    if(!is.null(opt2$file_feature_priorities)) crispieRLibopt$file_feature_priorities <- opt2$file_feature_priorities else crispieRLibopt$file_feature_priorities <- "data/symbol_ids_table.csv"
  })
  
  dir.create(dirname(crispieRLibopt$infiles), recursive = T, showWarnings = F)
  fwrite(sg, crispieRLibopt$infiles, sep = "\t")
  
  # Call crispieR-lib
  command <- foreach(i = 1:length(crispieRLibopt)) %do% {
    paste0(" --", names(crispieRLibopt)[i], " ", crispieRLibopt[i])
  } %>% unlist()
  
  command <- paste("crispieR-lib.R", paste0(command, collapse = ""))
  system(command)
  
  # Return to main loop (reannotateCts)
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
  opt$cts <- "cts/BenslimaneHarrington2020_NALM6/author/cts/BenslimaneHarrington2020.GSE150232_Nalm6_sgRNA_read_counts.txt"
  opt$lib <- NULL
  opt$out <- "~/bio/Sandbox/out.tsv"
  opt$sgRNA <- "1"
}

if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-c", "--cts"), type = "character", default = NULL,
                help = "Path to original counts file(s) containing sgRNA sequences and authors' gene symbols.", metavar = "character"),
    make_option(opt_str = c("-l", "--lib"), type = "character", default = NULL,
                help = "Optional. Path to crispieR library file(s) that you would like to use for the reannotation.", metavar = "character"),
    make_option(opt_str = c("-o", "--out"), type = "character", default = NULL,
                help = "Path to save the new counts file(s) after re-annotation", metavar = "character"),
    make_option(opt_str = c("-g", "--sgRNA"), type = "character", default = NULL,
                help = "Column number in the counts file that specifies the sgRNAs", metavar = "character"),
    make_option(opt_str = c("-n", "--symbol_col"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --lib is not specified. Column number(s) containing authors' gene symbols", metavar = "character"),
    make_option(opt_str = c("-q", "--mode"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --lib is not specified. CRISPR screen type: KO (knockout), a (activation), i (inhibition)", metavar = "character"),
    make_option(opt_str = c("-a", "--species"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --lib is not specified. Was the study done on Human or Mouse?", metavar = "character"),
    make_option(opt_str = c("-z", "--pam"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --lib is not specified. Protoadjacent motif to append to 3' ends of sgRNA sequences.", metavar = "character"),
    make_option(opt_str = c("-v", "--file_genome"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --lib is not specified. Path to 2bit assembly file for BLAT. See example.", metavar = "character"),
    make_option(opt_str = c("-w", "--file_exons"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --lib is not specified. Path to exons file. See example.", metavar = "character"),
    make_option(opt_str = c("-y", "--file_feature_priorities"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --lib is not specified. Path to feature priorities file. See example.", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

logfile <- paste0(dirname(opt$out), "/logfile_crispieR-cts_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
dir.create(dirname(logfile), recursive = T, showWarnings = F)
log_appender(appender_tee(logfile))
start_time <- proc.time()
log_info("Welcome to crispieR-cts, version ", version, ".")

# Find out how many times we need to run the re-annotation
n_analyses <- max(sapply(opt, FUN = length))

# Vectorise inputs, warn if any inputs are ignored
opt2 <- lapply(names(opt), FUN = function(x) vectoriseArgs(x))
names(opt2) <- names(opt)

# Handle file not found errors
files <- unique(opt2$cts, opt2$lib)
for (file in files) {
  checkIfFileExists(Sys.glob(paste0(file, "*"))[1])
}

# Change input to --sgRNA to integer, error if fails
sgRNAcols <- as.integer(opt2$sgRNA)
if(any(is.na(sgRNAcols))) {
  log_error("Non-numeric input found in --sgRNA. Quitting.")
  panic()
} else {
  opt2$sgRNA <- sgRNAcols
}

if(!is.null(opt2$symbol_col)) {
  symbolcols <- as.integer(opt2$symbol_col)
  if(any(is.na(symbolcols))) {
    log_error("Non-numeric input found in --symbol_col. Quitting.")
    panic()
  } else {
    opt2$symbol_col <- symbolcols
  }
}

opt2 <- as_tibble(opt2)

# Reannotate
for (i in 1:n_analyses) {
  commandArgs <- paste(paste0("--", names(opt2[i, ])), opt2[i, ], collapse = " ")
  log_info("Running crispieR-cts, ", i, " of ", n_analyses, "...")
  log_info("cwd: ", getwd())
  log_info("Command: crispieR-cts.R ", commandArgs)
  reannotateCts(opt2[i, ])
}

end_time <- proc.time()
log_info("crispieR-cts process completed in ", (end_time - start_time)[[3]], " seconds.")