#!/usr/bin/env Rscript
#
# author: Dr Simon Lam
# contact: sl681@cam.ac.uk
#
# Synopsis:
# Update gene symbols using an external mapping file. If no mapping file can be
# provided, then one will be generated using crispieR-lib.
#
# Description:
#
# This program expects a counts or library file containing a column with gene
# symbols needing to be updated. If an accompanying mapping file is provided,
# then the symbols are updated according to the mapping file. If not, then a
# new mapping file will be created using the crispieR-lib method, which requires
# a vector of sgRNA sequences, a vector of original symbols, a genome in 2bit
# format, a file describing exon genomic regions on the same co-ordines as the
# 2bit file, and a symbol priorities file.

#### INIT ####

library(data.table)
library(readxl)
library(dplyr)
library(optparse)
library(logger)
library(foreach)

#### VERSION HISTORY ####

# Version   Timestamp             Description
# 0.1       2023-02-20T16:59:44   Initial release
# 0.11      2023-02-27T12:33:25   Added globbed filename functionality

version <- 0.11

##### FUNCTIONS ####

# fread with hanging glob
fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

# Main reannotation function
crispieRup <- function(opt2) {
  
  if(suppressWarnings(is.null(opt2$mapping))) {
    log_info("No mapping supplied. Sending to crispieR-lib using sgRNAs as an ad-hoc library...")
    crispieRlibWrapper(opt2)
    opt2$mapping <- paste0(dirname(opt2$outfile), "/crispieRLib/6b-crispieRLib_inferred.tsv")
  }
  
  log_info("Attaching gene symbol mapping ", opt2$mapping, "...")
  mapping <- fread(opt2$mapping)
  mapping <- mapping %>% transmute(Original_symbol, Symbol)
  
  log_info("Attaching infile ", opt2$infile, "...")
  infile <- fread(opt2$infile)
  infile_symbol_col <- names(infile)[opt2$symbol_col]
  
  log_info("Using ", infile_symbol_col, " to join authors' symbols to reannotated symbols...")
  infile <- infile %>% left_join(mapping, by = setNames(c("Original_symbol"), infile_symbol_col), suffix = c("_cts", ""))
  infile <- infile %>% unique()
  infile <- infile %>% transmute(infile, guide = paste0("ID_", 1:nrow(infile)), Symbol)
  infile <- infile %>% transmute(guide, Symbol, infile)
  
  # If Symbol is empty, then retain original annotation
  infile <- infile %>% mutate(Symbol = case_when(Symbol == "" | is.null(Symbol) | is.na(Symbol) ~ select(infile, all_of(infile_symbol_col))[[1]],
                                                 T ~ Symbol))
  
  log_info("Writing new annotations to ", opt2$outfile)
  fwrite(infile, opt2$outfile, sep = "\t")
  
}

# Call crispieR-lib.R if no mapping is supplied
crispieRlibWrapper <- function(opt2) {
  
  # Obtain sgRNAs
  sg <- foreach(i = 1:length(opt2$infile), .combine = "rbind", .final = function(x) unique(x)) %do% {
    fread(opt2$infile[i]) %>% select(all_of(c(opt2$grna_col[i], opt2$symbol_col[i])))
  }
  
  if(any(grepl("[^ACGTacgt]", unlist(sg[,1])))) {
    log_error("Non-nucleotide (ACTG) character found in input. Exiting.")
    panic()
  }
  
  # Prepare crispieR-lib function call
  crispieRLibopt <- list()
  crispieRLibopt$infiles <- paste0(dirname(opt2$outfile), "/crispieRLib/0-crispieRLib_sgRNAs.tsv")
  crispieRLibopt$kind <- "txt"
  crispieRLibopt$grna_col <- "1"
  crispieRLibopt$symbol_col <- "2"
  crispieRLibopt$project_name <- "crispieRLib"
  crispieRLibopt$outdir <- dirname(opt2$outfile)
  suppressWarnings({
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
  
  # Return to main loop (crispieRup)
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
  opt$infile <- "/Users/simonlam/Library/CloudStorage/OneDrive-UniversityofCambridge/Sandbox/crispieR-update/BenslimaneHarrington2020.GSE150232_Nalm6_sgRNA_read_counts.txt"
  opt$outfile <- "/Users/simonlam/Library/CloudStorage/OneDrive-UniversityofCambridge/Sandbox/crispieR-update/BenslimaneHarrington2020.GSE150232_Nalm6_sgRNA_read_counts_crispieR.txt"
  opt$mapping <- "/Users/simonlam/Library/CloudStorage/OneDrive-UniversityofCambridge/Sandbox/crispieR-update/crispieRLib/6b-crispieRLib_inferred.tsv"
  opt$symbol_col <- "2"
  opt$grna_col <- "1"
  opt$file_genome <- "/Users/simonlam/Library/CloudStorage/OneDrive-UniversityofCambridge/Projects/dev/crispieR/data/hg38.2020-09-22.2bit"
  opt$file_exons <- "/Users/simonlam/Library/CloudStorage/OneDrive-UniversityofCambridge/Projects/dev/crispieR/data/hg38.refseq.exons.tsv"
  opt$file_feature_priorities <- "/Users/simonlam/Library/CloudStorage/OneDrive-UniversityofCambridge/Projects/dev/crispieR/data/symbol_ids_table.csv"
}

if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-i", "--infile"), type = "character", default = NULL,
                help = "Path to original counts/library file containing authors' gene symbols.", metavar = "character"),
    make_option(opt_str = c("-o", "--outfile"), type = "character", default = NULL,
                help = "Path to save the new counts/library file after updating gene symbols.", metavar = "character"),
    make_option(opt_str = c("-n", "--symbol_col"), type = "character", default = NULL,
                help = "Column number in the counts/library file that specifies the authors' gene symbols.", metavar = "character"),
    make_option(opt_str = c("-m", "--mapping"), type = "character", default = NULL,
                help = "Optional. Path to a crispieR inferred mapping file. If this is not specified, then an ad-hoc library will be generated for which --grna_col, --mode, --species, --file_genome, --file_exons, and --file_feature_priorities must be specified.", metavar = "character"),
    make_option(opt_str = c("-g", "--grna_col"), type = "character", default = NULL,
                help = "Optional. Column number in the counts/library file that specifies the sgRNAs. If this is not specified, then --mapping must be specified.", metavar = "character"),
    make_option(opt_str = c("-q", "--mode"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --mapping is not specified. CRISPR screen type: KO (knockout), a (activation), i (inhibition)", metavar = "character"),
    make_option(opt_str = c("-a", "--species"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --mapping is not specified. Was the study done on Human or Mouse?", metavar = "character"),
    make_option(opt_str = c("-z", "--pam"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --mapping is not specified. Protoadjacent motif to append to 3' ends of sgRNA sequences.", metavar = "character"),
    make_option(opt_str = c("-v", "--file_genome"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --mapping is not specified. Path to 2bit assembly file for BLAT. See example.", metavar = "character"),
    make_option(opt_str = c("-w", "--file_exons"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --mapping is not specified. Path to exons file. See example.", metavar = "character"),
    make_option(opt_str = c("-y", "--file_feature_priorities"), type = "character", default = NULL,
                help = "Argument passed to crispieR-lib when --mapping is not specified. Path to feature priorities file. See example.", metavar = "character")

  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

logfile <- paste0(dirname(opt$outfile), "/logfile_crispieR-up_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
dir.create(dirname(logfile), recursive = T, showWarnings = F)
log_appender(appender_tee(logfile))
start_time <- proc.time()
log_info("Welcome to crispieR-up, version ", version, ".")

# Find out how many times we need to run the re-annotation
n_analyses <- max(sapply(opt, FUN = length))

# Vectorise inputs, warn if any inputs are ignored
opt2 <- lapply(names(opt), FUN = function(x) vectoriseArgs(x))
names(opt2) <- names(opt)

# Handle file not found errors
files <- unique(opt2$infile)
for (file in files) {
  checkIfFileExists(Sys.glob(paste0(file, "*"))[1])
}

# Change input to --sgRNA to integer, error if fails
sgRNAcols <- as.integer(opt2$grna_col)
if(any(is.na(sgRNAcols))) {
  log_error("Non-numeric input found in --sgRNA. Quitting.")
  panic()
} else {
  opt2$grna_col <- sgRNAcols
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

opt2$help <- NULL
if(length(opt2$grna_col) == 0) { opt2$grna_col <- NULL }
opt2 <- as_tibble(opt2)

# Update gene symbols
for (i in 1:n_analyses) {
  commandArgs <- paste(paste0("--", names(opt2[i, ])), opt2[i, ], collapse = " ")
  log_info("Running crispieR-up, ", i, " of ", n_analyses, "...")
  log_info("cwd: ", getwd())
  log_info("Command: crispieR-up.R ", commandArgs)
  crispieRup(opt2[i, ])
}

end_time <- proc.time()
log_info("crispieR-up process completed in ", (end_time - start_time)[[3]], " seconds.")