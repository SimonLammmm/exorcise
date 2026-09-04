#!/usr/bin/env Rscript

# Exorcise: exome-guided reannotation of nucleotide sequences.
#
# Author:  Dr Simon Lam, University of Cambridge <sl681@cam.ac.uk>
# Source:  https://github.com/SimonLammmm/exorcise
# Licence: Creative Commons Zero v1.0 Universal
#
# The release history lives in CHANGELOG.md.

EXORCISE_VERSION <- "3.0.1"

# Load order matters. S4Vectors, which arrives with GenomicRanges, exports
# generics that share names with dplyr verbs (rename among them), so the
# Bioconductor stack goes first and the tidyverse masks it rather than the other
# way round.
suppressWarnings(suppressMessages({
  library(optparse)
  library(logger)
  library(GenomicRanges)
  library(rtracklayer)
  library(Biostrings)
  library(data.table)
  library(tidyr)
  library(dplyr)
}))

options(scipen = 999)


# Find the R/ directory belonging to this installation. Earlier versions hard
# coded /exorcise/bin, which meant the Docker image worked and a plain clone did
# not. Set EXORCISE_HOME to override.
exorcise_lib_dir <- function() {
  home <- Sys.getenv("EXORCISE_HOME", unset = "")
  candidates <- character()
  if (nzchar(home)) {
    candidates <- c(candidates, file.path(home, "R"))
  }

  invocation <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(invocation) > 0) {
    script <- normalizePath(sub("^--file=", "", invocation[[1]]), mustWork = FALSE)
    here <- dirname(script)
    candidates <- c(candidates, file.path(here, "..", "R"), file.path(here, "R"), here)
  }
  candidates <- c(candidates, "/exorcise/R", file.path(getwd(), "R"), getwd())

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "common.R"))) {
      return(normalizePath(candidate))
    }
  }
  stop("Could not find the exorcise R/ directory. Set EXORCISE_HOME to the root ",
       "of your exorcise installation.", call. = FALSE)
}

EXORCISE_MODULES <- c(
  "common.R", "options.R", "preprocess.R", "blat.R", "exome.R", "baseedit.R",
  "output.R", "pipeline.R"
)

EXORCISE_LIB <- exorcise_lib_dir()
for (module in EXORCISE_MODULES) {
  source(file.path(EXORCISE_LIB, module))
}


exorcise_banner <- function() {
  paste0("

 ▓█████ ▒██   ██▒ ▒█████   ██▀███   ▄████▄   ██▓  ██████ ▓█████
 ▓█   ▀ ▒▒ █ █ ▒░▒██▒  ██▒▓██ ▒ ██▒▒██▀ ▀█  ▓██▒▒██    ▒ ▓█   ▀
 ▒███   ░░  █   ░▒██░  ██▒▓██ ░▄█ ▒▒▓█    ▄ ▒██▒░ ▓██▄   ▒███
 ▒▓█  ▄  ░ █ █ ▒ ▒██   ██░▒██▀▀█▄  ▒▓▓▄ ▄██▒░██░  ▒   ██▒▒▓█  ▄
 ░▒████▒▒██▒ ▒██▒░ ████▓▒░░██▓ ▒██▒▒ ▓███▀ ░░██░▒██████▒▒░▒████▒
 ░░ ▒░ ░▒▒ ░ ░▓ ░░ ▒░▒░▒░ ░ ▒▓ ░▒▓░░ ░▒ ▒  ░░▓  ▒ ▒▓▒ ▒ ░░░ ▒░ ░
 ░ ░  ░░░░   ░▒ ░  ░ ▒ ▒░   ░▒ ░ ▒░  ░  ▒    ▒ ░░ ░▒  ░ ░ ░ ░  ░
 ░    ░    ░  ░ ░ ░ ▒    ░░   ░ ░         ▒ ░░  ░  ░     ░
 ░  ░ ░    ░      ░ ░     ░     ░ ░       ░        ░     ░  ░
 ░

 exorcise ", EXORCISE_VERSION, "

 Author: Dr Simon Lam, University of Cambridge
 GitHub: https://github.com/SimonLammmm/exorcise

================================================================
\n")
}

exorcise_citation <- function() {
  paste0("
 exorcise ", EXORCISE_VERSION, " was developed by Dr Simon Lam, University of Cambridge.
 https://github.com/SimonLammmm/exorcise

 If exorcise was useful in your work, please cite:

 * Lam S, Thomas JC, Jackson SP, 2024. Genome-aware annotation of CRISPR guides
   validates targets in variant cell lines and enhances discovery in screens.
   Genome Med 16(139). doi:10.1186/s13073-024-01414-4

 * Kent WJ, 2002. BLAT: the BLAST-like alignment tool.
   Genome Res 12(4):656-664. doi:10.1101/gr.229202

=======================================================================\n")
}


# Log to the console and, once an output directory is known, to a file inside it.
start_logging <- function(outdir) {
  if (length(outdir) == 0 || !nzchar(outdir[[1]])) {
    return(invisible(NULL))
  }
  logfile <- file.path(
    outdir[[1]],
    paste0("logfile_exorcise_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
  )
  dir.create(dirname(logfile), recursive = TRUE, showWarnings = FALSE)
  log_appender(appender_tee(logfile))
  invisible(logfile)
}


main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_exorcise_args(args)

  if (isTRUE(opt$version)) {
    cat(EXORCISE_VERSION, "\n", sep = "")
    return(invisible(0L))
  }
  cat(exorcise_banner())
  if (isTRUE(opt$ref)) {
    cat(exorcise_citation())
    return(invisible(0L))
  }

  require_gzip_support()

  start_logging(commasplit(opt$outdir))
  log_info("exorcise ", EXORCISE_VERSION, " starting in ", getwd(), ".")
  log_info("Command: exorcise ", paste(args, collapse = " "))

  opt <- validate_options(opt)
  log_info("Chemistry: ", describe_mode(opt), ".")

  started <- proc.time()
  run_exorcise(opt)
  log_info("Exorcism took ", round((proc.time() - started)[["elapsed"]], 1), " seconds.")
  invisible(0L)
}


if (!interactive()) {
  main()
}
