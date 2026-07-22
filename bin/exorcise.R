#!/usr/bin/env Rscript

# Exorcise
#
# https://github.com/SimonLammmm/exorcise
#
# author: Dr Simon Lam
# contact: sl681@cam.ac.uk
#
# Synopsis:
# Annotate sequences by aligning to the exome.
#


#### VERSION HISTORY ####
# version       datestamp             description
# 0.9           2023-07-20T17-30-00   evaluation of full release
# 0.9.1         2023-07-21T14-37-00   improved checkpointing
# 0.9.2         2023-07-21T17-00-00   various fixes
# 0.9.3         2023-07-21T17-10-00   various fixes
# 0.9.4         2023-07-24T09-40-00   fix duplicated guide ids for same-locus off-targets
# 0.9.5         2023-07-24T12-20-00   fix multiple control types behaviour
# 0.9.6         2023-07-24T14-14-00   fix multiple control types behaviour
# 0.9.7         2023-07-26T10-17-00   enable harmonisation and control reannotation with exorcised libraries
# 0.9.7.1       2023-07-26T10-17-00   enable harmonisation and control reannotation with exorcised libraries
# 0.9.8         2023-07-26T11-15-00   enable explicit non-harmonisation from pre-exorcised library
# 0.9.8.1       2023-07-26T11-15-00   enable explicit non-harmonisation from pre-exorcised library
# 0.9.9         2023-07-26T14-37-00   switch to NCBI Dataset Gene downloadable feature priority lists
# 0.9.9.1       2023-07-26T14-37-00   switch to NCBI Dataset Gene downloadable feature priority lists
# 1.0           2023-07-28T10-50-00   tested full release
# 1.0.1         2023-07-31T11-00-00   improve handling i/o
# 1.0.2         2023-07-31T11-45-00   warn on potentially incorrect checkpointed psl file
# 1.0.2.1       2023-07-31T11-45-00   warn on potentially incorrect checkpointed psl file
# 1.0.2.2       2023-07-31T11-45-00   warn on potentially incorrect checkpointed psl file
# 1.0.2.3       2023-07-31T11-45-00   warn on potentially incorrect checkpointed psl file
# 1.0.2.4       2023-09-24T15-07-11   enable support for exome files with comment headers
# 1.0.2.5       2023-10-06T18:09:20   enable support for exome files with comment headers
# 1.0.2.6       2023-10-06T18:18:13   enable support for exome files with comment headers
# 1.1           2023-10-11T23-48-30   add exo_id_harm to fix edge case where nonunique sequences specified
# 1.2           2023-11-09T14:08:55   enforce stricter exome column naming according to UCSC Table Browser format
# 1.2.1         2023-11-15T18:19:15   fix output columns in post-hoc exorcise
# 1.3           2024-02-06T16:55:32   add CRISPRi/a support
# 1.4           2024-02-06T22:19:34   speed up harmonisation by vectorising
# 1.4.1         2024-02-07T10:38:30   speed up harmonisation by vectorising
# 1.4.2         2024-02-08T11:14:39   relax CRISPRi/a distance restraint, ignore genomic hits when seqnames not in exome
# 1.4.3         2024-02-09T10:38:17   add arbitrary CRISPR chemistry mode with user-input CRISPR effect range
# 1.5           2024-03-01T15:04:12   add base editor mode
# 1.5.1         2024-03-19T17:47:13   enable inheriting arbitrary columns in exome
# 1.5.2         2024-04-29T13:06:59   add custom base editor mode
# 1.5.3         2024-10-21T17:23:06   fix guide ids in BE mode
# 1.6           2025-10-29T19:35:35   support bystander edits in BE mode
# 2.0           2026-07-21T16:19:30   restructure, enable expression mask, remove support for harmonisation and post-hoc mode
# 2.0.1         2026-07-22T11:48:10   allow expression file to be optional; remove support for globbed inputs

ver <- "2.0.1"

#### INIT ####
suppressWarnings(suppressMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(readxl)
  library(tidyr)
  library(foreach)
  library(logger)
  library(R.utils)
  library(GenomicRanges)
  library(plyranges)
  library(rtracklayer)
  library(stringi)
  library(Biostrings)
}))

options(scipen=999)

#### Init ####

# Testing
if (interactive()) {
  opt <- list()
  opt$infile <-            NULL
  opt$outdir <-            NULL
  opt$seq <-               NULL
  opt$pam <-               NULL
  opt$genome <-            NULL
  opt$exome <-             NULL
  #opt$priorities <-       NULL
  opt$mode <-              NULL
  opt$harm <-              NULL
  opt$control <-           NULL
  opt$control_type <-      NULL
  #opt$library <-          NULL
  opt$ref <-               NULL
  opt$expression <-        NULL
  opt$exprcutoff <-        NULL
  source("common.R")
  source("reannotate.R")
  source("preprocess.R")
  source("psl.R")
  source("exome.R")
  source("baseedit.R")
  source("output.R")
  source("validate.R")
} else {
  source("/exorcise/bin/common.R")
  source("/exorcise/bin/reannotate.R")
  source("/exorcise/bin/preprocess.R")
  source("/exorcise/bin/psl.R")
  source("/exorcise/bin/exome.R")
  source("/exorcise/bin/baseedit.R")
  source("/exorcise/bin/output.R")
  source("/exorcise/bin/validate.R")
}

## Execution
logo <- "


 ▓█████ ▒██   ██▒ ▒█████   ██▀███   ▄████▄   ██▓  ██████ ▓█████ 
 ▓█   ▀ ▒▒ █ █ ▒░▒██▒  ██▒▓██ ▒ ██▒▒██▀ ▀█  ▓██▒▒██    ▒ ▓█   ▀ 
 ▒███   ░░  █   ░▒██░  ██▒▓██ ░▄█ ▒▒▓█    ▄ ▒██▒░ ▓██▄   ▒███   
 ▒▓█  ▄  ░ █ █ ▒ ▒██   ██░▒██▀▀█▄  ▒▓▓▄ ▄██▒░██░  ▒   ██▒▒▓█  ▄ 
 ░▒████▒▒██▒ ▒██▒░ ████▓▒░░██▓ ▒██▒▒ ▓███▀ ░░██░▒██████▒▒░▒████▒
 ░░ ▒░ ░▒▒ ░ ░▓ ░░ ▒░▒░▒░ ░ ▒▓ ░▒▓░░ ░▒ ▒  ░░▓  ▒ ▒▓▒ ▒ ░░░ ▒░ ░
 ░ ░  ░░░   ░▒ ░  ░ ▒ ▒░   ░▒ ░ ▒░  ░  ▒    ▒ ░░ ░▒  ░ ░ ░ ░  ░
 ░    ░    ░  ░ ░ ░ ▒    ░░   ░ ░         ▒ ░░  ░  ░     ░   
 ░  ░ ░    ░      ░ ░     ░     ░ ░       ░        ░     ░  ░
 ░                            


 VERSION"

info <- "

 Author: Dr Simon Lam, University of Cambridge
 GitHub: https://github.com/SimonLammmm/exorcise

================================================================
"


cat(logo, ver, info, "\n\n")

if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-i", "--infile"), type = "character", default = NULL,
                help = "File to be exorcised.", metavar = "character"),
    make_option(opt_str = c("-o", "--outdir"), type = "character", default = NULL,
                help = "Path to destination directory containing output files.", metavar = "character"),
    make_option(opt_str = c("-g", "--seq"), type = "character", default = NULL,
                help = "Sequence column number.", metavar = "character"),
    make_option(opt_str = c("-z", "--pam"), type = "character", default = NULL,
                help = "(optional) PAM sequence.", metavar = "character"),
    make_option(opt_str = c("-q", "--mode"), type = "character", default = NULL,
                help = "(optional) CRISPR chemistry: ko (knockout [default]), a (activation), i (inhibition)", metavar = "character"),
    make_option(opt_str = c("-l", "--library"), type = "character", default = NULL,
                help = "No longer used.", metavar = "character"),
    make_option(opt_str = c("-v", "--genome"), type = "character", default = NULL,
                help = "2bit genome.", metavar = "character"),
    make_option(opt_str = c("-w", "--exome"), type = "character", default = NULL,
                help = "Exome.", metavar = "character"),
    # make_option(opt_str = c("-j", "--id"), type = "character", default = NULL,
    #             help = "(optional, ignored if --library specified) ID column number.", metavar = "character"),
    make_option(opt_str = c("-n", "--harm"), type = "character", default = NULL,
                help = "(optional) Existing annotation column number.", metavar = "character"),
    make_option(opt_str = c("-y", "--priorities"), type = "character", default = NULL,
                help = "No longer used.", metavar = "character"),
    make_option(opt_str = c("-x", "--expression"), type = "character", default = NULL,
                help = "(optional) File with expression values.", metavar = "character"),
    make_option(opt_str = c("-k", "--exprcutoff"), type = "character", default = 10,
                help = "(optional) Value indicating low expression (default: 10)", metavar = "character"),
    make_option(opt_str = c("-c", "--control"), type = "character", default = NULL,
                help = "(optional) Pattern indicating a control guide (comma-separated list).", metavar = "character"),
    make_option(opt_str = c("-d", "--control_type"), type = "character", default = NULL,
                help = "(optional) List of control guide types. Must be the same length as --control_strings (comma-separated list).", metavar = "character"),
    # make_option(opt_str = c("-r", "--just_go"), action = "store", default = F,
    #             help = "(optional) Execute without asking for confirmation.", metavar = "character"),
    make_option(opt_str = c("--ref"), action = "store_true", default = NULL,
                help = "Show citation and quit.", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

if(!is.null(opt$outdir)) {
  logfile <- paste0(opt$outdir, "/", opt$project_name, "/logfile_exorcise_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
  dir.create(dirname(logfile), recursive = T, showWarnings = F)
  log_appender(appender_tee(logfile))
}

# Parse arguments
opt <- fixOpts(opt)

start_time <- proc.time()
log_info("Welcome to exorcise, version ", ver, ".")
log_info("cwd: ", getwd())
command <- foreach(o = opt, .final = function(x) setNames(x, names(opt))) %do% { o }
command <- paste0("--", names(opt), " ", command)
command <- paste0(command, collapse = " ")
log_info("Command: exorcise ", command)

# Call main loop
reannotateLib(opt)
end_time <- proc.time()
log_info("exorcism took ", (end_time - start_time)[[3]], " seconds.")
