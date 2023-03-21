#!/usr/bin/env Rscript
library(dplyr)
library(data.table)
library(optparse)

fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-i", "--infile"), type = "character", default = NULL,
                help = "Path to original counts/library file containing authors' gene symbols.", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

infile <- opt$infile
data <- fread(infile)
data <- data %>% mutate(guide = seq)
fwrite(data, infile, sep = "\t")
