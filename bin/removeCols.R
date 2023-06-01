#!/usr/bin/env Rscript
# Simple column mover for MAGeck
#
# Inputs:
# -i path to file

library(optparse)
library(dplyr)
library(data.table)

removeCols <- function(opt) {
  file <- fread(opt$infile)
  if("gene.orig" %in% names(file)) {
    file <- file %>% select(-gene.orig)
  }
  if("guide.orig" %in% names(file)) {
    file <- file %>% select(-guide.orig)
  }
  if("gene" %in% names(file)) {
    file <- file %>% relocate(gene)
  }
  if("guide" %in% names(file)) {
    file <- file %>% relocate(guide)
  }
 
  return(file)
}

if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-i", "--infile"), type = "character", default = NULL,
                help = "Path to infile.", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

file <- removeCols(opt)
fwrite(file, opt$infile, sep = "\t")