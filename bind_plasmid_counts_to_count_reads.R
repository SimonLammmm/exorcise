#!/usr/bin/env Rscript
library(dplyr)
library(data.table)
library(optparse)

fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-i", "--infile"), type = "character", default = NULL,
                help = "Path to original counts file", metavar = "character"),
    make_option(opt_str = c("-p", "--plasmids"), type = "character", default = NULL,
                help = "Path to plasmids counts file to attach to original counts file", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

infile <- opt$infile
plasmidsfile <- opt$plasmids
data <- fread(infile)
plasmids <- fread(plasmidsfile, header = F)
plasmids_sampleName <- sub("^.+\\.(.+?)\\.rawcount.txt$", "\\1", plasmidsfile)
names(plasmids) <- c("guide", plasmids_sampleName)
data <- data %>% left_join(plasmids, by = "guide")
fwrite(data, infile, sep = "\t")
