#!/usr/bin/env Rscript
# Generate LFCs, compute FDRs

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(optparse)


if(!interactive()) {
  option_list <- list(
    make_option(opt_str = c("-c", "--ctsfile"), type = "character", default = NULL,
                help = "Path to counts file. Must contain guide and gene columns. Rest of columns must be numeric counts data, one column per replicate.", metavar = "character"),
    make_option(opt_str = c("-d", "--detfile"), type = "character", default = NULL,
                help = "Path to details Excel sheet in crispr_tools/DDRcs format.", metavar = "character"),
    make_option(opt_str = c("-o", "--outdir"), type = "character", default = NULL,
                help = "Path to folder to output the LFC results. One file will be created for each comparison defined in the detfile.", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
}

infile_cts <- opt$ctsfile
infile_det <- opt$detfile
outdir <- sub("/*$", "/", opt$outdir)
dir.create(outdir, showWarnings = F)

cts <- fread(infile_cts)
det_samp <- read_excel(infile_det, sheet = "Sample details")
det_ctrl <- read_excel(infile_det, sheet = "Control groups")

# Long data format
cts <- cts %>% select(guide, gene, all_of(names(cts)[which(sapply(cts, function(x) class(x)) %in% c("integer", "float", "double", "numeric"))]))
cts <- cts %>% pivot_longer(-c(guide, gene), names_to = "Replicate", values_to = "count")
cts <- cts %>% left_join(det_samp %>% select(Replicate, Sample), by = "Replicate")

# Size factor normalisation
gmeans <- cts %>% group_by(guide) %>% summarise(gmean = exp(mean(log(count))))
cts <- cts %>% left_join(gmeans, by = "guide")
cts <- cts %>% mutate(sfn_count = count / gmean)

# Null distribution
null <- cts %>% filter(gene == "Non-targeting")
# Let central limit theorem apply: assuming normal distribution
null_valid <- null %>% filter(sfn_count != Inf & !is.na(sfn_count))
null_mu <- mean(null_valid$sfn_count, na.rm = T)
null_rho <- sd(null_valid$sfn_count, na.rm = T)

# Fold change, zscore, pval, padj
sfn_counts <- cts %>% group_by(gene, Sample) %>% summarise(sfn_count = mean(sfn_count))
lfcs <- det_ctrl %>% select(`Control sample`, `Test sample`) %>%
  left_join(sfn_counts, by = c("Control sample" = "Sample")) %>%
  left_join(sfn_counts, by = c("Test sample" = "Sample", "gene"), suffix = c(".control", ".test")) %>%
  mutate(FC = sfn_count.test / sfn_count.control,
         LFC = log(FC),
         zscore = (sfn_count.test - sfn_count.control) / null_rho,
         pval = pnorm(zscore, null_mu, null_rho))
lfcs$padj <- p.adjust(lfcs$pval)

# Loop save
dir.create(outdir, showWarnings = F, recursive = T)
for (i in 1:nrow(det_ctrl)) {
  this_ctrl <- det_ctrl$`Control sample`[i]
  this_test <- det_ctrl$`Test sample`[i]
  this_lfcs <- lfcs %>% filter(`Control sample` == this_ctrl & `Test sample` == this_test) %>% select(gene, LFC)
  outfile <- paste0(outdir, "lfc.", this_ctrl, "-", this_test, ".tsv")
  fwrite(this_lfcs, outfile, sep = "\t")
}

