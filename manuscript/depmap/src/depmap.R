#!/usr/bin/env Rscript
# Exorcism of Depmap

# Exorcise data/AvanaRawReadcounts.tsv and data/KYRawReadcounts.csv against GRCh38 and the transcriptome-determined cell line-specific exome
# For each sample, determine the cell line by joining data/ScreenSequenceMap.csv and then Model.csv
# For each cell line, reconstruct the exome from the transcriptome by joining OmicsExpressionTranscriptsTPMLogp1Profile.csv (accept 1 TPM and higher) and exon boundaries in ensembl_109_transcriptome_exon_boundaries.txt (convert from 1-based half-open to 0-based half-open)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(tidyr)
  library(writexl)
})

fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

file.not.exist.or.zero <- function(x) !file.exists(Sys.glob(paste0(x, "*"))[1]) | file.exists(Sys.glob(paste0(x, "*"))[1]) & file.size(Sys.glob(paste0(x, "*"))[1]) == 0

# Load counts
# Download from <https://depmap.org/portal/data_page/?tab=allData>
cts_avana <- fread("data/AvanaRawReadcounts.csv")
cts_ky <- fread("data/KYRawReadcounts.csv")

# Load metadata
cts_coldata <- fread("data/ScreenSequenceMap.csv")
meta_modeldata <- fread("data/OmicsProfiles.csv")

# Load transcriptomes
tpm <- fread("data/OmicsExpressionTranscriptsTPMLogp1Profile.csv")

# Load exon boundaries

boundaries <- fread("data/ensembl_109_transcriptome_exon_boundaries.txt")


#### Obtain cell lines ####

# Obtain run names
runs <- c(names(cts_avana)[-1], names(cts_ky)[-1])

# All run names are found in the metadata - good
runsToRemove <- which(!(runs %in% cts_coldata$SequenceID))
runsToRemove

# Obtain samples of interest from runs of interest
cts_coldata <- cts_coldata %>% filter(SequenceID %in% runs)
models <- cts_coldata$ModelID

# These samples were not found in the metadata - remove
modelsToRemove <- which(!(models %in% meta_modeldata$ModelID))
modelsToRemove
models <- models[-modelsToRemove]

# Obtain cell lines of interest from samples of interest
meta_modeldata <- meta_modeldata %>% filter(ModelID %in% models & Datatype == "rna")
cellLines <- meta_modeldata$ProfileID


#### Determine transcriptomes ####
# Fix column names in the transcriptomes
names(tpm) <- sub("^.+(ENST\\d+).+$", "\\1", names(tpm))
names(tpm)[1] <- "ProfileID"

# Obtain transcriptomes of interest from cell lines of interest
tpm <- tpm %>% filter(ProfileID %in% cellLines)

# Any expression below this (TPM) is considered not expressed
tpmCutOff <- 1

# Binarise TPMs to expressed or not based on the cutoff
tpm <- bind_cols(tpm %>% select(ProfileID),
                 tpm %>% select(-ProfileID) %>% sapply(function(x) x > tpmCutOff))

#### Determine exomes #### 
boundaries <- boundaries %>% transmute(`#chrom` = `Chromosome/scaffold name`,
                                       strand = case_when(boundaries$Strand == 1 ~ "+",                              # Fix sense
                                                          boundaries$Strand == -1 ~ "-"),
                                       exonStarts = `Exon region start (bp)` - 1,                                    # Make boundary coordinatess 0-based
                                       exonEnds = `Exon region end (bp)`,
                                       name = case_when(boundaries$`Gene name` == "" ~ boundaries$`Gene stable ID`,  # Accept symbol otherwise Ensembl gene ID
                                                        boundaries$`Gene name` != "" ~ boundaries$`Gene name`),
                                       `Transcript stable ID`,
                                       `Locus_group` = case_when(grepl("protein.*coding.*", `Gene type`) ~ "protein-coding gene",
                                                                 grepl(".*pseudogene.*", `Gene type`) ~ "pseudogene",
                                                                 grepl(".*RNA.*", `Gene type`) ~ "non-coding RNA",
                                                                 T ~ "other"))

dir.create("data/priorities/", showWarnings = F)
priorities <- boundaries %>% transmute(`Approved_symbol` = name, Locus_group) %>% unique()
prioritiesFile <- "data/priorities/priorities.tsv"
fwrite(priorities, prioritiesFile)

# Loop through each cell line and construct the exome
dir.create("data/exome/", showWarnings = F)

for (i in 1:length(tpm$ProfileID)) {
  thisExomeName <- tpm$ProfileID[i]
  thisExomeFile <- paste0("data/exome/", thisExomeName, ".tsv")
  if(file.not.exist.or.zero(thisExomeFile)) {
    thisExome <- tpm %>%
      filter(ProfileID == thisExomeName) %>%
      pivot_longer(-ProfileID, names_to = "Transcript stable ID", values_to = "Expressed") %>%
      filter(Expressed) %>%
      left_join(boundaries, by = "Transcript stable ID") %>%
      mutate(`#chrom` = paste0("chr", `#chrom`)) %>%
      select(-ProfileID, -Locus_group, -`Transcript stable ID`, -Expressed)
    fwrite(thisExome, thisExomeFile, sep = "\t")
  }
}

#### Run batch exorcise ####
dir.create("data/blat/", showWarnings = F)
blatGenome <- "data/hg38.2020-09-22.2bit"  # Download 2bit genome files from <https://hgdownload.soe.ucsc.edu/downloads.html>
blatFa <- "data/blat/seq.fa"
blatPsl <- "data/blat/blat.psl"

# Run blat once, then use these results again for each exorcism
if(file.not.exist.or.zero(blatPsl)) {
  seq <- unique(c(cts_avana$sgRNA, cts_ky$sgRNA))
  seq <- tibble(V1 = paste0(">exorcise_", seq), V2 = paste0(seq, "NGG"))
  fwrite(seq, blatFa, sep = "\n", col.names = F)
  
system(paste("blat", blatGenome, blatFa, blatPsl, "-stepSize=4 -tileSize=10 -fine -repMatch=2000000 -minScore=20 -minIdentity=100"))
}

# Join cts
cts <- full_join(cts_avana, cts_ky, by = "sgRNA")

# Determine experiments of interest that we will recapitulate in exorcise
experiments <- unique(cts_coldata$ScreenID)
experiments <- experiments[which(experiments != "pDNA")]

# RefSeq reference files
# Download exome annotations for the same assembly from <https://genome.ucsc.edu/cgi-bin/hgTables>, choose NCBI RefSeq All, filter chrom doesn't match *_*
# Download symbol priorities from <https://www.ncbi.nlm.nih.gov/datasets/gene/>, select Symbol and Gene Type
RefSeqExonsFile <- "data/hg38.refseq.exons.tsv"
RefSeqPrioritiesFile <- "data/hsa.priorities.tsv"

# Batch perform exorcise
for (i in 1:length(experiments)) {
  thisExorciseExperiment <- experiments[i]
  dir.create(paste0("runs/", thisExorciseExperiment, "/cts/RefSeq"), recursive = T, showWarnings = F)
  dir.create(paste0("runs/", thisExorciseExperiment, "/cts/transcriptome"), recursive = T, showWarnings = F)
  thisExorciseCtsFile <- paste0("runs/", thisExorciseExperiment, "/cts/cts.tsv")
  
  # Determine counts for this experiment
  if(file.not.exist.or.zero(thisExorciseCtsFile)) {
    thisExorciseColdata <- cts_coldata %>% filter(ScreenID == thisExorciseExperiment)
    thisExorciseRuns <- c(paste0("pDNA_batch_", unique(thisExorciseColdata$pDNABatch)), thisExorciseColdata$SequenceID)
    thisExorciseCts <- cts %>% select(sgRNA, all_of(thisExorciseRuns))
    thisExorciseCts <- thisExorciseCts %>% filter(!is.na(thisExorciseCts[,2]))
    fwrite(thisExorciseCts, thisExorciseCtsFile, sep = "\t")
  }
  
  # Determine exome for this experiment
  thisExorciseModel <- cts_coldata %>% filter(ScreenID == thisExorciseExperiment) %>% select(ModelID) %>% unique() %>% unlist() %>% head(1)
  thisExorciseExome <- meta_modeldata %>% filter(ModelID == thisExorciseModel) %>% select(ProfileID) %>% unique() %>% unlist() %>% head(1)
  thisExorciseExomeFile <- paste0("data/exome/", thisExorciseExome, ".tsv")
  
  # Run exorcise on counts: RefSeq
  thisExorciseOutCtsRefSeq <- paste0("runs/", thisExorciseExperiment, "/cts/RefSeq/exorcise.tsv")
  if(file.not.exist.or.zero(thisExorciseOutCtsRefSeq)) {
    system(paste0("cp ", blatPsl, " runs/", thisExorciseExperiment, "/cts/RefSeq/exorcise.2-hg38.2020-09-22.2bit_BLAT.psl"))
    system(paste0("exorcise -i ", thisExorciseCtsFile, " -o runs/", thisExorciseExperiment, "/cts/RefSeq/ -g 1 -z NGG -v ", blatGenome, " -w ", RefSeqExonsFile))
    thisExorciseCtsRefSeq <- fread(thisExorciseOutCtsRefSeq)
    thisExorciseCtsRefSeq <- thisExorciseCtsRefSeq %>% select(exo_id, exo_symbol, all_of(grep("^exo_|^seqnames$|^start$|^end$|^width$|^strand$|^assembly$|^exome$", names(thisExorciseCtsRefSeq), value = T, invert = T)))
    fwrite(thisExorciseCtsRefSeq, thisExorciseOutCtsRefSeq, sep = "\t")
  }
  
  # Run exorcise on counts: exomes
  thisExorciseOutCtsExome <- paste0("runs/", thisExorciseExperiment, "/cts/transcriptome", thisExorciseExome, ".tsv")
  if(file.not.exist.or.zero(thisExorciseOutCtsExome)) {
    system(paste0("cp ", blatPsl, " runs/", thisExorciseExperiment, "/cts/exorcise.2-hg38.2020-09-22.2bit_BLAT.psl"))
    system(paste0("exorcise -i ", thisExorciseCtsFile, " -o runs/", thisExorciseExperiment, "/cts/transcriptome/ -g 1 -z NGG -v ", blatGenome, " -w ", thisExorciseExomeFile))
    thisExorciseCtsExome <- fread(thisExorciseOutCtsExome)
    thisExorciseCtsExome <- thisExorciseCtsExome %>% select(exo_id, exo_symbol, all_of(grep("^exo_|^seqnames$|^start$|^end$|^width$|^strand$|^assembly$|^exome$", names(thisExorciseCtsExome), value = T, invert = T)))
    fwrite(thisExorciseCtsExome, thisExorciseOutCtsExome, sep = "\t")
  }
}

#### Run batch crispr_pipeline ####
for (i in 1:length(experiments)) {
  thisPipelineExperiment <- experiments[i]
  dir.create(paste0("runs/", thisPipelineExperiment, "/det"), recursive = T, showWarnings = F)
  
  thisPipelineCtsFile <- paste0("runs/", thisPipelineExperiment, "/cts/cts.tsv")
  thisPipelineCts <- fread(thisPipelineCtsFile, nrow = 0) %>% names()
  thisPipelineSeqPDNA <- grep("pDNA", thisPipelineCts, value = T)
  thisPipelineSeqCts <- grep("pDNA", thisPipelineCts, invert = T, value = T)
  thisPipelineSeqCts <- thisPipelineSeqCts[which(thisPipelineSeqCts != "sgRNA")]
  
  # Make details sheet: RefSeq
  thisPipelineCtsRefSeqFile <- paste0("runs/", thisPipelineExperiment, "/cts/RefSeq/exorcise.tsv")
  thisPipelineDetRefSeqFile <- paste0("runs/", thisPipelineExperiment, "/det/det_RefSeq.xlsx")
  
  if(file.not.exist.or.zero(thisPipelineDetRefSeqFile)) {
    thisPipelineDetSheet1 <- tibble(V1 = c("Experiment name", "analysis_version", "file_prefix"),
                                    V2 = c(paste0("runs/", thisPipelineExperiment, "/"), "res/RefSeq", paste0("runs/", thisPipelineExperiment, "/")))
    thisPipelineDetSheet2 <- tibble(Replicate = c(thisPipelineSeqPDNA, thisPipelineSeqCts),
                                    Sample = c(rep("Plasmid", length(thisPipelineSeqPDNA)), rep("Donor", length(thisPipelineSeqCts))))
    thisPipelineDetSheet3 <- tibble(Group = "endpoints",
                                    `Control sample` = "Plasmid",
                                    `Test sample` = "Donor")
    thisPipelineDetSheet4 <- tibble(`Control group` = "endpoints",
                                    Paired = "FALSE",
                                    Method = "drugz",
                                    `Counts file` = "exorcise.tsv",
                                    `Add pseudocount` = 5,
                                    Arguments = "")
    write_xlsx(list("Experiment details" = thisPipelineDetSheet1,
                    "Sample details" = thisPipelineDetSheet2,
                    "Control groups" = thisPipelineDetSheet3,
                    "Analyses" = thisPipelineDetSheet4),
               thisPipelineDetRefSeqFile)
  }
  
  # Make details sheet: exome
  thisPipelineModel <- cts_coldata %>% filter(ScreenID == thisPipelineExperiment) %>% select(ModelID) %>% unique() %>% unlist() %>% head(1)
  thisPipelineExome <- meta_modeldata %>% filter(ModelID == thisPipelineModel) %>% select(ProfileID) %>% unique() %>% unlist() %>% head(1)
  thisPipelineCtsExomeFile <- paste0("runs/", thisPipelineExperiment, "/cts/transcriptome/exorcise.tsv")
  thisPipelineDetExomeFile <- paste0("runs/", thisPipelineExperiment, "/det/det_", thisPipelineExome, ".xlsx")
  
  if(file.not.exist.or.zero(thisPipelineDetExomeFile)) {
    thisPipelineDetSheet1 <- tibble(V1 = c("Experiment name", "analysis_version", "file_prefix"),
                                    V2 = c(paste0("runs/", thisPipelineExperiment, "/"), paste0("res/", thisPipelineExome), paste0("runs/", thisPipelineExperiment, "/")))
    thisPipelineDetSheet2 <- tibble(Replicate = c(thisPipelineSeqPDNA, thisPipelineSeqCts),
                                    Sample = c(rep("Plasmid", length(thisPipelineSeqPDNA)), rep("Donor", length(thisPipelineSeqCts))))
    thisPipelineDetSheet3 <- tibble(Group = "endpoints",
                                    `Control sample` = "Plasmid",
                                    `Test sample` = "Donor")
    thisPipelineDetSheet4 <- tibble(`Control group` = "endpoints",
                                    Paired = "FALSE",
                                    Method = "drugz",
                                    `Counts file` = "exorcise.tsv",
                                    `Add pseudocount` = 5,
                                    Arguments = "")
    write_xlsx(list("Experiment details" = thisPipelineDetSheet1,
                    "Sample details" = thisPipelineDetSheet2,
                    "Control groups" = thisPipelineDetSheet3,
                    "Analyses" = thisPipelineDetSheet4),
               thisPipelineDetExomeFile)
  }
  
  # Run crispr_pipeline.py
  if(file.not.exist.or.zero(paste0("runs/", thisPipelineExperiment, "/res/RefSeq/drugz/files/result.Plasmid-Donor.tsv"))) {
    system(paste0("crispr_pipeline.py --counts runs/", thisPipelineExperiment, "/cts/RefSeq/ --output-dir . ", thisPipelineDetRefSeqFile))
  }
  if(file.not.exist.or.zero(paste0("runs/", thisPipelineExperiment, "/res/", thisPipelineExome, "/drugz/files/result.Plasmid-Donor.tsv"))) {
    system(paste0("crispr_pipeline.py --counts runs/", thisPipelineExperiment, "/cts/transcriptome/ --output-dir . ", thisPipelineDetExomeFile))
  }
}
