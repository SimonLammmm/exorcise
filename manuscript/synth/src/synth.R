#!/usr/bin/env Rscript
# Synthetic CRISPR screen set-up, simulation of common mis-annotations, and analysis
# Supports the Exorcise paper by Lam et al., 2024
# Requires Exorcise (https://github.com/SimonLammmm/exorcise)
# Requires crispr_tools (https://github.com/SimonLammmm/crispr_tools)
# Requires rtracklayer (included in this repo, synth/env/rtracklayer.yaml)

library(data.table)
library(dplyr)
library(tidyselect)
library(Biostrings)
library(foreach)
library(writexl)
library(ggplot2)
library(gridExtra)
library(plotly)
library(DescTools)

fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

outdir <- getwd()
outdir <- sub("/*$", "/", outdir)

dir.create(paste0(outdir, "cts"), showWarnings = F)
dir.create(paste0(outdir, "data"), showWarnings = F)
dir.create(paste0(outdir, "runs"), showWarnings = F)

#### Make fake genes and conditions ####
nconditions <- round(runif(1, 12, 12)) # set number of conditions
conditions <- NULL
for (i in 1:nconditions) {
  keep_generating <- T
  while(keep_generating) {
    initial <- sample(LETTERS, size = 1)
    mid <- sample(c("az", "us", "ish", "op", "ik", "aur", "oz", "ur", "art", "equot", "opan", "im", "etr", "uz", "os", "imp", "his", "gop", "ip", "yr", "yn", "ap", "er", "et", "ov", "iv", "ouf", "pis", "is"), size = sample(1:3, size = 1))
    end <- sample(c("ab", "mab", "tab", "cab", "ib", "nib", "lib", "tib", "ide", "ate", "am", "ram", "pam", "lam", "tam", "in", "lin", "bin", "tin", "fin", "ir", "bir", "tir", "pir", "dir"), size = 1)
    condition <- paste0(initial, paste0(mid, collapse = ""), paste0(end, collapse = ""), collapse = "")
    if (!(condition %in% conditions)) {
      keep_generating <- F
      conditions <- c(conditions, condition)
    }
  }
}
conditions <- c(conditions, "NonTreated")

ngenes <- round(runif(1, 4000, 4000)) # set number of genes
genes <- NULL
for (i in 1:ngenes) {
  keep_generating <- T
  while(keep_generating) {
    initial <- sample(LETTERS, size = 1)
    mid <- sample(LETTERS, size = sample(2:4, size = 1))
    end <- sample(0:9, size = sample(1:3, size = 1))
    gene <- paste0(initial, paste0(mid, collapse = ""), paste0(end, collapse = ""), collapse = "")
    if (!(gene %in% genes)) {
      keep_generating <- F
      genes <- c(genes, gene)
    }
  }
}
ngenes <- ngenes + 1
genes <- c(genes, "NonTargeting")

#### Define construct layout ####
groundTruthGuides <- tibble(                        # Assign guides to genes
  locus = foreach(i = 1:ngenes, .combine = "c") %do% rep(genes[i], 10)) # set number of guides per gene
nguides <- length(groundTruthGuides$locus)
groundTruthGuides <- groundTruthGuides %>%
  mutate(inExon = rep(c(1,0,0,0,0,0,0,0,1,1), nguides/10), # set proportion of intron and exon guides
         gene = case_when(inExon == 1 ~ locus,
                          inExon == 0 ~ "NonTargeting"))

groundTruthEffects <- tibble(gene = genes[1:ngenes-1],                                                                    # Assign ground truth effect sizes per drug for each gene but not the NonTargeting gene
                             isEssential = sample(c(0, 1), ngenes-1, replace = T, prob = c(0.95, 0.05)))                  # Randomly designate some genes as essential. Set probability of non-essential and essential assignment here
for (condition in conditions) {
  groundTruthEffects[condition] = case_when(groundTruthEffects$isEssential == 1 ~ 0.1,
                                            T ~ (1+sqrt(5)/2) ^ rnorm(ngenes-1, 1, 1))
                                            #T ~ sample(c(0.1, 1, 10), ngenes-1, replace = T, prob = c(0.1, 0.8, 0.1))) # Randomly sample ground truth chemogenetic interaction. Set shape of distribution here
}
groundTruthEffects["NonTreated"] = case_when(groundTruthEffects$isEssential == 1 ~ 0.1,
                                             T ~ 1)                                             # Assign the NonTreated ground truth effect
groundTruthEffects[nrow(groundTruthEffects)+1,2:length(groundTruthEffects)] <- 1                # Now add in the NonTargeting gene drug effect
groundTruthEffects$isEssential[nrow(groundTruthEffects)] <- 0                                   # and set the NonTargeting gene non-essential
groundTruthEffects$gene[nrow(groundTruthEffects)] <- "NonTargeting"
fwrite(groundTruthEffects, paste0(outdir, "/data/groundTruthEffects.tsv"), sep = "\t")

groundTruthGuideAnnotations <-                      # Designate essential genes
  left_join(groundTruthGuides, groundTruthEffects, by = "gene") %>%
  group_by(gene) %>%
  mutate(guide = paste0(gene, "_", 1:n()),
         seq = "",
         boundaryStart = "",
         boundaryEnd = "") %>%
  ungroup(gene)

#### Generate DNA ####

# Define settings
bases <- c("A", "C", "G", "T")                                 # Nucleotide alphabet
intervals <- length(groundTruthGuideAnnotations$guide) + 1     # Flanking regions should be one greater than the number of guides to generate
intervalSizeMin <- 31                                          # Flanking region minimum size. Should be at least 19 otherwise guides will overlap
intervalSizeMax <- 120                                         # Flanking region maximum size

# Initialise the chromosome
chrSeq <- NULL
count <- 0

# Generate DNA
for (i in 1:intervals) {
  if (i %% 1000 == 0) { cat(paste("Generating DNA...", i, "out of", intervals, "pieces done.\n")) }
  # Upper flanks also generate guide sequences. Populate the ground truth table.
  if (i < intervals) {
    thisIntervalSize <- round(runif(1, intervalSizeMin, intervalSizeMax))      # Randomise flank size
    groundTruthGuideAnnotations$boundaryStart[i] <- count + 1                  # Record start position of this flank
    groundTruthGuideAnnotations$boundaryEnd[i] <- count + thisIntervalSize + 2 # Record end position of this flank
    count <- count + thisIntervalSize + 2                                      # Update current position in chromosome we're writing to
    
    thisSeq <- sample(bases, thisIntervalSize, replace = T)                    # Generate DNA
    thisSeq <- c(thisSeq, "G", "G")                                            # Append the protospacer adjacent motif (PAM)
    chrSeq <- c(chrSeq, thisSeq)                                               # Commit generated DNA to the chromosome
    
    guideSeq <- thisSeq[(thisIntervalSize-18):(thisIntervalSize-1)]            # Extract the 21mer guide sequence
    guideSeq <- paste0(guideSeq, collapse = "")                                # Commit the 21mer to the ground truth table
    groundTruthGuideAnnotations$seq[i] <- guideSeq     
  # Final flank is after the last guide and so does not go into the ground truth table.
  } else {
    thisIntervalSize <- round(runif(1, intervalSizeMin, intervalSizeMax))      # Randomise flank size
    thisSeq <- sample(bases, thisIntervalSize, replace = T)                    # Generate DNA
    thisSeq                                                                    # Commit generated DNA to the chromosome
  }
}

chrSeq <- paste0(chrSeq, collapse = "") %>% DNAStringSet() %>% setNames("chr1")   # Finalise the chromosome sequence
save(chrSeq, file = paste0(outdir, "data/genome.Rdata"))
system(paste0("~/mambaforge/envs/rtracklayer/bin/R -e 'library(rtracklayer); load(\"", outdir, "data/genome.Rdata\"); export.2bit(chrSeq, \"", outdir, "data/genome.2bit\")'")) # Save to disk in 2bit format
fwrite(groundTruthGuideAnnotations, paste0(outdir, "data/groundTruthGuideAnnotations.tsv"), sep = "\t")

#### Generate exome schemes ####
bmin <- 1
breg <- 10 # Number of introns and exons per gene
btot <- 3 * breg
bmax <- btot - 2 * bmin
overlapsSchemaHelper <- NULL
for (ba in bmin:(btot - breg - 1)) {
  for (bb in max(1, (breg - ba + 1)):(btot - ba - 1)) {
    bc <- btot - ba - bb
    bthis <- c(ba, bb, bc)
    overlapsSchemaHelper <- c(overlapsSchemaHelper, list(c(bthis)))
  }
}

exomesSchema <- tibble(
  guide = groundTruthGuideAnnotations$guide,
  groundTruth = groundTruthGuideAnnotations$gene,                                                    # Define a ground-truth (well-annotated) exome
  overlaps = c(rep(genes[1:ngenes-1], unlist(sample(overlapsSchemaHelper, ceiling(ngenes/3), replace = T))[1:ngenes-1]),
               rep("NonTargeting", nguides))[1:nguides],                                            # Define a misannotated exome: boundaries in the wrong positions
  random = sample(x = genes, size = nguides, replace = T),                                           # Define a misannotated exome: random
  missedtargets = case_when(sample(c(T,F), nguides, replace = T) ~ groundTruthGuideAnnotations$gene,  # Define a misannotated exome: with missed targets
                            T ~ groundTruthGuideAnnotations$locus),
  falseneg = case_when(sample(c(T,F), nguides, replace = T) ~ groundTruthGuideAnnotations$gene,      # Define a misannotated exome: with false negative controls
                       T ~ "NonTargeting")
) %>% mutate(overlaps = case_when(groundTruth == "NonTargeting" ~ "NonTargeting",                  # In overlaps schema, retain NonTargeting guides from the ground truth schema
                                  T ~ overlaps))
fwrite(exomesSchema, paste0(outdir, "data/exomeSchema.tsv"), sep = "\t")

# Export individual exomes
exomesFull <- groundTruthGuideAnnotations %>% left_join(exomesSchema, by = "guide")

exomes <- foreach (x = names(exomesSchema)[which(names(exomesSchema) != "guide")], .final = function(y) setNames(y, names(exomesSchema)[which(names(exomesSchema) != "guide")])) %do% {
  exomesFull %>% transmute(chrom = "chr1",
                           exonStarts = boundaryStart,
                           exonEnds = boundaryEnd,
                           strand = "*",
                           name2 = unlist(exomesFull[,which(names(exomesFull) == x)])) %>%
    group_by(chrom, name2) %>%
    summarise(exonStarts = paste0(exonStarts, collapse = ","),
              exonEnds = paste0(exonEnds, collapse = ","),
              strand = paste0(strand, collapse = ",")) %>%
    filter(name2 != "NonTargeting")
}

for (i in 1:length(exomes)) {
  fwrite(exomes[[i]], paste0(outdir, "data/", names(exomes)[i], ".exome.tsv"), sep = "\t")
}

#### Generate priorities ####
priorities <- tibble(
  Symbol = genes,
  `Gene Type` = "protein-coding gene"
)

fwrite(priorities, paste0(outdir, "data/priorities.tsv"), sep = "\t")

#### Generate synthetic counts data ####
baseline <- 200         # Background expression
randomEffect <- 150     # Range of expression variation due to random effects
multiplier <- 1000      # Number of guides expressed per unit ground truth effect size
counts <- groundTruthGuideAnnotations
for (condition in conditions) {
  counts[paste0("counts_", condition)] <- round(baseline + runif(nguides, -randomEffect, randomEffect) + multiplier * counts[condition]) # Guide counts of drug N depend on the ground-truth dependency of drug N with each guide plus a random effect
}
counts <- counts %>% select(seq, all_of(grep("counts_", names(counts))))
fwrite(counts, paste0(outdir, "cts/cts.tsv"), sep = "\t")

#### Run exorcise ####
for (exome in names(exomes)) {
  system(paste0("exorcise.R -i ", outdir, "cts/cts.tsv -g 1 -o ", outdir, "runs/", exome, "/exorcise/ -v ", outdir, "data/genome.2bit -w ", outdir, "data/", exome, ".exome.tsv -y ", outdir, "data/priorities.tsv -z NGG"))
  thisCts <- fread(paste0(outdir, "runs/", exome, "/exorcise/exorcise.tsv"))
  thisCts <- thisCts %>% select(guide = exo_id, gene = exo_symbol, all_of(grep("counts", names(thisCts), value = T)))
  fwrite(thisCts, paste0(outdir, "runs/", exome, "/exorcise/exorcise.tsv"), sep = "\t")
}

#### Create details.xlsx ####
for (exome in names(exomes)) {
  dir.create(paste0(outdir, "runs/", exome, "/det/"), showWarnings = F)
  sheet1 <- tibble(V1 = c("Experiment name"),
                   V2 = c(paste0("res/")))
  sheet2 <- tibble(Replicate = names(counts)[2:length(names(counts))],
                   Sample = names(counts)[2:length(names(counts))])
  sheet3 <- tibble(Group = rep("endpoints", length(names(counts))-2),
                   `Control sample` = rep(names(counts)[length(names(counts))], length(names(counts))-2),
                   `Test sample` = names(counts)[2:(length(names(counts))-1)])
  sheet4 <- tibble(`Control group` = "endpoints",
                   Paired = "FALSE",
                   Method = "drugz,mageck",
                   `Counts file` = paste0("exorcise.tsv"),
                   `Add pseudocount` = 5,
                   Arguments = "")
  write_xlsx(list("Experiment details" = sheet1,
               "Sample details" = sheet2,
               "Control groups" = sheet3,
               "Analyses" = sheet4),
             paste0(outdir, "runs/", exome, "/det/det.xlsx"))
}

#### Run differential analysis ####
for (exome in names(exomes)) {
  system(paste0("cd ", outdir, "; crispr_pipeline.py --counts ", outdir, "runs/", exome, "/exorcise/ --output-dir ", outdir, "runs/", exome, "/ ", outdir, "runs/", exome, "/det/det.xlsx"))
}

#### Run LFC calculator ####
for (exome in names(exomes)) {
  system(paste0("cd ", outdir, "; lfc.R -c ", outdir, "runs/", exome, "/exorcise/exorcise.tsv -d ", outdir, "runs/", exome, "/det/det.xlsx -o ", outdir, "runs/", exome, "/res/lfc/"))
}

################################################# Evaluation ###################################################

schema <- c("groundTruth", "random", "overlaps", "missedtargets", "falseneg")

dir.create(paste0(outdir, "plots/dz"), showWarnings = F, recursive = T)
dir.create(paste0(outdir, "plots/mag"), showWarnings = F, recursive = T)
dir.create(paste0(outdir, "biplots/dz"), showWarnings = F, recursive = T)
dir.create(paste0(outdir, "biplots/mag"), showWarnings = F, recursive = T)
dir.create(paste0(outdir, "super/dz"), showWarnings = F, recursive = T)
dir.create(paste0(outdir, "super/mag"), showWarnings = F, recursive = T)
dir.create(paste0(outdir, "roc-auc/dz"), showWarnings = F, recursive = T)
dir.create(paste0(outdir, "roc-auc/mag"), showWarnings = F, recursive = T)

################################################# DrugZ ###################################################

#### Plot ####
plotter <- function(scheme, condition) {
  dz <- fread(paste0(outdir, "runs/", scheme, "/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv"))
  lfc <- fread(paste0(outdir, "runs/", scheme, "/res/lfc/lfc.counts_NonTreated-counts_", condition, ".tsv"))
  plotdata <- left_join(dz, lfc, by = c("GENE" = "gene")) %>% left_join(groundTruthEffects, by = c("GENE" = "gene"))
  plotdata <- plotdata %>%
    mutate(groundTruthEffectSize = plotdata[[condition]],
           logEffectSize = log10(groundTruthEffectSize))
  
  fwrite(plotdata, paste0(outdir, "plots/dz/", scheme, ".", condition, ".tsv.gz"), sep = "\t")
  
  p <- ggplot(plotdata, aes(gene = GENE,
                            lfc = LFC,
                            x = rank_synth,
                            y = normZ,
                            colour = groundTruthEffectSize * (1-isEssential),
                            alpha = exp(abs(logEffectSize)) * (1-isEssential),
                            size = exp(abs(logEffectSize)) * (1-isEssential))) +
    geom_point() +
    theme_classic() +
    scale_y_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
    scale_colour_gradient2(name = condition, low = "blue", mid = "grey", high = "red", midpoint = 1, trans = "log", breaks = c(10, 1, 0.1)) +
    theme(legend.position = "none") +
    ggtitle(paste0(condition, ", ", scheme)) +
    xlab("") +
    ylab("")
  return (p)
}

plots <- foreach(scheme = schema, .final = function(y) setNames(y, schema)) %:% {
  foreach(condition = conditions[1:nconditions], .final = function(x) setNames(x, conditions[1:nconditions])) %do% {
    plotter(scheme, condition)
  }
}

plots <- unlist(plots, recursive = F)

pdf(paste0(outdir, "plots/dz/all.pdf"), width = 3 * nconditions, height = 12)
do.call("grid.arrange", c(plots, nrow = length(schema)))
dev.off()

for (plot in names(plots)) {
  htmlwidgets::saveWidget(ggplotly(plots[[plot]]), paste0(outdir, "plots/dz/", plot, ".html"))
}

for (plot in names(plots)) {
  pdf(paste0(outdir, "plots/dz/", plot, ".pdf"), width = 3, height = 3)
  print(plots[[plot]])
  dev.off()
}

#### Biplot ####
biplotter <- function(scheme_a, scheme_b, condition) {
  dz_a <- fread(paste0(outdir, "runs/", scheme_a, "/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv")) %>% select(GENE, normZ) %>% setNames(c("GENE", "x"))
  dz_b <- fread(paste0(outdir, "runs/", scheme_b, "/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv")) %>% select(GENE, normZ) %>% setNames(c("GENE", "y"))
  plotdata <- left_join(dz_a, dz_b, by = "GENE") %>% left_join(groundTruthEffects, by = c("GENE" = "gene"))
  plotdata <- plotdata %>%
    mutate(groundTruthEffectSize = plotdata[[condition]],
           logEffectSize = log10(groundTruthEffectSize))
  
  fwrite(plotdata, paste0(outdir, "biplots/dz/", condition, ".", scheme_a, ".", scheme_b, ".tsv.gz"), sep = "\t")
  
  q <- ggplot(plotdata, aes(gene = GENE,
                            x = x,
                            y = y,
                            colour = groundTruthEffectSize * (1-isEssential),
                            alpha = exp(abs(logEffectSize)) * (1-isEssential),
                            size = exp(abs(logEffectSize)) * (1-isEssential))) +
    geom_point() +
    theme_classic() +
    scale_x_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
    scale_y_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
    scale_colour_gradient2(name = condition, low = "blue", mid = "grey", high = "red", midpoint = 1, trans = "log", breaks = c(10, 1, 0.1)) +
    theme(legend.position = "none") +
    ggtitle(paste0(condition, ", ", scheme_a, " vs. ", scheme_b)) +
    xlab(scheme_a) +
    ylab(scheme_b)
}

biplots <- foreach(condition = conditions[1:nconditions], .final = function(x) setNames(x, conditions[1:nconditions])) %do% {
  foreach(scheme_a = schema, .final = function(y) setNames(y, schema)) %do% {
    foreach(scheme_b = schema, .final = function(x) setNames(x, schema)) %do% {
      biplotter(scheme_a, scheme_b, condition)
    }
  }
}

biplots <- biplots %>% unlist(recursive = F) %>% unlist(recursive = F)

for (biplot in names(biplots)) {
  htmlwidgets::saveWidget(ggplotly(biplots[[biplot]]), paste0(outdir, "biplots/dz/", biplot, ".html"))
}

for (biplot in names(biplots)) {
  pdf(paste0(outdir, "biplots/dz/", biplot, ".pdf"), width = 3, height = 3)
  print(biplots[[biplot]])
  dev.off()
}


#### Super- and subnumerary analysis ####
for (condition in conditions[1:length(conditions)-1]) {
  
  # Obtain super- and subnumerary hits data for this condition
  subnumerary <- fread(paste0(outdir, "runs/falseneg/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv")) %>% select(GENE, normZ, numObs) %>% setNames(c("GENE", "y", "n"))
  supernumerary <- fread(paste0(outdir, "runs/missedtargets/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv")) %>% select(GENE, normZ, numObs) %>% setNames(c("GENE", "y", "n"))
  normnumerary <- fread(paste0(outdir, "runs/groundTruth/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv")) %>% select(GENE, normZ) %>% setNames(c("GENE", "x"))
  plotdata <- bind_rows(subnumerary, supernumerary) %>% left_join(normnumerary, by = "GENE") %>% left_join(groundTruthEffects, by = c("GENE" = "gene"))
  plotdata <- plotdata %>%
    mutate(groundTruthEffectSize = plotdata[[condition]],
           logEffectSize = log10(groundTruthEffectSize))
  
  fwrite(plotdata, paste0(outdir, "super/dz/", condition, "_", "biplot", ".tsv.gz"), sep = "\t")
  
  # # Plot rank plots for each number of guides per gene
  for (myN in min(plotdata$n):max(plotdata$n)) {
    plotdata2 <- plotdata %>% filter(n == myN)
    
    fwrite(plotdata2, paste0(outdir, "super/dz/", condition, "_", myN, "gpg.tsv.gz"), sep = "\t")
    
    p <- ggplot(plotdata2, aes(gene = GENE,
                               x = rank(y),
                               y = y,
                               colour = groundTruthEffectSize * (1-isEssential),
                               alpha = exp(abs(logEffectSize)) * (1-isEssential),
                               size = exp(abs(logEffectSize)) * (1-isEssential))) +
      geom_point() +
      theme_classic() +
      scale_y_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
      scale_colour_gradient2(name = condition, low = "blue", mid = "grey", high = "red", midpoint = 1, trans = "log", breaks = c(10, 1, 0.1)) +
      theme(legend.position = "none") +
      ggtitle(paste0(condition, ", ", myN, " guides per gene")) +
      xlab("") +
      ylab("")
    
    pdf(paste0(outdir, "super/dz/", condition, "_", myN, "gpg.pdf"), 3, 3)
    print(p)
    dev.off()
    htmlwidgets::saveWidget(ggplotly(p), paste0(outdir, "super/dz/", condition, "_", myN, "gpg.html"))
  }
  
  # Plot biplots of super- and subnumerary hits versus ground truth
  q <- ggplot(plotdata, aes(gene = GENE,
                            x = x,
                            y = y,
                            colour = factor(n),
                            alpha = exp(abs(logEffectSize)) * (1-isEssential),
                            size = exp(abs(logEffectSize)) * (1-isEssential))) +
    geom_point() +
    theme_classic() +
    scale_alpha(guide = "none") +
    scale_size(guide = "none") +
    scale_colour_discrete(name = "") +
    ggtitle(paste0(condition)) +
    xlab("Ground truth") +
    ylab("Mis-annotated")
  
  pdf(paste0(outdir, "super/dz/", condition, "_", "biplot", ".pdf"), 3, 3)
  print(q)
  dev.off()
  htmlwidgets::saveWidget(ggplotly(q), paste0(outdir, "super/dz", condition, "_", "biplot", ".html"))
  
}


#### Sensitivity plots ####
scheme_b = "groundTruth"

data <- tibble(condition = character(0),
               scheme = character(0),
               auc = double(0),
               precision = double(0),
               recall = double(0))

for (condition in conditions) {
  for (scheme_a in schema) {
    if (scheme_a == "groundTruth" | condition == "NonTreated") {
      next
    }
    thisRocB = fread(paste0(outdir, "runs/", scheme_b, "/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv")) %>%
      filter(!grepl("Non-targeting", GENE)) %>%
      rowwise() %>%
      mutate(fdr = min(fdr_synth, fdr_supp)) %>%
      ungroup() %>%
      mutate(hit = case_when(fdr <= 0.5 ~ T, T ~ F), sign_gt = sign(normZ)) %>%
      select(GENE, hit_gt = hit, sign_gt)
    thisRocA = fread(paste0(outdir, "runs/", scheme_a, "/res/drugz/files/result.counts_NonTreated-counts_", condition, ".tsv")) %>%
      filter(!grepl("Non-targeting", GENE)) %>%
      rowwise() %>%
      mutate(fdr = min(fdr_synth, fdr_supp)) %>%
      ungroup() %>%
      mutate(hit = case_when(fdr <= 0.5 ~ T, T ~ F)) %>%
      left_join(thisRocB, by = "GENE") %>%
      filter(hit) %>%
      mutate(correct = case_when(hit == hit_gt ~ 1 & sign(normZ) == sign_gt, T ~ 0), incorrect = case_when(correct == 1 ~ 0, T ~ 1)) %>%
      arrange(fdr) %>%
      mutate(n = 1:n(), truepos = cumsum(correct), falsepos = cumsum(incorrect), tpr = truepos / sum(correct), fpr = falsepos / sum(incorrect))
    auc <- AUC(x = thisRocA$fpr, y = thisRocA$tpr)
    gt <- sum(thisRocB$hit_gt)
    tp <- sum(thisRocA$correct)
    fp <- sum(thisRocA$incorrect)
    p <- ggplot(thisRocA, aes(x = fpr, y = tpr)) +
      geom_line(linewidth = 1) +
      geom_area(fill = "grey", alpha = 0.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_text(x = 0.98, y = 0.02, label = paste0("AUC: ", round(auc,4), "\nTrue positives: ", tp, "\nFalse positives: ", fp, "\nActual positives: ", gt), hjust = 1, vjust = 0) +
      theme_classic() +
      xlab("False positive rate") +
      ylab("True positive rate") +
      ggtitle(paste0(condition, " vs control: ", scheme_a))
    pdf(paste0(outdir, "roc-auc/dz/", condition, "_", scheme_a, "_rocauc", ".pdf"), 3.2, 3.2)
    print(p)
    dev.off()
    htmlwidgets::saveWidget(ggplotly(p), paste0(outdir, "roc-auc/dz/", condition, "_", scheme_a, "_rocauc", ".html"))
    data <- bind_rows(data,
                      tibble(condition = condition,
                             scheme = scheme_a,
                             auc = auc,
                             precision = tp / (tp + fp),
                             recall = tp / gt))
  }
}
fwrite(data, paste0(outdir, "roc-auc/dz/roc-auc.tsv.gz"), sep = "\t")

################################################# MAGeCK ###################################################

#### Plot ####
plotter <- function(scheme, condition) {
  mag <- fread(paste0(outdir, "runs/", scheme, "/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt"))
  lfc <- fread(paste0(outdir, "runs/", scheme, "/res/lfc/lfc.counts_NonTreated-counts_", condition, ".tsv"))
  plotdata <- left_join(mag, lfc, by = c("id" = "gene")) %>% left_join(groundTruthEffects, by = c("id" = "gene"))
  plotdata <- plotdata %>%
    mutate(groundTruthEffectSize = plotdata[[condition]],
           logEffectSize = log10(groundTruthEffectSize)) %>%
    rowwise() %>%
    mutate(logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                 `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`))) %>%
  ungroup() %>%
    mutate(rank = rank(logAbsMag))
  
  fwrite(plotdata, paste0(outdir, "plots/mag/", scheme, ".", condition, ".tsv.gz"), sep = "\t")
  
  p <- ggplot(plotdata, aes(gene = id,
                            lfc = LFC,
                            x = rank,
                            y = logAbsMag,
                            colour = groundTruthEffectSize * (1-isEssential),
                            alpha = exp(abs(logEffectSize)) * (1-isEssential),
                            size = exp(abs(logEffectSize)) * (1-isEssential))) +
    geom_point() +
    theme_classic() +
    scale_y_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
    scale_colour_gradient2(name = condition, low = "blue", mid = "grey", high = "red", midpoint = 1, trans = "log", breaks = c(10, 1, 0.1)) +
    theme(legend.position = "none") +
    ggtitle(paste0(condition, ", ", scheme)) +
    xlab("") +
    ylab("")
  return (p)
}

plots <- foreach(scheme = schema, .final = function(y) setNames(y, schema)) %do% {
  foreach(condition = conditions[1:nconditions], .final = function(x) setNames(x, conditions[1:nconditions])) %do% {
    plotter(scheme, condition)
  }
}

plots <- unlist(plots, recursive = F)

pdf(paste0(outdir, "plots/mag/all.pdf"), width = 3 * nconditions, height = 12)
do.call("grid.arrange", c(plots, nrow = length(schema)))
dev.off()

for (plot in names(plots)) {
  htmlwidgets::saveWidget(ggplotly(plots[[plot]]), paste0(outdir, "plots/mag/", plot, ".html"))
}

for (plot in names(plots)) {
  pdf(paste0(outdir, "plots/mag/", plot, ".pdf"), width = 3, height = 3)
  print(plots[[plot]])
  dev.off()
}

#### Biplot ####
biplotter <- function(scheme_a, scheme_b, condition) {
  mag_a <- fread(paste0(outdir, "runs/", scheme_a, "/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt")) %>% rowwise() %>% transmute(id, logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                                                                                                                                                                                `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`))) %>% setNames(c("id", "x"))
  mag_b <- fread(paste0(outdir, "runs/", scheme_b, "/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt")) %>% rowwise() %>% transmute(id, logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                                                                                                                                                                               `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`))) %>% setNames(c("id", "y"))
  plotdata <- left_join(mag_a, mag_b, by = "id") %>% left_join(groundTruthEffects, by = c("id" = "gene"))
  plotdata <- plotdata %>%
    ungroup() %>%
    mutate(groundTruthEffectSize = plotdata[[condition]],
           logEffectSize = log10(groundTruthEffectSize))
  
  fwrite(plotdata, paste0(outdir, "biplots/mag/", condition, ".", scheme_a, ".", scheme_b, ".tsv.gz"), sep = "\t")
  
  q <- ggplot(plotdata, aes(gene = id,
                            x = x,
                            y = y,
                            colour = groundTruthEffectSize * (1-isEssential),
                            alpha = exp(abs(logEffectSize)) * (1-isEssential),
                            size = exp(abs(logEffectSize)) * (1-isEssential))) +
    geom_point() +
    theme_classic() +
    scale_x_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
    scale_y_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
    scale_colour_gradient2(name = condition, low = "blue", mid = "grey", high = "red", midpoint = 1, trans = "log", breaks = c(10, 1, 0.1)) +
    theme(legend.position = "none") +
    ggtitle(paste0(condition, ", ", scheme_a, " vs. ", scheme_b)) +
    xlab(scheme_a) +
    ylab(scheme_b)
}

biplots <- foreach(condition = conditions[1:nconditions], .final = function(x) setNames(x, conditions[1:nconditions])) %do% {
  foreach(scheme_a = schema, .final = function(y) setNames(y, schema)) %do% {
    foreach(scheme_b = schema, .final = function(x) setNames(x, schema)) %do% {
      biplotter(scheme_a, scheme_b, condition)
    }
  }
}

biplots <- biplots %>% unlist(recursive = F) %>% unlist(recursive = F)

for (biplot in names(biplots)) {
  htmlwidgets::saveWidget(ggplotly(biplots[[biplot]]), paste0(outdir, "biplots/mag/", biplot, ".html"))
}

for (biplot in names(biplots)) {
  pdf(paste0(outdir, "biplots/mag/", biplot, ".pdf"), width = 3, height = 3)
  print(biplots[[biplot]])
  dev.off()
}


#### Super- and subnumerary analysis ####
for (condition in conditions[1:length(conditions)-1]) {
  
  # Obtain super- and subnumerary hits data for this condition
  subnumerary <- fread(paste0(outdir, "runs/falseneg/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt")) %>% rowwise() %>% transmute(id, logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                                                                                                                                                                                `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`)), num) %>% setNames(c("GENE", "y", "n"))
  supernumerary <- fread(paste0(outdir, "runs/missedtargets/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt")) %>% rowwise() %>% transmute(id, logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                                                                                                                                                                                       `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`)), num) %>% setNames(c("GENE", "y", "n"))
  normnumerary <- fread(paste0(outdir, "runs/groundTruth/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt")) %>% rowwise() %>% transmute(id, logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                                                                                                                                                                                    `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`))) %>% setNames(c("GENE", "x"))
  plotdata <- bind_rows(subnumerary, supernumerary) %>% left_join(normnumerary, by = "GENE") %>% left_join(groundTruthEffects, by = c("GENE" = "gene"))
  plotdata <- plotdata %>%
    ungroup() %>%
    mutate(groundTruthEffectSize = plotdata[[condition]],
           logEffectSize = log10(groundTruthEffectSize))
  
  fwrite(plotdata, paste0(outdir, "super/mag/", condition, "_", "biplot", ".tsv.gz"), sep = "\t")
  
  # # Plot rank plots for each number of guides per gene
  for (myN in min(plotdata$n):max(plotdata$n)) {
    plotdata2 <- plotdata %>% filter(n == myN)
    
    fwrite(plotdata2, paste0(outdir, "super/mag/", condition, "_", myN, "gpg.tsv.gz"), sep = "\t")
    
    p <- ggplot(plotdata2, aes(gene = GENE,
                               x = rank(y),
                               y = y,
                               colour = groundTruthEffectSize * (1-isEssential),
                               alpha = exp(abs(logEffectSize)) * (1-isEssential),
                               size = exp(abs(logEffectSize)) * (1-isEssential))) +
      geom_point() +
      theme_classic() +
      scale_y_continuous(limits = c(-10, 10), oob = scales::oob_squish) +
      scale_colour_gradient2(name = condition, low = "blue", mid = "grey", high = "red", midpoint = 1, trans = "log", breaks = c(10, 1, 0.1)) +
      theme(legend.position = "none") +
      ggtitle(paste0(condition, ", ", myN, " guides per gene")) +
      xlab("") +
      ylab("")
    
    pdf(paste0(outdir, "super/mag/", condition, "_", myN, "gpg.pdf"), 3, 3)
    print(p)
    dev.off()
    htmlwidgets::saveWidget(ggplotly(p), paste0(outdir, "super/mag/", condition, "_", myN, "gpg.html"))
  }
  
  # Plot biplots of super- and subnumerary hits versus ground truth
  q <- ggplot(plotdata, aes(gene = GENE,
                            x = x,
                            y = y,
                            colour = factor(n),
                            alpha = exp(abs(logEffectSize)) * (1-isEssential),
                            size = exp(abs(logEffectSize)) * (1-isEssential))) +
    geom_point() +
    theme_classic() +
    scale_alpha(guide = "none") +
    scale_size(guide = "none") +
    scale_colour_discrete(name = "") +
    ggtitle(paste0(condition)) +
    xlab("Ground truth") +
    ylab("Mis-annotated")
  
  pdf(paste0(outdir, "super/mag/", condition, "_", "biplot", ".pdf"), 3, 3)
  print(q)
  dev.off()
  htmlwidgets::saveWidget(ggplotly(q), paste0(outdir, "super/mag/", condition, "_", "biplot", ".html"))
  
}


#### Sensitivity plots ####
scheme_b = "groundTruth"

data <- tibble(condition = character(0),
               scheme = character(0),
               auc = double(0),
               precision = double(0),
               recall = double(0))

for (condition in conditions) {
  for (scheme_a in schema) {
    if (scheme_a == "groundTruth" | condition == "NonTreated") {
      next
    }
    thisRocB = fread(paste0(outdir, "runs/", scheme_b, "/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt")) %>%
      filter(!grepl("Non-targeting", id)) %>%
      rowwise() %>%
      mutate(fdr = min(abs(`neg|fdr`), abs(`pos|fdr`))) %>%
      mutate(logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                   `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`))) %>%
      ungroup() %>%
      mutate(hit = case_when(fdr <= 0.95 ~ T, T ~ F), sign_gt = sign(logAbsMag)) %>%
      select(id, hit_gt = hit, sign_gt)
    thisRocA = fread(paste0(outdir, "runs/", scheme_a, "/res/mageck/files/result.counts_NonTreated-counts_", condition, ".gene_summary.txt")) %>%
      filter(!grepl("Non-targeting", id)) %>%
      rowwise() %>%
      mutate(fdr = min(abs(`neg|fdr`), abs(`pos|fdr`))) %>%
      mutate(logAbsMag = case_when(`neg|score` == min(`neg|score`, `pos|score`) ~ log10(`neg|score`),
                                   `pos|score` == min(`neg|score`, `pos|score`) ~ -log10(`pos|score`))) %>%
      ungroup() %>%
      mutate(hit = case_when(fdr <= 0.95 ~ T, T ~ F)) %>%
      left_join(thisRocB, by = "id") %>%
      filter(hit) %>%
      mutate(correct = case_when(hit == hit_gt ~ 1 & sign(logAbsMag) == sign_gt, T ~ 0), incorrect = case_when(correct == 1 ~ 0, T ~ 1)) %>%
      arrange(fdr) %>%
      mutate(n = 1:n(), truepos = cumsum(correct), falsepos = cumsum(incorrect), tpr = truepos / sum(correct), fpr = falsepos / sum(incorrect))
    auc <- AUC(x = thisRocA$fpr, y = thisRocA$tpr)
    gt <- sum(thisRocB$hit_gt)
    tp <- sum(thisRocA$correct)
    fp <- sum(thisRocA$incorrect)
    p <- ggplot(thisRocA, aes(x = fpr, y = tpr)) +
      geom_line(linewidth = 1) +
      geom_area(fill = "grey", alpha = 0.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      geom_text(x = 0.98, y = 0.02, label = paste0("AUC: ", round(auc,4), "\nTrue positives: ", tp, "\nFalse positives: ", fp, "\nActual positives: ", gt), hjust = 1, vjust = 0) +
      theme_classic() +
      xlab("False positive rate") +
      ylab("True positive rate") +
      ggtitle(paste0(condition, " vs control: ", scheme_a))
    pdf(paste0(outdir, "roc-auc/mag/", condition, "_", scheme_a, "_rocauc", ".pdf"), 3.2, 3.2)
    print(p)
    dev.off()
    htmlwidgets::saveWidget(ggplotly(p), paste0(outdir, "roc-auc/mag/", condition, "_", scheme_a, "_rocauc", ".html"))
    data <- bind_rows(data,
                      tibble(condition = condition,
                             scheme = scheme_a,
                             auc = auc,
                             precision = tp / (tp + fp),
                             recall = tp / gt))
  }
}
fwrite(data, paste0(outdir, "roc-auc/mag/roc-auc.tsv.gz"), sep = "\t")


