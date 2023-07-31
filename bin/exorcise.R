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

ver <- "1.0.2.2"

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
}))

options(scipen=999)


#### FUNCTIONS ####

# Create directory tree and write file
dirwrite <- function(x, dir, ...) {
  if (!file.exists(dir)) {
    dir.create(dirname(dir), recursive = T, showWarnings = F)
  }
  fwrite(x, dir, ...)
  invisible()
}

# Split strings
commasplit <- function(x) unlist(strsplit(x, ","))

# Return true if file doesn't exist or has zero size
file.not.exist.or.zero <- function(x) !file.exists(Sys.glob(paste0(x, "*"))[1]) | file.exists(Sys.glob(paste0(x, "*"))[1]) & file.size(Sys.glob(paste0(x, "*"))[1]) == 0

# fread with hanging glob
# fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)


#### MAIN FUNCTION ####

reannotateLib <- function(opt) {
  
  if(opt$adhoc) {
    log_info("Using ad-hoc mode: genome and exome specified.")
    
    authors <- importAuthorsLib(opt)
    blats <- list(file_genome = opt$genome,
                  file_exons = opt$exome,
                  file_feature_priorities = opt$priorities,
                  file_sgRNAs = paste0(opt$outdir, "/exorcise.1-seq.fa"),
                  file_psl = paste0(opt$outdir, "/exorcise.2-",  sub(".+/(.+?)$", "\\1", opt$genome), "_BLAT.psl"),
                  file_genomic_ranges = paste0(opt$outdir, "/exorcise.3-", sub(".+/(.+?)$", "\\1", opt$genome), "_genomicRanges.tsv"),
                  file_genomic_seqSpecs = paste0(opt$outdir, "/exorcise.3-", sub(".+/(.+?)$", "\\1", opt$genome), "_genomicSeqSpecs.tsv"),
                  file_genomic_seqs = paste0(opt$outdir, "/exorcise.3-", sub(".+/(.+?)$", "\\1", opt$genome), "_genomicSeqs.fa"),
                  file_genomic_ranges_matched = paste0(opt$outdir, "/exorcise.4-", sub(".+/(.+?)$", "\\1", opt$exome), "_exonHits.tsv"),
                  file_genomic_ranges_distances = paste0(opt$outdir, "/exorcise.5-", sub(".+/(.+?)$", "\\1", opt$exome), "_exonDist.tsv"),
                  file_exorcise_master_out = paste0(opt$outdir, "/exorcise.tsv"))
    
    extractGuides(opt, authors, blats)
    runBlat(opt, authors, blats)
    runPtgr(blats)
    runGrtem(blats)
    
    # generate mapping
    genome_hits <- fread(blats$file_genomic_ranges, colClasses = "character") %>% transmute(exo_seq, exo_target, exo_cut)
    exome_hits <- fread(blats$file_genomic_ranges_matched, colClasses = "character") %>% transmute(exo_seq, exo_cut, exo_symbol) %>% unique()
    all_mappings <- left_join(genome_hits, exome_hits, by = c("exo_seq", "exo_cut"), relationship = "many-to-many") %>% unique()
    
    premaster <- left_join(authors, all_mappings, by = "exo_seq", relationship = "many-to-many")
    master <- exorcisemaster(premaster, opt)
    
    outfile_m_tsv <- blats$file_exorcise_master_out
    dirwrite(master, outfile_m_tsv, sep = "\t")
    log_info("Wrote master library to ", outfile_m_tsv)
    
  } else {
    log_info("Using post-hoc mode: exorcised library specified.")
    reannotateExisting(opt)
  }
}

### OTHER FUNCTIONS ####

# Import authors' original library files, extracting original symbols, detecting/creating a primary key, and removing adapters where necessary
importAuthorsLib <- function(opt) {
  log_info("Opening file ", opt$infile, "...")
  authors <- fread(opt$infile)
  authors <- authors %>% mutate(exo_id = paste0("exorcise_", 1:n()))
  authors <- authors %>%                               # read authors' file
    as_tibble() %>%
    relocate(exo_id, exo_seq = opt$seq, exo_orig = opt$harm) %>%     # select appropriate columns
    unique()
  return(authors)
}

# Use authors' original library to extract guide RNA sequences for alignment
extractGuides <- function(opt, authors, blats) {
  log_info("Extracting guide sequences...")
  sgRNAs <- paste0(">", authors$exo_id, "\n", authors$exo_seq, opt$pam)
  dirwrite(list(sgRNAs), blats$file_sgRNAs, sep = "\n", quote = F)
  return()
}

# Run BLAT
runBlat <- function(opt, authors, blats) {
  if(length(opt$pam) > 0) {
    blatPam <- gsub("[^ATCG]", "", opt$pam)
    minScore <- min(nchar(authors$exo_seq)) + nchar(blatPam)
  } else {
    minScore <- min(nchar(authors$exo_seq))
  }
  blat_command <- "blat"                                                   # external scripts and common parameters
  blat_params <- paste0("-stepSize=4 -tileSize=10 -fine -repMatch=2000000 -minScore=", minScore, " -minIdentity=100")
  
  if(file.not.exist.or.zero(blats$file_psl)) {
    dir.create(dirname(blats$file_psl), showWarnings = F)
    run_command <- paste(blat_command, blats$file_genome, blats$file_sgRNAs, blats$file_psl, blat_params)
    log_info("Sending to BLAT: ", run_command, " ...")
    system(run_command)
  } else {
    log_info("Not running BLAT: results already exist, ", blats$file_psl, ".")
  }
  
  return()
}

# Convert psl to genomic ranges
runPtgr <- function(blats) {
  
  if(file.not.exist.or.zero(blats$file_genomic_ranges)) {
    log_info("Finding alignment coordinates...")
    infile <- blats$file_psl
    outfile <- blats$file_genomic_ranges
    
    psl <- suppressWarnings(fread(infile))
    if (nrow(psl) < 6) {
      log_error("No alignments found between ", opt$infile, " and ", opt$genome, ". Did you specify the correct --guide for the --infile? Did you specify the correct --genome? Quitting.")
      stop("FATAL: Quitting due to unrecoverable error.", call. = F)
    }
    
    psl <- fread(infile, skip=5)                                                                                                                # ignore first 5 rows with nothing in them
    psl <- lapply(psl, function(x) gsub(",$", "", x))                                                                                           # fix trailing commas
    psl <- as_tibble(psl)
    names(psl) <- c("match", "mismatch", "rep. match", "N's", "Q gap count", "Q gap bases", "T gap count", "T gap bases", "strand", "Q name",   # fix names
                    "Q size", "Q start", "Q end", "T name", "T size", "T start", "T end", "block count", "blockSizes", "qStarts", "tStarts")
    psl <- psl %>%
      filter(`Q end` == `Q size`) %>% filter(`block count` == 1) %>% filter(`Q start` == 0)                                                     # we're only interested in perfect matches
    ranges <- GRanges(seqnames = psl$`T name`,                                                                                                  # verify valid genomic ranges by making a GRanges object
                      ranges = IRanges(start = as.numeric(psl$`T start`),
                                       end = as.numeric(psl$`T end`)),
                      strand = psl$strand)
    out <- tibble(seqnames = as.character(rep(ranges@seqnames@values, ranges@seqnames@lengths)),                                                # evaluate RLE to get the seqnames
                  guideBegin = ranges@ranges@start,
                  guideFinal = ranges@ranges@start + ranges@ranges@width - 1,                                                                          # end is start + width - 1
                  strand = rep(ranges@strand@values, ranges@strand@lengths),
                  assembly = blats$file_genome)
    out <- out %>% transmute(seqnames, guideBegin, guideFinal, strand, assembly,
                             start = case_when(strand == "+" ~ guideFinal - (3 + nchar(opt$pam)), # -3 to -4 upstream from PAM is the cut site
                                               strand == "-" ~ guideBegin + (3 + nchar(opt$pam))),
                             end = case_when(strand == "+" ~ guideFinal - (3 + nchar(opt$pam)),
                                             strand == "-" ~ guideBegin + (3 + nchar(opt$pam))),
                             seqSpec = case_when(strand == "+" ~ paste0(seqnames, ":", guideBegin, "-", guideFinal - nchar(opt$pam)),
                                                 strand == "-" ~ paste0(seqnames, ":", guideBegin + nchar(opt$pam), "-", (guideFinal))),
                             exo_target = paste0(seqnames, ":", guideBegin, "-", guideFinal, "_", strand),
                             exo_cut = paste0(seqnames, ":", start))
    
    
    fwrite(as.list(out$seqSpec), blats$file_genomic_seqSpecs, sep = "\n", col.names = F)
    log_info("Verifying sequence of genome hits...")
    system(paste0("twoBitToFa ", blats$file_genome, " -seqList=", blats$file_genomic_seqSpecs, " ", blats$file_genomic_seqs))
    
    seq <- fread(blats$file_genomic_seqs, header = F)
    seq <- tibble(exo_seq = seq %>% filter(!grepl("^>", V1)) %>% unlist())
    
    revcom <- function(x) {
      r <- foreach(i = x, .combine = "c") %do% {
        y = chartr("ACGT", "TGCA", toupper(i))
        y = intToUtf8(rev(utf8ToInt(y)))
      }
      return(r)
    }
    
    out <- out %>% mutate(exo_seq = seq$exo_seq,
                          exo_seq = case_when(strand == "+" ~ toupper(exo_seq),
                                              strand == "-" ~ revcom(exo_seq)))
    
    blatIn <- fread(blats$file_sgRNAs, header = F) %>%
      filter(!grepl("^>", V1)) %>%
      transmute(seq = sub(paste0(opt$pam, "$"), "", V1))
    
    nFoundByBlatNotInput <- length(which(!(unique(out$exo_seq) %in% unique(blatIn$seq))))
    nBlatIn <- length(unique(blatIn$seq))
    nFound <- length(unique(out$exo_seq))
    nFoundByBlat <- nFound - nFoundByBlatNotInput
    
    
    log_info("BLAT aligned ", nFoundByBlat, " of ", nBlatIn, " sequences given.")
    if (nFoundByBlatNotInput > 0) {
      log_warn("BLAT results ", blats$file_psl, " aligned ", nFoundByBlatNotInput," sequences that were not found in the infile. Did you checkpoint with the correct .psl file? Continuing.")
    }
    
    dir.create(dirname(outfile), recursive = T, showWarnings = F)							
    fwrite(out, outfile, sep = "\t") 
  } else {
    log_info("Not recalculating genomic hit coordinates: results already exist, ", blats$file_genomic_ranges, ".")
  }
  
  return()
}

# Convert genomic ranges to exon hits
inferExomeCols <- function(file_exons) {
  exome_headers <- fread(file_exons, nrows = 0) %>% names() %>% tolower()
  exome_cols <- list()
  exome_cols$chr <- grep("chr|seqname", exome_headers)
  exome_cols$start <- grep("start", exome_headers)
  exome_cols$end <- grep("end", exome_headers)
  exome_cols$strand <- "*"                                                      # uncomment to search both strands (default)
  # exome_cols$strand <- grep("str", exome_headers)                             # uncomment to search sense strand only
  exome_cols$symbol <- grep("symbol|name", exome_headers)
  return(exome_cols)
}

inferHitsCols <- function(file_genomic_ranges) {
  hits_headers <- fread(file_genomic_ranges, nrows = 0) %>% names() %>% tolower()
  hits_cols <- list()
  hits_cols$chr <- grep("chr|seqname", hits_headers)
  hits_cols$start <- grep("start", hits_headers)
  hits_cols$end <- grep("end", hits_headers)
  hits_cols$strand <- grep("str", hits_headers)
  hits_cols$exo_seq <- grep("^exo_seq$", hits_headers)
  return(hits_cols)
}

exon_melter <- function(exons, exome_cols) {
  exons[[exome_cols$start]] <- strsplit(exons[[exome_cols$start]], ",")
  exons[[exome_cols$end]] <- strsplit(exons[[exome_cols$end]], ",")
  exons <- exons %>%
    unnest(cols = c(all_of(exome_cols$start), all_of(exome_cols$end))) %>%
    unique()
  exons[[exome_cols$start]] <- as.numeric(exons[[exome_cols$start]])
  exons[[exome_cols$end]] <- as.numeric(exons[[exome_cols$end]])
  return(exons)
}

import_exome <- function(file_exons, exome_cols) {
  
  exome <- fread(file_exons)
  if (any(grepl(",", exome[[exome_cols$start]]))) { # check if we need to melt exons
    exome <- exon_melter(exome, exome_cols)
  }
  seqnames <- exome[[as.numeric(exome_cols$chr)]]
  start <- exome[[as.numeric(exome_cols$start)]]
  end <- exome[[as.numeric(exome_cols$end)]]
  if (exome_cols$strand == "*") {                    # use * if in strandless mode
    strand <- rep("*", nrow(exome))
  } else {
    strand <- exome[[as.numeric(exome_cols$strand)]]
  }
  symbol <- exome[[as.numeric(exome_cols$symbol)]] %>% gsub(pattern = ",$", replacement = "")
  
  exome <- tibble(seqnames = seqnames, start = start, end = end, strand = strand, exo_symbol = symbol) %>%
    filter(!is.na(seqnames) & !is.na(start) & !is.na(end) & end >= start) %>%
    mutate(ranges = paste(start, end, sep = "-"))
  exome <- GRanges(seqnames = exome$seqnames, ranges = exome$ranges, strand = exome$strand, exo_symbol = exome$exo_symbol)
  return(exome)
}

import_hits <- function(file_genomic_ranges, hits_cols) {
  
  hits <- fread(file_genomic_ranges)
  seqnames <- hits[[as.numeric(hits_cols$chr)]]
  start <- hits[[as.numeric(hits_cols$start)]]
  end <- hits[[as.numeric(hits_cols$end)]]
  exo_seq = hits[[as.numeric(hits_cols$exo_seq)]]
  if (hits_cols$strand == "*") {                    # use * if in strandless mode
    strand <- rep("*", nrow(hits))
  } else {
    strand <- hits[[as.numeric(hits_cols$strand)]]
  }
  
  hits <- tibble(seqnames = seqnames, start = start, end = end, strand = strand, exo_seq = exo_seq) %>%
    filter(!is.na(seqnames) & !is.na(start) & !is.na(end) & end >= start) %>%
    mutate(ranges = paste(start, end, sep = "-"))
  hits <- GRanges(seqnames = hits$seqnames, ranges = hits$ranges, strand = hits$strand, exo_seq = as.character(hits$exo_seq))
  return(hits)
}

runGrtem <- function(blats) {
  
  if(file.not.exist.or.zero(blats$file_genomic_ranges_matched) | file.not.exist.or.zero(blats$file_genomic_ranges_distances)) {
    log_info("Determining exonic hits... ", blats$file_exons, "...")
    
    exome_cols <- inferExomeCols(blats$file_exons)                         # Infer column identities in the exome file
    hits_cols <- inferHitsCols(blats$file_genomic_ranges)                  # Infer column identities in the hits file
    
    exons <- import_exome(blats$file_exons, exome_cols)                    # Make a gRanges object from the exome file (with Symbol column to provide re-annotations)
    hits <- import_hits(blats$file_genomic_ranges, hits_cols)              # Make a gRanges object from the hits file (with ID column indicating guides in the library to be reannotated)
    
    exon_hits <- suppressWarnings(findOverlaps(hits, exons))                  # Make a gRanges Hits object indicating the pairs of gRanges that overlap between query (hits) and subject (exome)
    exon_hitsRanges <- suppressWarnings(findOverlapPairs(hits, exons))        # Make a gRanges Pairs object indicating the genomic ranges of pairs of gRanges that overlap between the first (hits) and second (exome)
    sgrna_hits <- GRanges(exon_hitsRanges@first,                              # Make a gRanges object which contains the genomic ranges of the hits and the Symbols from the exome
                          exo_symbol = exons$exo_symbol[exon_hits@to],
                          exo_seq = as.character(hits$exo_seq[exon_hits@from]))
    mapping <- as_tibble(sgrna_hits) %>%
      mutate(exo_cut = paste0(seqnames, ":", start)) %>%
      unique()
    mapping$assembly <- blats$file_genome
    mapping$exome <- blats$file_exons
    
    distances <- suppressWarnings(distanceToNearest(hits, exons, select = "arbitrary")) # Make a gRanges Hits object showing the distances from each hit (queryHits) to the nearest exon in the exome (subjectHits)
    sgrna_distances <- GRanges(hits[distances@from],                          # Make a gRanges object annotated with the gene of the nearest exon and the distance to that exon
                               nearestGene = exons$exo_symbol[distances@to],
                               distance = distances@elementMetadata$distance,
                               exo_seq = as.character(hits$exo_seq[distances@from]))
    
    final_distances <- as_tibble(sgrna_distances) %>% 
      mutate(exo_cut = paste0(seqnames, ":", start)) %>%
      unique()
    final_distances$assembly <- blats$file_genome
    final_distances$exome <- blats$file_exons
    
    dir.create(dirname(blats$file_genomic_ranges_matched), recursive = T, showWarnings = F)
    dir.create(dirname(blats$file_genomic_ranges_distances), recursive = T, showWarnings = F)
    fwrite(mapping, blats$file_genomic_ranges_matched, sep = "\t")
    fwrite(final_distances, blats$file_genomic_ranges_distances, sep = "\t")
  } else {
    log_info("Not recalculating exonic hits: results already exist, ", blats$file_genomic_ranges_matched, ".")
  }
  return()
}

# Write master mapping file
exorcisemaster <- function(premaster, opt) {
  
  
  for (i in 1:length(premaster)) {
    premaster[[i]][which(is.na(premaster[[i]]) | premaster[[i]] == "" | grepl("^exo_Non-targeting_", premaster[[i]]))] <- "X"
  }
  
  if(is.null(opt$control)) {
    ncontrols <- length(premaster$exo_symbol[which(premaster$exo_symbol == "X")])
    premaster$exo_symbol[which(premaster$exo_symbol == "X")] <- paste0("exo_Non-targeting", 1:ncontrols)
  }
  
  # Fix controls
  if (!is.null(opt$control)) {
    for (s in 1:length(opt$control)) {
      log_info("Finding control sequences: replacing ", opt$control[s], " with ", opt$control_type[s], ".")
      control_guides <- rep(F, nrow(premaster))
      if ("exo_orig" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control[s], premaster$exo_orig)  # search for the control string in authors' symbols and IDs if those columns exist
      
      nThisControl <- length(premaster$exo_symbol[control_guides & premaster$exo_symbol == "X"])
      premaster$exo_symbol[control_guides & premaster$exo_symbol == "X"] <- paste0(opt$control_type[s], "_", 1:nThisControl)   # Map control guides to control annotations unless there is an approved Symbol column
    }
    
    ncontrols <- length(premaster$exo_symbol[which(premaster$exo_symbol == "X")])
    premaster$exo_symbol[which(premaster$exo_symbol == "X")] <- paste0("exo_Non-targeting_", 1:ncontrols) # Catch the remaining non-targeting guides
  }
  
  master <- premaster %>% unique()
  
  # Harmonise
  
  
  if ("exo_orig" %in% names(master)) {
    simple <- exorciseinferTargets(master, opt)
    master <- left_join(master, simple, by = "exo_orig") %>%
      mutate(exo_harm = case_when(is.na(exo_harm) ~ exo_symbol, # Accept harmonised symbol if found
                                  T ~ exo_harm))                # Accept exorcised symbol if no harmonised symbol
    
    if (length(opt$control) > 0) { # For control guides, retain control types in the harmonised column, ie. without exorcism
      for (i in 1:length(opt$control)) {
        thisControlPattern <- opt$control[i]
        thisControlType <- opt$control_type[i]
        master <- master %>%
          mutate(exo_harm = case_when(grepl(thisControlPattern, exo_orig) ~ paste0(thisControlType, "_", 1:n()),
                                      T ~ exo_harm))
      }
    }
    
    master <- master %>% relocate(exo_id, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut)
  } else {
    master <- master %>% relocate(exo_id, exo_seq, exo_symbol, exo_target, exo_cut)
  }
  
  master <- master %>% mutate(exo_id = paste0("exorcise_", exo_seq, "_", exo_target, "_", exo_symbol))
  
  
  return(master)
}

# Infer intended target per authors' original symbols and write
exorciseinferTargets <- function(master, opt) {
  
  simple <- master %>%
    dplyr::select(exo_symbol, exo_orig) %>% unique()
  
  # Load in gene symbol annotations for ranking. This file contains all possible approved symbols and gene classes
  annot <- fread(opt$priorities)
  names(annot) <- c("Symbol", "Gene Type")
  annot$`Gene Type` <- factor(annot$`Gene Type`, ordered = T, levels = c("PROTEIN_CODING", "ncRNA", "PSEUDO", "rRNA", "tRNA", "snRNA", "snoRNA", "scRNA", "BIOLOGICAL_REGION", "OTHER", "")) 
  
  # Infer
  
  simple_dup <- simple %>%
    left_join(annot, by = c("exo_symbol" = "Symbol"))                  # Join the multi-mapped original symbols to annotations for their mapping target
  simple_dup$`Gene Type`[which(is.na(simple_dup$`Gene Type`))] <- "OTHER"
  
  simple_dup2 <- simple_dup %>% filter(!grepl(paste0(paste0("(^", c("X", opt$control_type), "\\d+$)"), collapse = "|"), exo_symbol)) # Remove controls and unmapped guides when considering intended target
  
  dup <- unique(simple_dup2$exo_orig)
  
  # Fix multi-mapped original symbols
  j <- 1
  tot <- length(dup) 
  log_info("Harmonising... ")
  
  simple_fixed <- data.table(exo_harm = rep("", tot),
                             exo_orig = rep("", tot))        # Initialise a tibble to contain fixed, singly-mapped original symbols
  
  for (s in dup) {  
    q <- simple_dup2 %>% filter(exo_orig == s)
    cons <- q %>%                                                                         # Determine the candidates that appeared the most often
      group_by(exo_symbol) %>%
      summarise(consensus = n()) %>%
      mutate(consensus = case_when(consensus == max(consensus) ~ 1,
                                   T ~ 0))
    q <- q %>% left_join(cons, by = "exo_symbol") %>% filter(consensus == 1)                  # Eliminate candidates which did not appear the most times
    q <- q %>% mutate(fit_rank = case_when(exo_symbol == exo_orig ~ 0,                 # Rank remaining, tied 1st place candidates: identical names > protein-coding > ncRNA > pseudogene > others
                                           T ~ as.numeric(`Gene Type`)))
    accept <- q %>% filter(fit_rank == min(fit_rank)) %>% arrange(exo_symbol) %>% head(1)     # If the top rank is still tied, accept the first match alphabetically
    set(simple_fixed, as.integer(j), "exo_harm", accept$exo_symbol[1])
    set(simple_fixed, as.integer(j), "exo_orig", accept$exo_orig[1])
    j <- j + 1
    if(j %% 1000 == 0) {
      log_info("Harmonising... ", j, " out of ", tot, " groups done...")
    }
  }
  log_info("Harmonising... ", tot, " out of ", tot, " groups done.")
  
  # Combine
  simple <- simple_fixed %>% unique()
  
  return (simple)
  
}

reannotateExisting <- function(opt) {
  
  input <- fread(opt$infile) %>%
    relocate(exo_seq = opt$seq) %>%
    unique()
  
  library <- fread(opt$library)
  
  if(!is.null(opt$harm)) {
    input <- input %>% relocate(exo_orig = opt$harm, exo_seq)
    library <- library %>% select(exo_id, exo_seq, exo_symbol, exo_target, exo_cut) %>% unique()                         # harm specified: harmonise
    exorcised <- left_join(input, library, by = "exo_seq")
    exorcised <- exorcisemaster(exorcised, opt)
  } else {
    if("exo_harm" %in% names(library)) {
      log_info("Inheriting existing harmonisations from the library...")
      library <- library %>% select(exo_id, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut) %>% unique()   # harm not specified: inherit harmonisations if exist
    } else {
      library <- library %>% select(exo_id, exo_seq, exo_symbol, exo_target, exo_cut) %>% unique()
    }
    exorcised <- left_join(input, library, by = "exo_seq")
  }
  
  if("exo_harm" %in% names(exorcised)) {
    exorcised <- exorcised %>% relocate(exo_id, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut)
  } else {
    exorcised <- exorcised %>% relocate(exo_id, exo_seq, exo_symbol, exo_target, exo_cut)
  }
  
  fwrite(exorcised, paste0(opt$outdir, "/exorcise.tsv"), sep = "\t")
  
}

fixOpts <- function(opt) {
  
  # fix comma-separated argument inputs
  if(!is.null(opt$infile))       opt$infile       <- commasplit(opt$infile)
  if(!is.null(opt$outdir))       opt$outdir       <- commasplit(opt$outdir)
  if(!is.null(opt$seq))          opt$seq          <- commasplit(opt$seq)
  if(!is.null(opt$pam))          opt$pam          <- commasplit(opt$pam)
  #if(!is.null(opt$mode))         opt$mode         <- commasplit(opt$mode)
  if(!is.null(opt$library))      opt$library      <- commasplit(opt$library)
  if(!is.null(opt$genome))       opt$genome       <- commasplit(opt$genome)
  if(!is.null(opt$exome))        opt$exome        <- commasplit(opt$exome)
  if(!is.null(opt$priorities))   opt$priorities   <- commasplit(opt$priorities)
  #if(!is.null(opt$id))           opt$id           <- commasplit(opt$id)
  if(!is.null(opt$harm))         opt$harm         <- commasplit(opt$harm)
  if(!is.null(opt$control))      opt$control      <- commasplit(opt$control)
  if(!is.null(opt$control_type)) opt$control_type <- commasplit(opt$control_type)
  
  opt$adhoc <- F
  
  # check arguments
  errors <- list()
  warnings <- list()
  
  # check if number of inputs are correct
  lengthChecker <- function(got, need, label) {
    got <- length(got)
    need <- length(need)
    if(got != 0 & got != 1 & got != need) {
      return(paste0("Invalid number of inputs to --", label, ". Should be 1 or ", need, ", but got ", got, ". Using first value only."))
    } else {
      return(NULL)
    }
  }
  
  # check infile
  if(length(opt$infile) > 0) {
    if(length(opt$infile) > 1) {
      opt$infile <- opt$infile[1]
      warnings <- c(warnings, paste0("Warning: --infile received more than one argument. Using first value only: ", opt$infile, "."))
    }
    if(!file.exists(opt$infile)) {
      opt$infile_glob <- Sys.glob(paste0(opt$infile, "*"))
      if(length(opt$infile_glob) == 0) {
        errors <- c(errors, paste0("Error: --infile ", opt$infile, " not found."))
        opt$infile <- NULL
      } else {
        opt$infile <- opt$infile_glob
        warnings <- c(warnings, paste0("Warning: --infile accepted as globbed argument ", opt$infile, "."))
      }
    }
    if(!isFile(opt$infile)) {
      errors <- c(errors, paste0("Error: --infile ", opt$infile, " is not a file."))
      opt$infile <- NULL
    }
  } else {
    errors <- c(errors, paste0("Error: --infile not specified."))
  }
  
  # check infile headers
  if(length(opt$infile) > 0) {
    opt$infile_headers <- names(fread(opt$infile, nrows = 1))
    if (any(grepl("^exo_", opt$infile_headers))) {
      warnings <- c(warnings, paste0("Warning: --infile ", opt$infile, " contains exorcise-like column names (\"exo_\"). These might be clobbered."))
    }
  }
  
  # check outdir
  if(length(opt$outdir) > 0) {
    if(length(opt$outdir) > 1) {
      opt$outdir <- opt$outdir[1]
      warnings <- c(warnings, paste0("Warning: --outdir received more than one argument. Using first value only: ", opt$outdir, "."))
    }
    if(!dir.exists(opt$outdir)) {
      dir.create(opt$outdir, showWarnings = F, recursive = T)
      warnings <- c(warnings, paste0("Info: --outdir ", opt$outdir, " being recursively created."))
    }
  } else {
    errors <- c(errors, paste0("Error: --outdir not specified."))
  }
  
  # check seq
  if(length(opt$seq) > 0) {
    if(length(opt$seq) > 1) {
      opt$seq <- opt$seq[1]
      warnings <- c(warnings, paste0("Warning: --seq received more than one argument. Using first value only: ", opt$seq, "."))
    }
    if(!grepl("^\\d+$", opt$seq)) {
      errors <- c(errors, paste0("Error: --seq must be an integer, got ", opt$seq, "."))
    } else {
      opt$seq <- as.numeric(opt$seq)
      if (length(opt$infile) > 0) {
        if (opt$seq > length(opt$infile_headers)) {
          errors <- c(errors, paste0("Error: --seq is out of bounds, got ", opt$seq, " but --infile only contains ", length(opt$infile_headers), " columns."))
        }
      }
    }
  } else {
    errors <- c(errors, paste0("Error: --seq not specified."))
  }
  
  # check library
  if(length(opt$library) > 0) {
    if(length(opt$library) > 1) {
      opt$library <- opt$library[1]
      warnings <- c(warnings, paste0("Warning: --library received more than one argument. Using first value only: ", opt$library, "."))
    }
    if(!file.exists(opt$library)) {
      opt$library_glob <- Sys.glob(paste0(opt$library, "*"))
      if(length(opt$library_glob) == 0) {
        errors <- c(errors, paste0("Error: --library ", opt$library, " not found."))
      } else {
        opt$library <- opt$library_glob
        warnings <- c(warnings, paste0("Warning: --library accepted as globbed argument ", opt$library))
      }
    }
    if(file.exists(opt$library)) {
      opt$library_headers <- names(fread(opt$library, nrows = 0))
      if(!("exo_id" %in% opt$library_headers)) {
        errors <- c(errors, paste0("Error: --library ", opt$library, " doesn't look like an exorcised library. Expected an `exo_id` column."))
      }
      if(!("exo_seq" %in% opt$library_headers)) {
        errors <- c(errors, paste0("Error: --library ", opt$library, " doesn't look like an exorcised library. Expected an `exo_seq` column."))
      }
      if(!("exo_symbol" %in% opt$library_headers)) {
        errors <- c(errors, paste0("Error: --library ", opt$library, " doesn't look like an exorcised library. Expected an `exo_symbol` column."))
      }
    }
  } else {
    opt$adhoc <- T
    warnings <- c(warnings, paste0("Info: --library not specified. Using ad-hoc mode."))
  }
  
  # check pam
  if(length(opt$pam) > 0) {
    if(opt$adhoc) {
      if(length(opt$pam) > 1) {
        opt$pam <- opt$pam[1]
        warnings <- c(warnings, paste0("Warning: --pam received more than one argument. Using first value only: ", opt$pam, "."))
      }
      if(grepl("[^ACTGNactgn]", opt$pam)) {
        errors <- c(errors, paste0("Error: --pam ", opt$pam," contains non-nucleotide letters. Only [ACTGN] are accepted."))
      }
    } else {
      warnings <- c(warnings, paste0("Warning: --pam specified outside of ad-hoc mode. Ignoring"))
    }
  } else if(opt$adhoc) {
    opt$pam <- ""
    warnings <- c(warnings, paste0("Warning: --pam not specified. Not appending any PAM."))
  }
  
  # check mode
  # if(length(opt$mode) > 0) {
  #   if(length(opt$mode) > 1) {
  #     opt$mode <- opt$mode[1]
  #     warnings <- c(warnings, paste0("Warning: --mode received more than one argument. Using first value only: ", opt$mode, "."))
  #   }
  #   opt$mode <- toupper(opt$mode)
  #   if(!(opt$mode %in% c("KO", "A", "I"))) {
  #     warnings <- c(warnings, paste0("Warning: --mode received an invalid value: ", opt$mode, ". Falling back to KO."))
  #     opt$mode <- "KO"
  #   }
  # }
  
  # check genome
  if(length(opt$genome) > 0) {
    if(length(opt$genome) > 1) {
      opt$genome <- opt$genome[1]
      warnings <- c(warnings, paste0("Warning: --genome received more than one argument. Using first value only: ", opt$genome, "."))
    }
    if(!file.exists(opt$genome)) {
      opt$genome_glob <- Sys.glob(paste0(opt$genome, "*"))
      if(length(opt$genome_glob) == 0) {
        errors <- c(errors, paste0("Error: --genome ", opt$genome, " not found."))
      } else {
        opt$genome <- opt$genome_glob
        warnings <- c(warnings, paste0("Warning: --genome accepted as globbed argument ", opt$genome))
      }
    }
  } else if(opt$adhoc) {
    errors <- c(errors, paste0("Error: --genome not passed while in ad-hoc mode."))
  }
  
  # check exome
  if(length(opt$exome) > 0) {
    if(length(opt$exome) > 1) {
      opt$exome <- opt$exome[1]
      warnings <- c(warnings, paste0("Warning: --exome received more than one argument. Using first value only: ", opt$exome, "."))
    }
    if(!file.exists(opt$exome)) {
      opt$exome_glob <- Sys.glob(paste0(opt$exome, "*"))
      if(length(opt$exome_glob) == 0) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " not found."))
      } else {
        opt$exome <- opt$exome_glob
        warnings <- c(warnings, paste0("Warning: --exome accepted as globbed argument ", opt$exome))
      }
    }
    if(file.exists(opt$exome)) {
      opt$exome_headers <- names(fread(opt$exome, nrows = 0))
      if(!("#chrom" %in% opt$exome_headers)) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected a `#chrom` column."))
      }
      if(!("strand" %in% opt$exome_headers)) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected a `strand` column."))
      }
      if(!("exonStarts" %in% opt$exome_headers)) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected an `exonStarts` column."))
      }
      if(!("exonEnds" %in% opt$exome_headers)) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected an `exonEnds` column."))
      }
      if(!("name2" %in% opt$exome_headers)) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected a `name2` column."))
      }
    }
  } else if(opt$adhoc) {
    errors <- c(errors, paste0("Error: --exome not passed while in ad-hoc mode."))
  }
  
  # check id
  # if(length(opt$id) > 0) {
  #   if(opt$adhoc) {
  #     if(length(opt$id) > 1) {
  #       opt$id <- opt$id[1]
  #       warnings <- c(warnings, paste0("Warning: --id received more than one argument. Using first value only: ", opt$id, "."))
  #     }
  #     if(!grepl("^\\d+$", opt$id)) {
  #       errors <- c(errors, paste0("Error: --id must be an integer, got ", opt$id, "."))
  #     } else {
  #       opt$id <- as.numeric(opt$id)
  #       if (length(opt$infile) > 0) {
  #         if (opt$id > length(opt$infile_headers)) {
  #           errors <- c(errors, paste0("Error: --id is out of bounds, got ", opt$id, " but --infile only contains ", length(opt$infile_headers), " columns."))
  #         }
  #       }
  #     }
  #   } else {
  #     warnings <- c(warnings, paste0("Warning: --id specified outside of ad-hoc mode. Ignoring."))
  #   }
  # }
  
  # check harm
  if(length(opt$harm) > 0) {
    if(length(opt$harm) > 1) {
      opt$harm <- opt$harm[1]
      warnings <- c(warnings, paste0("Warning: --harm received more than one argument. Using first value only: ", opt$harm, "."))
    }
    if(!grepl("^\\d+$", opt$harm)) {
      errors <- c(errors, paste0("Error: --harm must be an integer, got ", opt$harm, "."))
    } else {
      opt$harm <- as.numeric(opt$harm)
      if (length(opt$infile) > 0) {
        if (opt$harm > length(opt$infile_headers)) {
          errors <- c(errors, paste0("Error: --harm is out of bounds, got ", opt$harm, " but --infile only contains ", length(opt$infile_headers), " columns."))
        }
      }
    }
  }
  
  # check priorities
  if(length(opt$priorities) > 0) {
    if(length(opt$priorities) > 1) {
      opt$priorities <- opt$priorities[1]
      warnings <- c(warnings, paste0("Warning: --priorities received more than one argument. Using first value only: ", opt$priorities, "."))
    }
    if(!file.exists(opt$priorities)) {
      opt$priorities_glob <- Sys.glob(paste0(opt$priorities, "*"))
      if(length(opt$priorities_glob) == 0) {
        errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " not found."))
      } else {
        opt$priorities <- opt$priorities_glob
        warnings <- c(warnings, paste0("Warning: --priorities accepted as globbed argument ", opt$priorities))
      }
    }
    if(file.exists(opt$priorities)) {
      opt$priorities_headers <- names(fread(opt$priorities, nrows = 0))
      if(!("Symbol" %in% opt$priorities_headers)) {
        errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " doesn't look like a feature priorities file. Expected a `Symbol` column."))
      }
      if(!("Gene Type" %in% opt$priorities_headers)) {
        errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " doesn't look like a feature priorities file. Expected a `Gene Type` column."))
      }
    }
    if(length(opt$harm) == 0) {
      warnings <- c(warnings, paste0("Warning: --priorities passed without --harm. Ignoring."))
    }
  } else if(length(opt$harm) > 0) {
    if(opt$harm != 0) {
      errors <- c(errors, paste0("Error: --priorities not passed while --harm passed."))
    }
  }
  
  # check control
  if(length(opt$control) > 0) {
    if(length(opt$harm) == 0) {
      opt$control <- NULL
      warnings <- c(warnings, paste0("Warning: --control passed without --harm. Ignoring."))
    }
  }
  
  # check control_type
  if(length(opt$control_type) > 0) {
    if(length(opt$control) == 0) {
      opt$control_type <- NULL
      warnings <- c(warnings, paste0("Warning: --control_types passed without --control. Ignoring."))
    } else if(length(opt$control) == 1) {
      if(length(opt$control_type) > 1) {
        opt$control_type <- opt$control_type[1]
        warnings <- c(warnings, paste0("Warning: --control_type received more than one argument. Using first value only: ", opt$control_type, "."))
      }
    } else if(length(opt$control) > 1) {
      if(length(opt$control_type) == 1) {
        opt$control_type <- rep(opt$control_type, length(opt$control))
      } else if(length(opt$control_type) != length(opt$control)) {
        errors <- c(errors, paste0("Error: --control and --control_type lengths are not equal."))
      }
    }
    if("exo_Non-targeting" %in% opt$control_type) {
      errors <- c(errors, paste0("Error: --control_type contains the disallowed value \"exo_Non-targeting\"."))
    }
    if(any(duplicated(opt$control_type))) {
      errors <- c(errors, paste0("Error: --control_type contains duplicated values."))
    }
  } else if(length(opt$control) > 0) {
    opt$control_type <- paste0("Non-targeting", 1:length(opt$control))
    warnings <- c(warnings, paste0("Warning: --control passed without --control_type. Assuming --control_type for ", opt$control, " is ", opt$control_type, "."))
  }
  
  # check ref
  if(length(opt$ref) > 0) {
    ref = paste0("
 exorcise version ", ver, " was developed by Dr Simon Lam, University of Cambridge
 https://github.com/SimonLammmm/exorcise

 If you found exorcise useful in your work, please cite:
 Lam S, exorcise [https://github.com/SimonLammmm/exorcise] (manuscript under preparation).

 If you used exorcise in ad-hoc mode, please cite:
 Kent WJ, 2002, BLAT--the BLAST-like alignment tool, Genome Res 12(4): 656-664, doi: 10.1101/gr.229202
 and the source of your genome assembly and exome annotations.

=======================================================================================================")
    cat(ref, "\n\n")
    options(show.error.messages = FALSE)
    stop()
  }
  
  # remove temporary options
  opt$infile_glob <- NULL
  opt$infile_headers <- NULL
  opt$library_glob <- NULL
  opt$library_headers <- NULL
  opt$genome_glob <- NULL
  opt$exome_glob <- NULL
  opt$exome_headers <- NULL
  opt$priorities_glob <- NULL
  opt$priorities_headers <- NULL
  
  # print warnings and errors
  if(length(warnings) > 0) {
    for (i in 1:length(warnings)) {
      log_warn(warnings[[i]])
    }
  }
  
  if(length(errors) > 0) {
    for (i in 1:length(errors)) {
      log_error(errors[[i]])
    }
    if(!is.null(opt$outdir)) {
      stop("FATAL: There were errors in the input. Quitting. Please check the logfile: ", logfile, ".")
    } else {
      stop("FATAL: There were errors in the input. Quitting.", call. = F)
    }
  }
  
  # if(length(warnings) > 0 & opt$just_go == F) {
  #   readline(prompt = "Accept these settings? Press [enter] to continue.")
  # }
  
  return(opt)
}

#### MAIN ####

# Test vector
if (interactive()) {
  opt <- list()
  opt$infile <-       NULL
  opt$outdir <-       NULL
  opt$seq <-          NULL
  opt$pam <-          NULL
  opt$genome <-       NULL
  opt$exome <-        NULL
  opt$priorities <-   NULL
  opt$harm <-         NULL
  opt$control <-      NULL
  opt$control_type <- NULL
  opt$library <-      NULL
  opt$ref <-          NULL
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
    # make_option(opt_str = c("-q", "--mode"), type = "character", default = NULL,
    #             help = "CRISPR screen type: KO (knockout), a (activation), i (inhibition), BE (base editing)", metavar = "character"),
    make_option(opt_str = c("-l", "--library"), type = "character", default = NULL,
                help = "(required if --genome, --exome, and --priorities not specified) Exorcised file to be used as library.", metavar = "character"),
    make_option(opt_str = c("-v", "--genome"), type = "character", default = NULL,
                help = "(required if --library not specified) 2bit genome.", metavar = "character"),
    make_option(opt_str = c("-w", "--exome"), type = "character", default = NULL,
                help = "(required if --library not specified) Exome.", metavar = "character"),
    # make_option(opt_str = c("-j", "--id"), type = "character", default = NULL,
    #             help = "(optional, ignored if --library specified) ID column number.", metavar = "character"),
    make_option(opt_str = c("-n", "--harm"), type = "character", default = NULL,
                help = "(optional) Existing annotation column number.", metavar = "character"),
    make_option(opt_str = c("-y", "--priorities"), type = "character", default = NULL,
                help = "(required if --harm specified) Priorities file.", metavar = "character"),
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
