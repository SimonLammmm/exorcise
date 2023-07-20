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
# 0.9           2022-12-09T12-28-16   evaluation of full release


ver <- 0.9

#### INIT ####
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(readxl)
  library(tidyr)
  library(foreach)
  library(logger)
  library(R.utils)
  library(GenomicRanges)
})


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
    authors <- importAuthorsLib(opt)
    blats <- list(file_genome = opt$genome,
                  file_exons = opt$exome,
                  file_feature_priorities = opt$priorities,
                  file_sgRNAs = paste0(opt$outdir, "/exorcise.1-seq.fa"),
                  file_psl = paste0(opt$outdir, "/exorcise.2-",  sub(".+/(.+?)$", "\\1", opt$genome), "_BLAT.psl"),
                  file_genomic_ranges = paste0(opt$outdir, "/exorcise.3-", sub(".+/(.+?)$", "\\1", opt$genome), "_genomicRanges.tsv"),
                  file_genomic_ranges_matched = paste0(opt$outdir, "/exorcise.4-", sub(".+/(.+?)$", "\\1", opt$exome), "_exonHits.tsv"),
                  file_genomic_ranges_distances = paste0(opt$outdir, "/exorcise.5-", sub(".+/(.+?)$", "\\1", opt$exome), "_exonDist.tsv"),
                  file_exorcise_master_out = paste0(opt$outdir, "/exorcise.tsv"))
  
    extractGuides(opt, authors, blats)
    runBlat(opt, authors, blats)
    runPtgr(blats)
    runGrtem(blats)
    
     # generate mapping
      genome_hits <- fread(blats$file_genomic_ranges, colClasses = "character") %>% transmute(exo_id, exo_target = paste0(seqnames, ":", guideBegin, "-", guideFinal, "_", strand), exo_cut = paste0(seqnames, ":", start))
      exome_hits <- fread(blats$file_genomic_ranges_matched, colClasses = "character") %>% transmute(exo_id, exo_symbol)
      all_mappings <- left_join(genome_hits, exome_hits, by = c("exo_id")) %>% unique()
    
    premaster <- left_join(authors, all_mappings, by = "exo_id")
    master <- exorcisemaster(premaster, blats, opt)
    
  } else {
    reannotateExisting(opt)
  }
 
}

### OTHER FUNCTIONS ####

# Import authors' original library files, extracting original symbols, detecting/creating a primary key, and removing adapters where necessary
importAuthorsLib <- function(opt) {
  log_info("Importing file ", opt$infile)
  authors <- fread(opt$infile)
  
  if(!is.null(opt$id)) {
    authors <- authors %>% mutate(exo_id = paste0("exorcise_", 1:n(), "_", .[[opt$id]], "_", .[[opt$harm]]))
  } else {
    authors <- authors %>% mutate(exo_id = paste0("exorcise_", 1:n(), "_", .[[opt$harm]]))
  }
  
  authors <- authors %>%                               # read authors' file
    as_tibble() %>%
    relocate(exo_id, exo_seq = opt$seq, exo_orig = opt$harm) %>%     # select appropriate columns
    unique()
  
  return(authors)
}

# Use authors' original library to extract guide RNA sequences for alignment
extractGuides <- function(opt, authors, blats) {
  log_info("Creating sequence FASTA")
  sgRNAs <- paste0(">", authors$exo_id, "\n", authors$exo_seq, opt$pam)
  dirwrite(list(sgRNAs), blats$file_sgRNAs, sep = "\n", quote = F)
  log_info("Wrote sgRNA FASTA to ", blats$file_sgRNAs)
  return()
}

# Run BLAT
runBlat <- function(opt, authors, blats) {
  if(length(opt$pam) == 1) {
    blatPam <- gsub("[^ATCG]", "", opt$pam)
    minScore <- min(nchar(authors$exo_seq)) + nchar(blatPam)
  } else {
    minScore <- min(nchar(authors$exo_seq))
  }
  blat_command <- "blat"                                                   # external scripts and common parameters
  blat_params <- paste0("-stepSize=4 -tileSize=10 -fine -repMatch=2000000 -minScore=", minScore, " -minIdentity=100")
  
  while(file.not.exist.or.zero(blats$file_psl)) {
      dir.create(dirname(blats$file_psl), showWarnings = F)
      run_command <- paste(blat_command, blats$file_genome, blats$file_sgRNAs, blats$file_psl, blat_params)
      log_info("Sending to BLAT: ", run_command)
      system(run_command)
  }
  
  return()
}

# Convert psl to genomic ranges
runPtgr <- function(blats) {
  
      while(file.not.exist.or.zero(blats$file_genomic_ranges)) {
      infile <- blats$file_psl
      outfile <- blats$file_genomic_ranges
      
      psl <- fread(infile, skip=5)                                                                                                                # ignore first 5 rows with nothing in them
      psl <- lapply(psl, function(x) gsub(",$", "", x))                                                                                           # fix trailing commas
      psl <- as_tibble(psl)
      names(psl) <- c("match", "mismatch", "rep. match", "N's", "Q gap count", "Q gap bases", "T gap count", "T gap bases", "strand", "Q name",   # fix names
                      "Q size", "Q start", "Q end", "T name", "T size", "T start", "T end", "block count", "blockSizes", "qStarts", "tStarts")
      log_info("Obtaining genomic hits... ", infile, "...")
      psl <- psl %>%
        filter(`Q end` == `Q size`) %>% filter(`block count` == 1) %>% filter(`Q start` == 0)                                                     # we're only interested in perfect matches
      ranges <- GRanges(seqnames = psl$`T name`,                                                                                                  # verify valid genomic ranges by making a GRanges object
                        ranges = IRanges(start = as.numeric(psl$`T start`),
                                         end = as.numeric(psl$`T end`)),
                        exo_id = psl$`Q name`,
                        strand = psl$strand)
      out <- tibble(seqnames = as.character(rep(ranges@seqnames@values, ranges@seqnames@lengths)),                                                # evaluate RLE to get the seqnames
                    guideBegin = ranges@ranges@start,
                    guideFinal = ranges@ranges@start + ranges@ranges@width - 1,                                                                          # end is start + width - 1
                    exo_id = ranges$exo_id,
                    strand = rep(ranges@strand@values, ranges@strand@lengths),
                    assembly = blats$file_genome)
      out <- out %>% transmute(seqnames, guideBegin, guideFinal, exo_id, strand, assembly,
                               start = case_when(strand == "+" ~ guideFinal - (3 + nchar(opt$pam)), # -3 to -4 upstream from PAM is the cut site
                                                 strand == "-" ~ guideBegin + (3 + nchar(opt$pam))),
                               end = case_when(strand == "+" ~ guideFinal - (3 + nchar(opt$pam)),
                                               strand == "-" ~ guideBegin + (3 + nchar(opt$pam))))
      dir.create(dirname(outfile), recursive = T, showWarnings = F)							
      fwrite(out, outfile, sep = "\t") 
      log_info("Wrote genomic hits file to ", outfile)
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
  hits_cols$id <- grep("^exo_id$", hits_headers)
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
  exo_id = hits[[as.numeric(hits_cols$id)]]
  if (hits_cols$strand == "*") {                    # use * if in strandless mode
    strand <- rep("*", nrow(hits))
  } else {
    strand <- hits[[as.numeric(hits_cols$strand)]]
  }
  
  hits <- tibble(seqnames = seqnames, start = start, end = end, strand = strand, exo_id = exo_id) %>%
    filter(!is.na(seqnames) & !is.na(start) & !is.na(end) & end >= start) %>%
    mutate(ranges = paste(start, end, sep = "-"))
  hits <- GRanges(seqnames = hits$seqnames, ranges = hits$ranges, strand = hits$strand, exo_id = as.character(hits$exo_id))
  return(hits)
}

runGrtem <- function(blats) {
  
    while(file.not.exist.or.zero(blats$file_genomic_ranges_matched) | file.not.exist.or.zero(blats$file_genomic_ranges_distances)) {
      log_info("Obtaining exonic hits... ", blats$file_exons, "...")
      
      exome_cols <- inferExomeCols(blats$file_exons)                         # Infer column identities in the exome file
      hits_cols <- inferHitsCols(blats$file_genomic_ranges)                  # Infer column identities in the hits file
      
      exons <- import_exome(blats$file_exons, exome_cols)                    # Make a gRanges object from the exome file (with Symbol column to provide re-annotations)
      hits <- import_hits(blats$file_genomic_ranges, hits_cols)              # Make a gRanges object from the hits file (with ID column indicating guides in the library to be reannotated)
      
      exon_hits <- suppressWarnings(findOverlaps(hits, exons))                  # Make a gRanges Hits object indicating the pairs of gRanges that overlap between query (hits) and subject (exome)
      exon_hitsRanges <- suppressWarnings(findOverlapPairs(hits, exons))        # Make a gRanges Pairs object indicating the genomic ranges of pairs of gRanges that overlap between the first (hits) and second (exome)
      sgrna_hits <- GRanges(exon_hitsRanges@first,                              # Make a gRanges object which contains the genomic ranges of the hits and the Symbols from the exome
                            exo_symbol = exons$exo_symbol[exon_hits@to],
                            exo_id = as.character(hits$exo_id[exon_hits@from]))
      mapping <- as_tibble(sgrna_hits) %>% unique()
      mapping$assembly <- blats$file_genome
      mapping$exome <- blats$file_exons
      
      distances <- suppressWarnings(distanceToNearest(hits, exons, select = "arbitrary")) # Make a gRanges Hits object showing the distances from each hit (queryHits) to the nearest exon in the exome (subjectHits)
      sgrna_distances <- GRanges(hits[distances@from],                          # Make a gRanges object annotated with the gene of the nearest exon and the distance to that exon
                                 nearestGene = exons$exo_symbol[distances@to],
                                 distance = distances@elementMetadata$distance,
                                 exo_id = as.character(hits$exo_id[distances@from]))
      
      final_distances <- as_tibble(sgrna_distances) %>% unique()
      final_distances$assembly <- blats$file_genome
      final_distances$exome <- blats$file_exons
      
      dir.create(dirname(blats$file_genomic_ranges_matched), recursive = T, showWarnings = F)
      dir.create(dirname(blats$file_genomic_ranges_distances), recursive = T, showWarnings = F)
      fwrite(mapping, blats$file_genomic_ranges_matched, sep = "\t")
      log_info("Wrote exonic hits to ", blats$file_genomic_ranges_matched)
      fwrite(final_distances, blats$file_genomic_ranges_distances, sep = "\t")
      log_info("Wrote exonic distances to ", blats$file_genomic_ranges_distances)
    }
  return()
}

# Write master mapping file
exorcisemaster <- function(premaster, blats, opt) {
  
  for (i in 1:length(premaster)) {
    premaster[[i]][which(is.na(premaster[[i]]) | premaster[[i]] == "")] <- "X"
  }
  
  if(is.null(opt$control)) {
    ncontrols <- length(premaster$exo_symbol[which(premaster$exo_symbol == "X")])
    premaster$exo_symbol[which(premaster$exo_symbol == "X")] <- paste0("Non-targeting", 1:ncontrols)
  }
  
  # Fix controls
  if (!is.null(opt$control)) {
    for (s in 1:length(opt$control)) {
      log_info("Replacing ", opt$control[s], " with ", opt$control_type[s])
      control_guides <- rep(F, nrow(premaster))
      if ("exo_orig" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control[s], premaster$exo_orig)
      if ("exo_id" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control[s], premaster$exo_id) # search for the control string in authors' symbols and IDs if those columns exist

        nThisControl <- length(premaster$exo_symbol[control_guides & premaster$exo_symbol == "X"])
        premaster$exo_symbol[control_guides & premaster$exo_symbol == "X"] <- paste0(opt$control_type[s], 1:nThisControl)   # Map control guides to control annotations unless there is an approved Symbol column
      }
    }
  
  master <- premaster %>% unique()
  
  # Output
  outfile_m_tsv <- blats$file_exorcise_master_out
  
  if ("exo_orig" %in% names(master)) {
    simple <- exorciseinferTargets(master, blats, opt)
    master <- left_join(master, simple, by = "exo_orig") %>%
      mutate(exo_harm = case_when(is.na(exo_harm) ~ exo_symbol,
                                  T ~ exo_harm)) %>%
      relocate(exo_id, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut)
  } else {
    master <- master %>% relocate(exo_id, exo_seq, exo_symbol, exo_target, exo_cut)
  }
  
  dirwrite(master, outfile_m_tsv, sep = "\t")
  log_info("Wrote master library to ", outfile_m_tsv)
  return(master)
}

# Infer intended target per authors' original symbols and write
exorciseinferTargets <- function(master, blats, opt) {
  
      simple <- master %>%
        dplyr::select(exo_symbol, exo_orig) %>% unique()

    # Load in gene symbol annotations for ranking. This file contains all possible approved symbols and gene classes
    annot <- fread(opt$priorities)
    names(annot) <- c("Approved_symbol", "Locus_group")
    annot$Locus_group <- factor(annot$Locus_group, ordered = T, levels = c("protein-coding gene", "non-coding RNA", "pseudogene", "other")) 
    
    # Infer
    log_info("Inferring intended targets per original symbol...")
    
    simple_dup <- simple %>%
      left_join(annot, by = c("exo_symbol" = "Approved_symbol"))                  # Join the multi-mapped original symbols to annotations for their mapping target
    simple_dup$Locus_group[which(is.na(simple_dup$Locus_group))] <- "other"
    
    simple_dup <- simple_dup %>% filter(!grepl(paste0(paste0("^", c("X", opt$control_type), "\\d+$"), collapse = "|"), exo_symbol)) # Remove controls and unmapped guides when considering intended target
    
    dup <- unique(simple_dup$exo_orig)
    
    # Fix multi-mapped original symbols
    j <- 1
    tot <- length(dup) 
    log_info("Inferring intended targets... ")
    
    simple_fixed <- data.table(exo_harm = rep("", tot),
                               exo_orig = rep("", tot),
                               Locus_group = rep("", tot),
                               consensus = rep("", tot),
                               fit_rank = rep("", tot))        # Initialise a tibble to contain fixed, singly-mapped original symbols
    
    for (s in dup) {  
      q <- simple_dup %>% filter(exo_orig == s)
      cons <- q %>%                                                                         # Determine the candidates that appeared the most often
        group_by(exo_symbol) %>%
        summarise(consensus = n()) %>%
        mutate(consensus = case_when(consensus == max(consensus) ~ 1,
                                     T ~ 0))
      q <- q %>% left_join(cons, by = "exo_symbol") %>% filter(consensus == 1)                  # Eliminate candidates which did not appear the most times
      q <- q %>% mutate(fit_rank = case_when(exo_symbol == exo_orig ~ 0,                 # Rank remaining, tied 1st place candidates: identical names > protein-coding > ncRNA > pseudogene > others
                                             T ~ as.numeric(Locus_group)))
      accept <- q %>% filter(fit_rank == min(fit_rank)) %>% arrange(exo_symbol) %>% head(1)     # If the top rank is still tied, accept the first match alphabetically
      set(simple_fixed, as.integer(j), "exo_harm", accept$exo_symbol[1])
      set(simple_fixed, as.integer(j), "exo_orig", accept$exo_orig[1])
      set(simple_fixed, as.integer(j), "Locus_group", accept$Locus_group[1])
      set(simple_fixed, as.integer(j), "consensus", accept$consensus[1])
      set(simple_fixed, as.integer(j), "fit_rank", accept$fit_rank[1])
      j <- j + 1
      if(j %% 1000 == 0) {
        log_info("Inferring intended targets... ", j, " out of ", tot, " done...")
      }
    }
    log_info("Inferring intended targets... ", tot, " out of ", tot, " done.")
    
    # Combine
      simple <- simple_fixed %>% dplyr::select(-Locus_group, -fit_rank, -consensus) %>% unique()
      
    return (simple)

}

reannotateExisting <- function(opt) {
  input <- fread(opt$infile) %>%
    relocate(exo_seq = opt$seq) %>%
    unique()
  library <- fread(opt$library)
  if("exo_harm" %in% names(library)) {
    library <- library %>% select(exo_id, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut) %>% unique()
  } else {
    library <- library %>% select(exo_id, exo_seq, exo_symbol, exo_target, exo_cut) %>% unique()
  }
  exorcised <- left_join(input, library, by = "exo_seq")
  if("exo_harm" %in% names(exorcised)) {
    exorcised <- exorcised %>% relocate(exo_id, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut)
  } else {
    exorcised <- exorcised %>% relocate(exo_id, exo_seq, exo_symbol, exo_target, exo_cut)
  }
  fwrite(exorcised, paste0(opt$outdir, "/exorcise.tsv"), sep = "\t")
}

fixOpts <- function(opt) {
  
  # fix comma-separated argument inputs
  if(!is.null(opt$infiles))      opt$infiles      <- commasplit(opt$infiles)
  if(!is.null(opt$outdir))       opt$outdir       <- commasplit(opt$outdir)
  if(!is.null(opt$seq))          opt$seq          <- commasplit(opt$seq)
  if(!is.null(opt$pam))          opt$pam          <- commasplit(opt$pam)
  #if(!is.null(opt$mode))         opt$mode         <- commasplit(opt$mode)
  if(!is.null(opt$library))      opt$library      <- commasplit(opt$library)
  if(!is.null(opt$genome))       opt$genome       <- commasplit(opt$genome)
  if(!is.null(opt$exome))        opt$exome        <- commasplit(opt$exome)
  if(!is.null(opt$priorities))   opt$priorities   <- commasplit(opt$priorities)
  if(!is.null(opt$id))           opt$id           <- commasplit(opt$id)
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
      } else {
        opt$infile <- opt$infile_glob
        warnings <- c(warnings, paste0("Warning: --infile accepted as globbed argument ", opt$infile))
      }
    }
  } else {
    errors <- c(errors, paste0("Error: --infile not specified."))
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
    }
  } else {
    errors <- c(errors, paste0("Error: --seq not specified."))
  }
  
  # check pam
  if(length(opt$pam) > 0) {
    if(length(opt$pam) > 1) {
      opt$pam <- opt$pam[1]
      warnings <- c(warnings, paste0("Warning: --pam received more than one argument. Using first value only: ", opt$pam, "."))
    }
    if(grepl("[^ACTGNactgn]", opt$pam)) {
      errors <- c(errors, paste0("Error: --pam ", opt$pam," contains non-nucleotide letters. Only [ACTGN] are accepted."))
    }
  } else {
    warnings <- c(warnings, paste0("Warning: --pam not specified. Not appending any PAM."))
  }
  
  # check mode
  if(length(opt$mode) > 0) {
    if(length(opt$mode) > 1) {
      opt$mode <- opt$mode[1]
      warnings <- c(warnings, paste0("Warning: --mode received more than one argument. Using first value only: ", opt$mode, "."))
    }
    opt$mode <- toupper(opt$mode)
    if(!(opt$mode %in% c("KO", "A", "I"))) {
      warnings <- c(warnings, paste0("Warning: --mode received an invalid value: ", opt$mode, ". Falling back to KO."))
      opt$mode <- "KO"
    }
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
      if(!("Approved_symbol" %in% opt$priorities_headers)) {
        errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " doesn't look like a feature priorities file. Expected an `Approved_symbol` column."))
      }
      if(!("Locus_group" %in% opt$priorities_headers)) {
        errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " doesn't look like a feature priorities file. Expected a `Locus_group` column."))
      }
    }
  } else if(opt$adhoc) {
    errors <- c(errors, paste0("Error: --priorities not passed while in ad-hoc mode."))
  }
  
  # check id
  if(length(opt$id) > 0) {
    if(length(opt$id) > 1) {
      opt$id <- opt$id[1]
      warnings <- c(warnings, paste0("Warning: --id received more than one argument. Using first value only: ", opt$id, "."))
    }
    if(!grepl("^\\d+$", opt$id)) {
      errors <- c(errors, paste0("Error: --id must be an integer, got ", opt$id, "."))
    } else {
      opt$id <- as.numeric(opt$id)
    }
  }
  
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
      } else if(length(control_type) != length(control)) {
        errors <- c(errors, paste0("Error: --control and --control_type lengths are not equal."))
      }
    }
  } else if(length(opt$control) > 0) {
    opt$control_type <- "Non-targeting"
    warnings <- c(warnings, paste0("Warning: --control passed without --control_type. Assuming --control_type is ", opt$control_type, "."))
  }
  
  # remove temporary options
  opt$infile_glob <- NULL
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
    stop("There were errors in the input. Please check the logfile: ", logfile)
  }
  
  return(opt)
}

#### MAIN ####

# Test vector
if (interactive()) {
  opt <- list()
  opt$infile <- "/Users/lam02/Library/CloudStorage/OneDrive-UniversityofCambridge/Projects/operation/crave/runs/Almu/NVS097/orig/cts/NVS097.counts.tsv.gz"
  opt$outdir <- "exorcise"
  opt$seq <- "1"
  opt$pam <- NULL
  opt$genome <- NULL
  opt$exome <- NULL
  opt$priorities <- NULL
  opt$id <- NULL
  opt$harm <- NULL
  opt$control <- NULL
  opt$control_type <- NULL
  opt$library <- "/Users/lam02/Downloads/exorcise/exorcisedsada.tsv"
}

## Execution
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
    #             help = "CRISPR screen type: KO (knockout), a (activation), i (inhibition)", metavar = "character"),
    make_option(opt_str = c("-l", "--library"), type = "character", default = NULL,
                help = "(required if --genome, --exome, and --priorities not specified) Exorcised file to be used as library.", metavar = "character"),
    make_option(opt_str = c("-v", "--genome"), type = "character", default = NULL,
                help = "(required if --library not specified) 2bit genome.", metavar = "character"),
    make_option(opt_str = c("-w", "--exome"), type = "character", default = NULL,
                help = "(required if --library not specified) Exome.", metavar = "character"),
    make_option(opt_str = c("-y", "--priorities"), type = "character", default = NULL,
                help = "(required if --library not specified) Priorities file.", metavar = "character"),
    make_option(opt_str = c("-j", "--id"), type = "character", default = NULL,
                help = "(optional, ignored if --library specified) ID column number.", metavar = "character"),
    make_option(opt_str = c("-n", "--harm"), type = "character", default = NULL,
                help = "(optional, ignored if --library specified) Existing annotation column number.", metavar = "character"),
    make_option(opt_str = c("-c", "--control"), type = "character", default = NULL,
                help = "(optional, ignored if --library specified) Pattern indicating a control guide (comma-separated list).", metavar = "character"),
    make_option(opt_str = c("-d", "--control_type"), type = "character", default = NULL,
                help = "(optional, ignored if --library specified) List of control guide types. Must be the same length as --control_strings (comma-separated list).", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
  if(!is.null(opt$project_name)) {
    logfile <- paste0(opt$outdir, "/", opt$project_name, "/logfile_exorcise_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
    dir.create(dirname(logfile), recursive = T, showWarnings = F)
    log_appender(appender_tee(logfile))
  }
  
  
  
  
}

# Parse arguments
opt <- fixOpts(opt)

start_time <- proc.time()
log_info("Welcome to exorcise, version ", ver, ".")
log_info("cwd: ", getwd())
command <- foreach(o = opt, .final = function(x) setNames(x, names(opt))) %do% { o }
command <- paste0("--", names(opt), " ", command)
command <- paste0(command, collapse = " ")
log_info("Command: exorcise.R ", command)

# Call main loop
reannotateLib(opt)
end_time <- proc.time()
log_info("exorcise process completed in ", (end_time - start_time)[[3]], " seconds.")

