#!/usr/bin/env Rscript
#
# author: Dr Simon Lam
# contact: sl681@cam.ac.uk
#
# Synopsis:
# Annotate guide RNA sequences by aligning to the exome.
#
# Description:
# This program accepts guide RNA sequences and returns all exonic matches to the
# desired exome. It accepts pre-existing guide annotations such as gene symbols,
# control guide designation, and arbitrary comments to facilitate your comparison
# with the output annotation. It reports on guides which hit exons, distances
# between missed-target guides and the nearest exon, and details of guides which
# hit multiple exons.
#
# This program expects CRISPRko guide sequences. Missed-target effects are
# guides that map to the genome but do not overlap exons by even one nucleotide;
# or guides that do not map to the genome at all.
# 

#### VERSION HISTORY ####
# version       datestamp             description
# 0.1           2022-12-09T12-28-16   initial release
# 0.11          2022-12-16T11-47-33   bugfixes and introduced new bugs
# 0.12          2023-02-01T11:52:40   enabled globbed reading from files (useful for reading from gzip)
# 0.13          2023-02-03T15:30:10   fixed erroneous candidate loss in the target inferring algorithm where candidate genes containing "X" were discarded
# 0.2           2023-02-20T17:00:30   fixed various bugs
# 0.21          2023-03-02T11:38:40   fixed issue with infile type inference
# 0.3           2023-03-10T13:54:14   refined to account for spCas9 cut site
# 0.31          2023-04-11T14:54:52   fixed to add genomic coordinates of guide hits outside of exomes
# 0.32          2023-05-25T10:55:03   enabled bypassing biomart mapping by setting the species to anything other than Human or Mouse
# 0.33          2023-05-25T14:17:41   increased stringency of BLAT. in a future version, these parameters will be tuneable

ver <- 0.33

#### INIT ####
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(readxl)
  library(writexl)
  library(tidyr)
  library(foreach)
  library(logger)
  library(R.utils)
  library(GenomicRanges)
  library(biomaRt)
})


#### FUNCTIONS ####

# Infer infile type and read it
mapper_read_infile <- function(infile, kind = NULL, skip = 0, header, sheet, nrow = Inf) {
  if(is.null(kind)) {
    if(grepl("\\.xlsx*$", infile)) kind <- "xls"
    else if(grepl("\\.tsv$", infile)) kind <- "tsv"
    else if(grepl("\\.csv$", infile)) kind <- "csv"
    else kind <- "txt"
  }
  if(grepl("xls", kind)) {
    return (read_excel(infile, skip = skip, sheet, col_names = header, n_max = nrow))
  } else if (grepl("tsv", kind)) {
    return (fread(infile, skip = skip, fill = T, header = header, nrow = nrow, sep = "\t"))
  } else if (grepl("csv", kind)) {
    return (fread(infile, skip = skip, fill = T, header = header, nrow = nrow, sep = ","))
  } else if (grepl("txt", kind)) {
    return (fread(infile, skip = skip, header = header, nrow = nrow))
  }
}

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
colonsplit <- function(x) unlist(strsplit(x, ":"))
rangesplit <- function(x) {
  x <- commasplit(x)
  x <- foreach(y = 1:length(x)) %do%
    if (grepl(":", x[y])) {
      colonsplit(x[y])[1]:colonsplit(x[y])[2]
    } else {
      as.numeric(x[y])
    }
  x
}

# Return true if file doesn't exist or has zero size
file.not.exist.or.zero <- function(x) !file.exists(Sys.glob(paste0(x, "*"))[1]) | file.exists(Sys.glob(paste0(x, "*"))[1]) & file.size(Sys.glob(paste0(x, "*"))[1]) == 0

# fread with hanging glob
fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

#### MAIN FUNCTION ####

reannotateLib <- function(opt) {
  authors <- importAuthorsLib(opt)
  
  blats <- list(file_genome = opt$file_genome,
                file_exons = opt$file_exons,
                file_feature_priorities = opt$file_feature_priorities,
                file_sgRNAs = paste0(opt$outdir, "/", opt$project_name, "/1-", opt$project_name, "_sgRNAs.fa"),
                file_psl = paste0(opt$outdir, "/", opt$project_name, "/2-", opt$project_name, "_",  sub(".+/(.+?)$", "\\1", opt$file_genome), "_BLAT.psl"),
                file_genomic_ranges = paste0(opt$outdir, "/", opt$project_name, "/3-", opt$project_name, "_", sub(".+/(.+?)$", "\\1", opt$file_genome), "_gR.tsv"),
                file_genomic_ranges_matched = paste0(opt$outdir, "/", opt$project_name, "/4a-", opt$project_name, "_", sub(".+/(.+?)$", "\\1", opt$file_exons), "_exon-hits.tsv"),
                file_genomic_ranges_distances = paste0(opt$outdir, "/", opt$project_name, "/4b-", opt$project_name, "_", sub(".+/(.+?)$", "\\1", opt$file_exons), "_exon-dist.tsv"),
                file_biomart_out = paste0(opt$outdir, "/", opt$project_name, "/5-", opt$project_name, "_biomart_out.tsv"),
                file_crispier_negctrl_out = paste0(opt$outdir, "/", opt$project_name, "/6c-", opt$project_name, "_negctrl.tsv"),
                file_crispier_master_out = paste0(opt$outdir, "/", opt$project_name, "/6a-", opt$project_name, "_master.tsv"),
                file_crispier_intended_targets_out = paste0(opt$outdir, "/", opt$project_name, "/6b-", opt$project_name, "_inferred.tsv"))
  
  extractGuides(opt, authors, blats)
  runBlat(blats)
  runPtgr(blats)
  runGrtem(blats)
  
  all_mappings <- foreach(b = 1:length(blats$file_genomic_ranges), .combine = "bind_rows") %do% { # generate mapping
    genome_hits <- fread(blats$file_genomic_ranges[b], colClasses = "character")
    exome_hits <- fread(blats$file_genomic_ranges_matched[b], colClasses = "character")
    left_join(genome_hits, exome_hits, by = c("seqnames", "start", "end", "ID", "strand", "assembly"))
  } %>%
  #filter(!is.na(Symbol)) %>%
  unique()
  
  genes <- all_mappings %>% dplyr::select(Symbol) %>% filter(Symbol != "") %>% unique() %>% unlist()
  premaster <- left_join(authors, all_mappings, by = "ID")
  
  if(opt$species == "Human") {
    biomart <- queryBiomaRtForHuman(genes, blats)
    premaster <- premaster %>% left_join(biomart, by = "Symbol") %>% unique()
  }
  if(opt$species == "Mouse") {
    biomart <- queryBiomaRtForMouse(genes, blats)
    premaster <- premaster %>% left_join(biomart, by = "Symbol") %>% unique()
  }
  
  master <- crispieRmaster(premaster, blats, opt)
  crispieRnegctrls(premaster, blats, opt)
  if(!is.null(opt$symbol_col)) { inferredTargets <- crispieRinferTargets(master, blats, opt) }
  
}

### OTHER FUNCTIONS ####

# Import authors' original library files, extracting original symbols, detecting/creating a primary key, and removing adapters where necessary
importAuthorsLib <- function(opt) {
  authors <- foreach(i = 1:length(opt$infiles), .combine = "rbind", .final = function(x) unique(x)) %do% {
    log_info("Importing file ", i, " of ", length(opt$infiles), ", ", opt$infiles[[i]])
    data <- mapper_read_infile(opt$infiles[[i]], opt$kind[[i]], opt$skip[[i]], opt$header[[i]], opt$sheet[[i]]) %>%                               # read authors' file
      dplyr::select(ID = opt$id_col[[i]], sgRNA = opt$grna_col[[i]], Original_symbol = opt$symbol_col[[i]], Notes = opt$comment_col[[i]]) %>%     # select appropriate columns
      unique()
    
    if(!is.null(opt$symbol_trim)) { # use regex to obtain symbols if specified to do so
      opt$symbol_trim[i] <- sub("^[\'\"]|[\'\"]$", "", opt$symbol_trim[i])
      log_info("Obtaining symbols using regex rule: ", opt$symbol_trim[i])
      data$Original_symbol <- gsub(opt$symbol_trim[i], "\\1", data$Original_symbol)
    }
    
    if(is.null(opt$id_col)) {                                                               # add id column if missing
      log_info("Adding sgRNA ID column")
      data <- cbind(tibble(ID = paste0(opt$project_name, "_", data$sgRNA)), data)
    } else if (any(duplicated(data$ID))) {
        log_info("Ignoring supplied sgRNA ID column due to duplicated IDs")
        log_info("Adding sgRNA ID column")
        data <- cbind(tibble(ID = paste0(opt$project_name, "_", data$sgRNA)), data %>% dplyr::select(-ID))
    } else {
      log_info("Accepting supplied sgRNA ID column")
    }
    
    if(sum(grepl("Original_symbol", names(data))) > 1) {
      log_info("Multiple original symbol columns detected. Melting")
      data <- data %>% pivot_longer(cols = grep("Original_symbol", names(data))) %>% dplyr::select(-name)           # melt additional symbol columns
      data <- data %>% filter(!is.na(value))                                                                        # assume all NAs arose by melting a non-rectangular symbol matrix and remove them
    }
    
    if(!is.null(opt$adapters)) { # remove adapters if specified
      for (adapter in opt$adapters) {
        log_info("Removing adapter from guide sequences: ", adapter)
        data$sgRNA <- gsub(paste0("^", adapter), "", data$sgRNA)
        data$sgRNA <- gsub(paste0(adapter, "$"), "", data$sgRNA)
      }
    }
    
    if(!is.null(opt$symbol_col)) {
      names(data)[1:3] <- c("ID", "sgRNA", "Original_symbol") # be sure that column names are correct
    } else {
      names(data)[1:2] <- c("ID", "sgRNA") # be sure that column names are correct
    }
    
    if(!is.null(opt$symbol_col)) {
      data <- as_tibble(lapply(data, as.character)) %>%
        filter(!is.na(Original_symbol)) %>%
        filter(Original_symbol != "") # remove NA author's symbols
    }
    
    if(nrow(data) > 0) { data } else { NULL }
  }
  return(authors)
}

# Use authors' original library to extract guide RNA sequences for alignment
extractGuides <- function(opt, authors, blats) {
  log_info("Creating sgRNA FASTA")
  
  if(is.null(opt$PAM)) opt$PAM <- ""
  
  sgRNAs <- foreach (PAM = opt$PAM, .combine = "c") %do% {
    paste0(">", authors$ID, "\n", authors$sgRNA, PAM) # convert to fasta
  }
  
  dirwrite(list(sgRNAs), blats$file_sgRNAs, sep = "\n", quote = F)
  log_info("Wrote sgRNA FASTA to ", blats$file_sgRNAs)
  return()
}

# Run BLAT
runBlat <- function(blats) {
  blat_command <- "blat"                                                   # external scripts and common parameters
  blat_params <- "-stepSize=4 -tileSize=10 -fine -repMatch=2000000 -minScore=20 -minIdentity=100"
  
  for (b in 1:length(blats$file_genome)) {
    while(file.not.exist.or.zero(blats$file_psl[b])) {
      dir.create(dirname(blats$file_psl[b]), showWarnings = F)
      run_command <- paste(blat_command, blats$file_genome[b], blats$file_sgRNAs, blats$file_psl[b], blat_params)
      log_info("Sending to BLAT: ", run_command)
      system(run_command)
    }
  }
  
  return()
}

# Convert psl to genomic ranges
runPtgr <- function(blats) {
  
  for (b in 1:length(blats$file_psl)) {
    while(file.not.exist.or.zero(blats$file_genomic_ranges[b])) {
      infile <- blats$file_psl[b]
      outfile <- blats$file_genomic_ranges[b]
      
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
                        ID = psl$`Q name`,
                        strand = psl$strand)
      out <- tibble(seqnames = as.character(rep(ranges@seqnames@values, ranges@seqnames@lengths)),                                                # evaluate RLE to get the seqnames
                    guideBegin = ranges@ranges@start,
                    guideFinal = ranges@ranges@start + ranges@ranges@width - 1,                                                                          # end is start + width - 1
                    ID = ranges$ID,
                    strand = rep(ranges@strand@values, ranges@strand@lengths),
                    assembly = blats$file_genome[b])
      out <- out %>% transmute(seqnames, guideBegin, guideFinal, ID, strand, assembly,
                               start = case_when(strand == "+" ~ guideFinal - 4, # -3 to -4 upstream from PAM is the cut site
                                                 strand == "-" ~ guideBegin + 3),
                               end = case_when(strand == "+" ~ guideFinal - 3,
                                               strand == "-" ~ guideBegin + 4))
      dir.create(dirname(outfile), recursive = T, showWarnings = F)							
      fwrite(out, outfile, sep = "\t") 
      log_info("Wrote genomic hits file to ", outfile)
    }
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
  exome_cols$symbol <- grep("symbol|name2", exome_headers)
  return(exome_cols)
}

inferHitsCols <- function(file_genomic_ranges) {
  hits_headers <- fread(file_genomic_ranges, nrows = 0) %>% names() %>% tolower()
  hits_cols <- list()
  hits_cols$chr <- grep("chr|seqname", hits_headers)
  hits_cols$start <- grep("start", hits_headers)
  hits_cols$end <- grep("end", hits_headers)
  hits_cols$strand <- grep("str", hits_headers)
  hits_cols$id <- grep("^id$", hits_headers)
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
  ranges <- IRanges(start = exome[[as.numeric(exome_cols$start)]],
                    end = exome[[as.numeric(exome_cols$end)]])
  if (exome_cols$strand == "*") {                    # use * if in strandless mode
    strand <- rep("*", nrow(exome))
  } else {
    strand <- exome[[as.numeric(exome_cols$strand)]]
  }
  symbol <- exome[[as.numeric(exome_cols$symbol)]] %>% gsub(pattern = ",$", replacement = "")
  
  exome <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand, Symbol = symbol)
  return(exome)
}

import_hits <- function(file_genomic_ranges, hits_cols) {
  
  hits <- fread(file_genomic_ranges)
  seqnames <- hits[[as.numeric(hits_cols$chr)]]
  ranges <- IRanges(start = hits[[as.numeric(hits_cols$start)]],
                    end = hits[[as.numeric(hits_cols$end)]])
  ID = hits[[as.numeric(hits_cols$id)]]
  if (hits_cols$strand == "*") {                    # use * if in strandless mode
    strand <- rep("*", nrow(hits))
  } else {
    strand <- hits[[as.numeric(hits_cols$strand)]]
  }
  
  hits <- GRanges(seqnames = seqnames, ranges = ranges, strand = strand, ID = as.character(ID))
  return(hits)
}

runGrtem <- function(blats) {
  
  for (b in 1:length(blats$file_genomic_ranges)) {
    while(file.not.exist.or.zero(blats$file_genomic_ranges_matched[b]) | file.not.exist.or.zero(blats$file_genomic_ranges_distances[b])) {
      log_info("Obtaining exonic hits... ", blats$file_exons[b], "...")
      
      exome_cols <- inferExomeCols(blats$file_exons[b])                         # Infer column identities in the exome file
      hits_cols <- inferHitsCols(blats$file_genomic_ranges[b])                  # Infer column identities in the hits file
      
      exons <- import_exome(blats$file_exons[b], exome_cols)                    # Make a gRanges object from the exome file (with Symbol column to provide re-annotations)
      hits <- import_hits(blats$file_genomic_ranges[b], hits_cols)              # Make a gRanges object from the hits file (with ID column indicating guides in the library to be reannotated)
      
      exon_hits <- suppressWarnings(findOverlaps(hits, exons))                  # Make a gRanges Hits object indicating the pairs of gRanges that overlap between query (hits) and subject (exome)
      exon_hitsRanges <- suppressWarnings(findOverlapPairs(hits, exons))        # Make a gRanges Pairs object indicating the genomic ranges of pairs of gRanges that overlap between the first (hits) and second (exome)
      sgrna_hits <- GRanges(exon_hitsRanges@first,                              # Make a gRanges object which contains the genomic ranges of the hits and the Symbols from the exome
                            Symbol = exons$Symbol[exon_hits@to],
                            ID = as.character(hits$ID[exon_hits@from]))
      mapping <- as_tibble(sgrna_hits) %>% unique()
      mapping$assembly <- blats$file_genome[b]
      mapping$exome <- blats$file_exons[b]
      
      distances <- suppressWarnings(distanceToNearest(hits, exons, select = "arbitrary")) # Make a gRanges Hits object showing the distances from each hit (queryHits) to the nearest exon in the exome (subjectHits)
      sgrna_distances <- GRanges(hits[distances@from],                          # Make a gRanges object annotated with the gene of the nearest exon and the distance to that exon
                                 nearestGene = exons$Symbol[distances@to],
                                 distance = distances@elementMetadata$distance,
                                 ID = as.character(hits$ID[distances@from]))
      
      final_distances <- as_tibble(sgrna_distances) %>% unique()
      final_distances$assembly <- blats$file_genome[b]
      final_distances$exome <- blats$file_exons[b]
      
      dir.create(dirname(blats$file_genomic_ranges_matched[b]), recursive = T, showWarnings = F)
      dir.create(dirname(blats$file_genomic_ranges_distances[b]), recursive = T, showWarnings = F)
      fwrite(mapping, blats$file_genomic_ranges_matched[b], sep = "\t")
      log_info("Wrote exonic hits to ", blats$file_genomic_ranges_matched[b])
      fwrite(final_distances, blats$file_genomic_ranges_distances[b], sep = "\t")
      log_info("Wrote exonic distances to ", blats$file_genomic_ranges_distances[b])
    }
  }
  return()
}

# Query biomaRt for human
queryBiomaRtForHuman <- function(genes, blats) {
  while(file.not.exist.or.zero(blats$file_biomart_out)) {
    try({
      
      mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl")
    
      res1 <- getBM(filters = c("entrezgene_accession", "chromosome_name"),
                    attributes = c("entrezgene_accession", "entrezgene_id", "hgnc_id", "ensembl_gene_id"),
                    values = list(entrezgene_accession = genes, chromosome_name = c(1:22, "X", "Y")), mart = mart)
      res2 <- getBM(filters = c("external_gene_name", "chromosome_name"),
                    attributes = c("external_gene_name", "entrezgene_id", "hgnc_id", "ensembl_gene_id"),
                    values = list(external_gene_name = genes, chromosome_name = c(1:22, "X", "Y")), mart = mart)
      res3 <- getBM(filters = c("hgnc_symbol", "chromosome_name"),
                    attributes = c("hgnc_symbol", "entrezgene_id", "hgnc_id", "ensembl_gene_id"),
                    values = list(hgnc_symbol = genes, chromosome_name = c(1:22, "X", "Y")), mart = mart)
      
      res <- tibble(Symbol = genes) %>%
        left_join(res1, by = c("Symbol" = "entrezgene_accession")) %>%
        left_join(res2, by = c("Symbol" = "external_gene_name")) %>%
        left_join(res3, by = c("Symbol" = "hgnc_symbol")) %>%
        transmute(Symbol,
                  Ensembl_gene_ID = case_when(!is.na(ensembl_gene_id.x) ~ ensembl_gene_id.x,
                                              !is.na(ensembl_gene_id.y) ~ ensembl_gene_id.y,
                                              !is.na(ensembl_gene_id)   ~ ensembl_gene_id),
                  NCBI_gene_ID = case_when(!is.na(entrezgene_id.x) ~ entrezgene_id.x,
                                           !is.na(entrezgene_id.y) ~ entrezgene_id.y,
                                           !is.na(entrezgene_id)   ~ entrezgene_id),
                  HGNC_ID = case_when(!is.na(hgnc_id.x) ~ hgnc_id.x,
                                      !is.na(hgnc_id.y) ~ hgnc_id.y,
                                      !is.na(hgnc_id)   ~ hgnc_id)) %>%
        unique()
      
      dir.create(dirname(blats$file_biomart_out), recursive = T, showWarnings = F)
      fwrite(res, blats$file_biomart_out, sep = "\t")
    }, silent = T)
  }
  
  biomart <- fread(blats$file_biomart_out)
  return(biomart)
}

# Query biomaRt for mouse
queryBiomaRtForMouse <- function(genes, blats) {
  while(file.not.exist.or.zero(blats$file_biomart_out)) {
    mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                       dataset = "mmusculus_gene_ensembl")
    
    res1 <- getBM(filters = c("mgi_symbol", "chromosome_name"),
                  attributes = c("mgi_symbol", "entrezgene_id", "ensembl_gene_id"),
                  values = list(mgi_symbol = genes, chromosome_name = c(1:19, "X", "Y")), mart = mart)
    res2 <- getBM(filters = c("ensembl_gene_id", "chromosome_name"),
                  attributes = c("ensembl_gene_id" , "hsapiens_homolog_ensembl_gene"),
                  values = list(ensembl_gene_id = res1$ensembl_gene_id, chromosome_name = c(1:19, "X", "Y")), mart = mart)
    
    hmart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL",
                        dataset = "hsapiens_gene_ensembl")
    
    res3 <- getBM(filters = c("ensembl_gene_id", "chromosome_name"),
                  attributes = c("hgnc_symbol", "entrezgene_id", "hgnc_id", "ensembl_gene_id"),
                  values = list(ensembl_gene_id = res2$hsapiens_homolog_ensembl_gene, chromosome_name = c(1:22, "X", "Y")), mart = hmart)
    
    res <- tibble(Symbol = genes) %>%
      left_join(res1, by = c("Symbol" = "mgi_symbol")) %>%
      left_join(res2, by = "ensembl_gene_id") %>%
      left_join(res3, by = c("hsapiens_homolog_ensembl_gene" = "ensembl_gene_id")) %>%
      transmute(Symbol,
                Ensembl_gene_ID = ensembl_gene_id,
                NCBI_gene_ID = as.character(entrezgene_id.x),
                Human_symbol = hgnc_symbol,
                Human_HGNC_ID = hgnc_id,
                Human_Ensembl_gene_ID = hsapiens_homolog_ensembl_gene,
                Human_NCBI_gene_ID = as.character(entrezgene_id.y)) %>% 
      unique()
    
    dir.create(dirname(blats$file_biomart_out), recursive = T, showWarnings = F)
    fwrite(res, blats$file_biomart_out, sep = "\t")
  }
  
  biomart <- fread(blats$file_biomart_out)
  return(biomart)
}

# Indicate negative controls
crispieRnegctrls <- function(premaster, blats, opt) {
  
  # Fix empties
  for (i in 1:length(premaster)) {
    premaster[[i]][which(is.na(premaster[[i]]) | premaster[[i]] == "")] <- "X"
  }
  
  # Create a report of negative controls which map to exons
  if (!is.null(opt$control_string)) {
    log_info("Creating a report of negative controls which map to exons")
    
    all_control_guides <- rep(F, nrow(premaster))
    
    neg_control_strings <- opt$control_string[opt$control_type == "Non-targeting"]
    
    for (s in 1:length(opt$control_string)) {
      log_info("Replacing ", opt$control_string[s], " with ", opt$control_type[s])
      control_guides <- rep(F, nrow(premaster))
      if ("Original_symbol" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control_string[s], premaster$Original_symbol)
      if ("ID" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control_string[s], premaster$ID) # search for the control string in authors' symbols and IDs if those columns exist
      if (!is.null(opt$comment_col)) {
        for (m in grep("^Notes", names(premaster))) {
          control_guides <- control_guides | grepl(opt$control_string[s], premaster[[m]]) # also search authors' notes if specified
        }      }
      
      all_control_guides <- all_control_guides | control_guides
      
    }
    
    negctrl_mappings <- premaster %>% filter(all_control_guides) # Output a report of whether negative controls had any hits
    negctrl_mappings <- negctrl_mappings %>% unique()
    
    outfile_negctrl_mappings_tsv <- blats$file_crispier_negctrl_out
    dirwrite(negctrl_mappings, outfile_negctrl_mappings_tsv, sep = "\t")
    log_info("Wrote negative control mapping report to ", outfile_negctrl_mappings_tsv)
    
  }
  
  return()
}

# Write master mapping file
crispieRmaster <- function(premaster, blats, opt) {
  
  for (i in 1:length(premaster)) {
    premaster[[i]][which(is.na(premaster[[i]]) | premaster[[i]] == "")] <- "X"
  }
  
  if(is.null(opt$control_string)) {
    premaster$Symbol[which(premaster$Symbol == "X")] <- "Non-targeting"
  }
  
  # Fix controls
  if (opt$control_retain) premaster <- premaster %>% mutate(Control_type = NA) # initialise new column to contain control types if specified
  
  if (!is.null(opt$control_string)) {
    for (s in 1:length(opt$control_string)) {
      log_info("Replacing ", opt$control_string[s], " with ", opt$control_type[s])
      control_guides <- rep(F, nrow(premaster))
      if ("Original_symbol" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control_string[s], premaster$Original_symbol)
      if ("ID" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control_string[s], premaster$ID) # search for the control string in authors' symbols and IDs if those columns exist
      if (!is.null(opt$comment_col)) {
        for (m in grep("^Notes", names(premaster))) {
          control_guides <- control_guides | grepl(opt$control_string[s], premaster[[m]]) # also search authors' notes if specified
        }
      }
      if (opt$control_retain) {
        premaster$Control_type[control_guides] <- opt$control_type[s]   # Map control guides to control annotations in a new column
      } else {
        premaster$Symbol[control_guides] <- opt$control_type[s]   # Map control guides to control annotations and overwrite approved Symbol column
        for (t in grep("Ensembl|NCBI|HGNC", names(premaster))) {
          premaster[[t]][control_guides] <- "X"                   # Remove external annotations
        }
      }
    }
  }
  
  master <- premaster %>% unique()
  
  # Output
  outfile_m_tsv <- blats$file_crispier_master_out
  dirwrite(master, outfile_m_tsv, sep = "\t")
  log_info("Wrote master library to ", outfile_m_tsv)
  return(master)
}

# Infer intended target per authors' original symbols and write
crispieRinferTargets <- function(master, blats, opt) {
  
  if(file.not.exist.or.zero(blats$file_crispier_intended_targets_out)) {
    
    if (opt$control_retain) {
      simple <- master %>%
        dplyr::select(grep("ymbol|Ensembl|NCBI|HGNC|Control_type", names(master)))
    } else {
      simple <- master %>%
        dplyr::select(grep("ymbol|Ensembl|NCBI|HGNC", names(master)))
    }
    
    # Load in gene symbol annotations for ranking. This file contains all possible approved symbols and gene classes
    annot <- fread(opt$file_feature_priorities)
    names(annot) <- c("Approved_symbol", "Locus_group")
    if (opt$species == "Human") {
      annot$Locus_group <- factor(annot$Locus_group, ordered = T, levels = c("protein-coding gene", "non-coding RNA", "pseudogene", "other")) 
    } else if (opt$species == "Mouse") {
      annot$Locus_group <- factor(annot$Locus_group, ordered = T, levels = c("protein coding gene", "pseudogene", "polymorphic pseudogene", "ribozyme gene", "rRNA gene", "scRNA gene", "snoRNA gene", "snRNA gene", "tRNA gene", "RNase MRP RNA gene", "RNase P RNA gene", "SRP RNA gene", "telomerase RNA gene", "miRNA gene", "antisense lncRNA gene", "sense intronic lncRNA gene", "sense overlapping lncRNA gene", "non-coding RNA gene", "bidirectional promoter lncRNA gene", "lncRNA gene", "lincRNA gene", "gene segment", "pseudogenic gene segment", "unclassified gene", "unclassified non-coding RNA gene", "heritable phenotypic marker", "other")) 
    }
    
    # Infer
    log_info("Inferring intended targets per original symbol...")
    
    simple_dup <- simple %>%
      left_join(annot, by = c("Symbol" = "Approved_symbol"))                  # Join the multi-mapped original symbols to annotations for their mapping target
    simple_dup$Locus_group[which(is.na(simple_dup$Locus_group))] <- "other"
    
    simple_dup <- simple_dup %>% filter(!grepl(paste0(c("^X$", opt$control_type), collapse = "|"), Symbol)) # Remove controls and unmapped guides when considering intended target
    
    dup <- unique(simple_dup$Original_symbol)
    
    # Fix multi-mapped original symbols
    simple_fixed <- tibble()                                                                # Initialise a tibble to contain fixed, singly-mapped original symbols
    j <- 1
    tot <- length(dup) 
    log_info("Inferring intended targets... ")
    for (s in dup) {  
      q <- simple_dup %>% filter(Original_symbol == s)
      cons <- q %>%                                                                         # Determine the candidates that appeared the most often
        group_by(Symbol) %>%
        summarise(consensus = n()) %>%
        mutate(consensus = case_when(consensus == max(consensus) ~ 1,
                                     T ~ 0))
      q <- q %>% left_join(cons, by = "Symbol") %>% filter(consensus == 1)                  # Eliminate candidates which did not appear the most times
      q <- q %>% mutate(fit_rank = case_when(Symbol == Original_symbol ~ 0,                 # Rank remaining, tied 1st place candidates: identical names > protein-coding > ncRNA > pseudogene > others
                                             T ~ as.numeric(Locus_group)))
      accept <- q %>% filter(fit_rank == min(fit_rank)) %>% arrange(Symbol) %>% head(1)     # If the top rank is still tied, accept the first match alphabetically
      simple_fixed <- rbind(simple_fixed, accept)
      j <- j + 1
      if(j %% 1000 == 0) {
        log_info("Inferring intended targets... ", j, " out of ", tot, " done...")
      }
    }
    log_info("Inferring intended targets... ", tot, " out of ", tot, " done.")
    
    # Combine
    if("Locus_group" %in% names(simple_fixed)) {
      simple <- simple_fixed %>% dplyr::select(-Locus_group, -fit_rank, -consensus)
    } else {
      simple <- simple_fixed
    }
    # }
    
    simple <- simple %>%
      dplyr::select(grep("ymbol|Ensembl|NCBI|HGNC", names(simple))) %>%
      unique()
    
    # Attach controls mappings without inference
    simple_controls <- master %>%
      dplyr::select(Original_symbol, Symbol) %>%
      unique() %>%
      filter(Symbol %in% opt$control_type) %>%
      unique()
    simple <- bind_rows(simple, simple_controls)
    
    outfile_s_tsv <- blats$file_crispier_intended_targets_out
    dirwrite(simple, outfile_s_tsv, sep = "\t")
    log_info("Wrote simplified library to ", outfile_s_tsv)
  }
  
}

fixOpts <- function(opt) {
  
  # fix comma-separated argument inputs
  if(!is.null(opt$infiles))        opt$infiles        <- commasplit(opt$infiles)
  if(!is.null(opt$kind))           opt$kind           <- commasplit(opt$kind)
  if(!is.null(opt$sheet))          opt$sheet          <- commasplit(opt$sheet) %>% as.numeric()
  if(!is.null(opt$id_col))         opt$id_col         <- commasplit(opt$id_col) %>% as.numeric()
  if(!is.null(opt$grna_col))       opt$grna_col       <- commasplit(opt$grna_col) %>% as.numeric()
  if(!is.null(opt$adapters))       opt$adapters       <- commasplit(opt$adapters)
  if(!is.null(opt$PAM))            opt$PAM            <- commasplit(opt$PAM)
  if(!is.null(opt$symbol_col))     opt$symbol_col     <- rangesplit(opt$symbol_col)
  if(!is.null(opt$comment_col))    opt$comment_col    <- rangesplit(opt$comment_col)
  if(!is.null(opt$skip)) {         opt$skip           <- commasplit(opt$skip) %>% as.numeric() } else opt$skip <- 0
  if(!is.null(opt$control_string)) opt$control_string <- commasplit(opt$control_string)
  if(!is.null(opt$control_type))   opt$control_type   <- commasplit(opt$control_type)
  if(!is.null(opt$symbol_trim))    opt$symbol_trim    <- commasplit(opt$symbol_trim)
  if(!is.null(opt$file_genome))    opt$file_genome    <- commasplit(opt$file_genome)
  if(!is.null(opt$file_exons))     opt$file_exons     <- commasplit(opt$file_exons)
  if(!is.null(opt$file_feature_priorities)) opt$file_feature_priorities <- commasplit(opt$file_feature_priorities)
  if(is.null(opt$control_retain))  opt$control_retain <- F
  
  # check arguments
  errors <- list()
  warnings <- list()
  # check if infiles exist
  for (x in opt$infiles) {
    if(!file.exists(Sys.glob(paste0(x, "*"))[1])) {
      errors <- c(errors, paste0("File not found: ", x, "."))
    }
  }
  
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
  
  warnings <- c(warnings, lengthChecker(opt$kind, opt$infiles, "kind"))
  warnings <- c(warnings, lengthChecker(opt$sheet, opt$infiles, "sheet"))
  warnings <- c(warnings, lengthChecker(opt$id_col, opt$infiles, "id_col"))
  warnings <- c(warnings, lengthChecker(opt$grna_col, opt$infiles, "grna_col"))
  warnings <- c(warnings, lengthChecker(opt$symbol_col, opt$infiles, "symbol_col"))
  warnings <- c(warnings, lengthChecker(opt$comment_col, opt$infiles, "comment_col"))
  warnings <- c(warnings, lengthChecker(opt$skip, opt$infiles, "skip"))
  warnings <- c(warnings, lengthChecker(opt$header, opt$infiles, "header"))
  warnings <- c(warnings, lengthChecker(opt$control_type, opt$control_string, "control_type"))
  
  # check if kind is valid
  for (x in opt$kind) {
    if(!(grepl("xls|sv|txt", x))) {
      errors <- c(errors, paste0("Invalid input to --kind: ", x, ". Expected xls, tsv, csv, or txt."))
    }
  }
  
  # check if mode is valid
  for (x in opt$mode) {
    if(!(grepl("KO$|a$|i$", x, ignore.case = T))) {
      errors <- c(errors, paste0("Invalid input to --mode: ", x, ". Expected KO, a, or i."))
    }
  }
  
  # check if sheet is valid
  for (x in opt$sheet) {
    if(!(grepl("^\\d+$", x))) {
      errors <- c(errors, paste0("Invalid input to --sheet: ", x, ". Must be a number."))
    }
  }
  
  # expand argument inputs of length one
  if(length(opt$kind)        != length(opt$infiles)) { opt$kind        <- rep(opt$kind[1], length(opt$infiles)) }
  if(length(opt$sheet)       != length(opt$infiles)) { opt$sheet       <- rep(opt$sheet[1], length(opt$infiles)) }
  if(length(opt$id_col)      != length(opt$infiles)) { opt$id_col      <- rep(opt$id_col[1], length(opt$infiles)) }
  if(length(opt$grna_col)    != length(opt$infiles)) { opt$grna_col    <- rep(opt$grna_col[1], length(opt$infiles)) }
  if(length(opt$symbol_col)  != length(opt$infiles)) { opt$symbol_col  <- rep(opt$symbol_col[1], length(opt$infiles)) }
  if(length(opt$comment_col) != length(opt$infiles)) { opt$comment_col <- rep(opt$comment_col[1], length(opt$infiles)) }
  if(length(opt$skip)        != length(opt$infiles)) { opt$skip        <- rep(opt$skip[1], length(opt$infiles)) }
  if(length(opt$header)      != length(opt$infiles)) { opt$header      <- rep(opt$header[1], length(opt$infiles)) }
  
  
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

## Execution (IDE)
if (interactive()) {
  opt <- list()
  opt$infiles <- "data/combo_library.1.3guides.tsv"
  opt$mode <- NULL
  opt$kind <- NULL
  opt$sheet <- NULL
  opt$id_col <- NULL
  opt$grna_col <- "1"
  opt$adapters <- NULL
  opt$PAM <- NULL
  opt$symbol_col <- NULL
  opt$symbol_trim <- NULL
  opt$comment_col <- "2"
  opt$skip <- NULL
  opt$header <- T
  opt$species <- "Human"
  opt$control_string <- "NonTargeting"
  opt$control_type <- "Non-targeting"
  opt$control_retain <- NULL
  opt$project_name <- "Panton"
  opt$file_genome <- "/Users/simonlam/bio/Projects/dev/crispieR/data/hg38.2020-09-22.2bit"
  opt$file_exons <- "/Users/simonlam/bio/Projects/dev/crispieR/data/hg38.refseq.exons.tsv"
  opt$file_feature_priorities <- "/Users/simonlam/bio/Projects/dev/crispieR/data/symbol_ids_table.csv"
  opt$outdir <- "crispieR/"
}

## Execution (commandline)
if (!interactive()) {
  
  option_list <- list(
    make_option(opt_str = c("-i", "--infiles"), type = "character", default = NULL,
                help = "Path to file(s) containing sgRNA sequences and authors' gene symbols.", metavar = "character"),
    make_option(opt_str = c("-q", "--mode"), type = "character", default = "KO",
                help = "CRISPR screen type: KO (knockout), a (activation), i (inhibition)", metavar = "character"),
    make_option(opt_str = c("-k", "--kind"), type = "character", default = NULL,
                help = "Are the infile(s) plaintext or Excel format?", metavar = "character"),
    make_option(opt_str = c("-e", "--sheet"), type = "character", default = NULL,
                help = "Numbers of sheets in xlsx infiles containing sgRNA sequences and authors' gene symbols (one value per infile, comma-separated).", metavar = "character"),
    make_option(opt_str = c("-j", "--id_col"), type = "character", default = NULL,
                help = "Column number(s) containing sgRNA IDs (one value per infile, comma-separated).", metavar = "character"),
    make_option(opt_str = c("-g", "--grna_col"), type = "character", default = NULL,
                help = "Column number(s) containing sgRNA sequences (one value per infile, comma-separated).", metavar = "character"),
    make_option(opt_str = c("-x", "--adapters"), type = "character", default = NULL,
                help = "Adapters to remove from 5' and/or 3' ends of oligo sequences to obtain the sgRNA (comma-separated).", metavar = "character"),
    make_option(opt_str = c("-z", "--PAM"), type = "character", default = NULL,
                help = "Protoadjacent motif to append to 3' ends of sgRNA sequences.", metavar = "character"),
    make_option(opt_str = c("-n", "--symbol_col"), type = "character", default = NULL,
                help = "Column number(s) containing authors' gene symbols (one value or range per infile, comma-separated)", metavar = "character"),
    make_option(opt_str = c("-t", "--symbol_trim"), type = "character", default = NULL,
                help = "Regex to obtain symbols from symbol column. Expects one match group.", metavar = "character"),
    make_option(opt_str = c("-m", "--comment_col"), type = "character", default = NULL,
                help = "Column number(s) containing authors' comments (one value per infile, comma-separated).", metavar = "character"),
    make_option(opt_str = c("-s", "--skip"), type = "character", default = NULL,
                help = "Number of rows to skip from the top of the file (one value per infile, comma-separated).", metavar = "character"),
    make_option(opt_str = c("-u", "--header"), type = "logical", default = T,
                help = "Does the infile have a header row?", metavar = "character"),
    make_option(opt_str = c("-a", "--species"), type = "character", default = NULL,
                help = "Was the study done on Human or Mouse?", metavar = "character"),
    make_option(opt_str = c("-c", "--control_string"), type = "character", default = NULL,
                help = "List of strings in the authors' gene symbols that correspond to control sgRNAs (comma-separated).", metavar = "character"),
    make_option(opt_str = c("-d", "--control_type"), type = "character", default = NULL,
                help = "List of control types corresponding to each string in --control_strings. Must be the same length as --control_strings (comma-separated).", metavar = "character"),
    make_option(opt_str = c("-b", "--control_retain"), type = "logical", default = F,
                help = "Retain mappings for control guides? Default: F.", metavar = "character"),
    make_option(opt_str = c("-p", "--project_name"), type = "character", default = NULL,
                help = "Name to append to output files.", metavar = "character"),
    make_option(opt_str = c("-v", "--file_genome"), type = "character", default = NULL,
                help = "Path to 2bit assembly file for BLAT. See example.", metavar = "character"),
    make_option(opt_str = c("-w", "--file_exons"), type = "character", default = NULL,
                help = "Path to exons file. See example.", metavar = "character"),
    make_option(opt_str = c("-y", "--file_feature_priorities"), type = "character", default = NULL,
                help = "Path to feature priorities file. See example.", metavar = "character"),
    make_option(opt_str = c("-o", "--outdir"), type = "character", default = NULL,
                help = "Path to place output files.", metavar = "character")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
  if(!is.null(opt$project_name)) {
    logfile <- paste0(opt$outdir, "/", opt$project_name, "/logfile_crispieR-lib_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
    dir.create(dirname(logfile), recursive = T, showWarnings = F)
    log_appender(appender_tee(logfile))
  }
  
  
  
  
}

# Parse arguments
opt <- fixOpts(opt)

start_time <- proc.time()
log_info("Welcome to crispieR-lib, version ", ver, ".")
log_info("cwd: ", getwd())
command <- foreach(o = opt, .final = function(x) setNames(x, names(opt))) %do% { o }
command <- paste0("--", names(opt), " ", command)
command <- paste0(command, collapse = " ")
log_info("Command: crispieR-lib.R ", command)

# Call main loop
reannotateLib(opt)
end_time <- proc.time()
log_info("crispieR-lib process completed in ", (end_time - start_time)[[3]], " seconds.")

