#### Map hits to exome ####

# Convert genomic ranges to exon hits
inferExomeCols <- function(file_exons) {
  exome_headers <- fread(file_exons, nrows = 0, skip = "Starts") %>% names()
  exome_cols <- list()
  exome_cols$chr <- grep("chr|seqname", exome_headers)
  exome_cols$start <- grep("exonStarts", exome_headers)
  exome_cols$end <- grep("exonEnds", exome_headers)
  exome_cols$strand <- "*"                                                      # uncomment to search both strands (default)
  # exome_cols$strand <- grep("str", exome_headers)                             # uncomment to search sense strand only
  exome_cols$symbol <- grep("symbol|name2", exome_headers)
  exome_cols$inherit <- grep("chr|seqname|exonStarts|exonEnds|str|symbol|name2|name$|cdsStart|cdsEnd|exonFrame|^ranges$|^seqlevels$|^seqlengths$|^isCircular$|^start$|^end$|^width$|^element$", exome_headers, value = T, invert = T)
  return(exome_cols)
}

inferHitsCols <- function(file_genomic_ranges) {
  hits_headers <- fread(file_genomic_ranges, nrows = 0) %>% names()
  hits_cols <- list()
  hits_cols$chr <- grep("chr|seqname", hits_headers)
  hits_cols$start <- grep("start", hits_headers)
  hits_cols$end <- grep("end", hits_headers)
  hits_cols$strand <- grep("str", hits_headers)
  hits_cols$exo_seq <- grep("^exo_seq$", hits_headers)
  hits_cols$exo_target <- grep("^exo_target$", hits_headers)
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

import_exome <- function(opt, file_exons, exome_cols) {
  
  exome <- fread(file_exons, skip = "Start")
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
  
  exome <- tibble(seqnames = seqnames, start = start, end = end, strand = strand, exo_symbol = symbol, exome[exome_cols$inherit]) %>%
    filter(!is.na(seqnames) & !is.na(start) & !is.na(end) & end >= start)
  
  if(opt$mode == "a" | opt$mode == "i") { # If in CRISPRi/a chemistry mode,
    exome <- exome %>%
      group_by(seqnames, exo_symbol, exome[exome_cols$inherit]) %>%
      reframe(seqnames, start = min(start) - 500, end = max(end) + 500, strand = "*", exo_symbol) %>% # enable upstream and downstream matches to default distance
      unique()
  } else if(is.numeric(opt$mode)) { # Else if in arbitrary chemistry mode,
    exome <- exome %>%
      group_by(seqnames, exo_symbol, exome[exome_cols$inherit]) %>% 
      reframe(seqnames, start = min(start) - as.numeric(opt$mode), end = max(end) + as.numeric(opt$mode), strand = "*", exo_symbol) %>% # enable upstream and downstream matches to arbitrary distance
      unique()
  }
  
  exome <- exome %>%
    mutate(ranges = paste(start, end, sep = "-"))
  
  exome <- GRanges(seqnames = exome$seqnames, ranges = exome$ranges, strand = exome$strand, exo_symbol = exome$exo_symbol, exome[exome_cols$inherit])
  
  return(exome)
}

import_hits <- function(file_genomic_ranges, hits_cols) {
  
  hits <- fread(file_genomic_ranges)
  seqnames <- hits[[as.numeric(hits_cols$chr)]]
  start <- hits[[as.numeric(hits_cols$start)]]
  end <- hits[[as.numeric(hits_cols$end)]]
  exo_seq = hits[[as.numeric(hits_cols$exo_seq)]]
  exo_target = hits[[as.numeric(hits_cols$exo_target)]]
  if (hits_cols$strand == "*") {                    # use * if in strandless mode
    strand <- rep("*", nrow(hits))
  } else {
    strand <- hits[[as.numeric(hits_cols$strand)]]
  }
  
  hits <- tibble(seqnames = seqnames, start = start, end = end, strand = strand, exo_seq = exo_seq, exo_target = exo_target) %>%
    filter(!is.na(seqnames) & !is.na(start) & !is.na(end) & end >= start) %>%
    mutate(ranges = paste(start, end, sep = "-"))
  hits <- GRanges(seqnames = hits$seqnames, ranges = hits$ranges, strand = hits$strand, exo_seq = as.character(hits$exo_seq), exo_target = as.character(hits$exo_target))
  return(hits)
}

runGrtem <- function(opt, blats) {
  
  if(file.not.exist.or.zero(blats$file_genomic_ranges_matched) | file.not.exist.or.zero(blats$file_genomic_ranges_distances)) {
    log_info("Determining exonic hits... ", blats$file_exons, "...")
    
    exome_cols <- inferExomeCols(blats$file_exons)                         # Infer column identities in the exome file
    hits_cols <- inferHitsCols(blats$file_genomic_ranges)                  # Infer column identities in the hits file
    
    exons <- import_exome(opt, blats$file_exons, exome_cols)                    # Make a gRanges object from the exome file (with Symbol column to provide re-annotations)
    hits <- import_hits(blats$file_genomic_ranges, hits_cols)                   # Make a gRanges object from the hits file (with ID column indicating guides in the library to be reannotated)
    
    exon_hits <- suppressWarnings(findOverlaps(hits, exons))                  # Make a gRanges Hits object indicating the pairs of gRanges that overlap between query (hits) and subject (exome)
    exon_hitsRanges <- suppressWarnings(findOverlapPairs(hits, exons))        # Make a gRanges Pairs object indicating the genomic ranges of pairs of gRanges that overlap between the first (hits) and second (exome)
    inherit <- grep("chr|seqname|exonStarts|exonEnds|str|symbol|name2|name$|cdsStart|cdsEnd|exonFrame|^ranges$|^seqlevels$|^seqlengths$|^isCircular$|^start$|^end$|^width$|^element$", names(as_tibble(exons)), value = T, invert = T)
    sgrna_hits <- GRanges(exon_hitsRanges@first,                              # Make a gRanges object which contains the genomic ranges of the hits and the Symbols from the exome
                          exo_symbol = exons$exo_symbol[exon_hits@to],
                          exo_seq = as.character(hits$exo_seq[exon_hits@from]),
                          exo_target = as.character(hits$exo_target[exon_hits@from]),
                          inherit = as_tibble(exons)[exon_hits@to, inherit])
    mapping <- as_tibble(sgrna_hits) %>%
      mutate(exo_cut = paste0(seqnames, ":", start)) %>%
      unique()
    mapping$assembly <- blats$file_genome
    mapping$exome <- blats$file_exons
    
    distances <- suppressWarnings(distanceToNearest(hits, exons, select = "arbitrary")) # Make a gRanges Hits object showing the distances from each hit (queryHits) to the nearest exon in the exome (subjectHits)
    sgrna_distances <- GRanges(hits[distances@from],                          # Make a gRanges object annotated with the gene of the nearest exon and the distance to that exon
                               nearestGene = exons$exo_symbol[distances@to],
                               distance = distances@elementMetadata$distance,
                               exo_seq = as.character(hits$exo_seq[distances@from]),
                               exo_target = as.character(hits$exo_target[distances@from]),
                               inherit = as_tibble(exons)[distances@to, inherit])
    
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
