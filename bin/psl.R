#### Send and parse BLAT ####

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
    
    psl <- suppressWarnings(fread(infile, fill = T))
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
    twoBitToFa_command <- "twoBitToFa"
    system(paste0(twoBitToFa_command, " ", blats$file_genome, " -seqList=", blats$file_genomic_seqSpecs, " ", blats$file_genomic_seqs))
    
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
    
    nBoth <- length(which(unique(out$exo_seq) %in% unique(blatIn$seq)))
    nBlatIn <- length(unique(blatIn$seq))
    nPsl <- length(unique(out$exo_seq))
    nPslOnly <- nPsl - nBoth
    
    
    log_info("BLAT aligned ", nBoth, " of ", nBlatIn, " sequences given.")
    if (nPslOnly > 0) {
      log_warn("BLAT results ", blats$file_psl, " aligned ", nPslOnly," sequences that were not found in the infile. Did you checkpoint with the correct .psl file? Continuing.")
    }
    
    dir.create(dirname(outfile), recursive = T, showWarnings = F)							
    fwrite(out, outfile, sep = "\t") 
  } else {
    log_info("Not recalculating genomic hit coordinates: results already exist, ", blats$file_genomic_ranges, ".")
  }
  
  return()
}
