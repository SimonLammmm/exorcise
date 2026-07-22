#### Main loop ####

reannotateLib <- function(opt) {
  
  #if(opt$adhoc) {
#    log_info("Using ad-hoc mode: genome and exome specified.")
    
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
                  file_base_edited = paste0(opt$outdir, "/exorcise.6-", sub(".+/(.+?)$", "\\1", opt$exome), "_baseEdited.tsv"),
                  file_exorcise_master_out = paste0(opt$outdir, "/exorcise.tsv"))
    
    # If specified, mask low expression genes
    if(length(opt$expression) > 0) {
      authors <- maskLowExpression(opt, authors)
    }
    
    extractGuides(opt, authors, blats)
    runBlat(opt, authors, blats)
    runPtgr(blats)
    runGrtem(opt, blats)
    
    # Additional behaviour for base editing
    if(opt$mode %in% c("cbe", "abe") | grepl("^be", opt$mode)) {
      baseEdit(opt, blats)
      exome_hits <- blats$file_base_edited
    } else {
      exome_hits <- blats$file_genomic_ranges_matched
    }
    
    # generate mapping
    genome_hits <- fread(blats$file_genomic_ranges, colClasses = "character")
    exome_hits <- fread(exome_hits, colClasses = "character")
    genome_hits <- genome_hits %>% filter(seqnames %in% unique(exome_hits$seqnames)) # ignore genome hits in chromosomes/variants not in the exome
    genome_hits <- genome_hits %>% transmute(exo_seq, exo_target, exo_cut)
    #exome_hits <- exome_hits %>% transmute(exo_seq, exo_cut, exo_symbol) %>% unique()
    all_mappings <- left_join(genome_hits, exome_hits, by = c("exo_seq", "exo_cut", "exo_target"), relationship = "many-to-many") %>% unique()
    
    premaster <- left_join(authors, all_mappings, by = "exo_seq", relationship = "many-to-many")
    master <- exorcisemaster(premaster, opt)
    
    outfile_m_tsv <- blats$file_exorcise_master_out
    dirwrite(master, outfile_m_tsv, sep = "\t")
    log_info("Wrote master library to ", outfile_m_tsv)
    
  #} else {
  #  log_info("Using post-hoc mode: exorcised library specified.")
  #  reannotateExisting(opt)
  #}
}