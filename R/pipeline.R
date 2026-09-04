#### Pipeline ####


# Every path exorcise reads or writes, in one place. Names are unchanged from
# earlier versions so that intermediate files from a previous run are still
# picked up as checkpoints.
exorcise_paths <- function(opt) {
  genome <- basename(opt$genome)
  exome <- basename(opt$exome)
  stage <- function(...) file.path(opt$outdir, paste0(...))
  list(
    file_genome = opt$genome,
    file_exons = opt$exome,
    file_sgRNAs = stage("exorcise.1-seq.fa"),
    file_psl = stage("exorcise.2-", genome, "_BLAT.psl"),
    file_genomic_ranges = stage("exorcise.3-", genome, "_genomicRanges.tsv"),
    file_genomic_seqSpecs = stage("exorcise.3-", genome, "_genomicSeqSpecs.tsv"),
    file_genomic_seqs = stage("exorcise.3-", genome, "_genomicSeqs.fa"),
    file_exon_hits = stage("exorcise.4-", exome, "_exonHits.tsv"),
    file_exon_distances = stage("exorcise.5-", exome, "_exonDist.tsv"),
    file_base_edited = stage("exorcise.6-", exome, "_baseEdited.tsv"),
    file_out = stage("exorcise.tsv")
  )
}


# Stage 1 to 6 each write a checkpoint and are skipped if that checkpoint is
# already present and non-empty. Only the final join is always recomputed, since
# it is cheap next to the alignment.
run_exorcise <- function(opt) {
  paths <- exorcise_paths(opt)

  authors <- import_input(opt)
  if (!is.null(opt$expression)) {
    authors <- mask_low_expression(opt, authors)
  }

  write_guide_fasta(opt, authors, paths)
  run_blat(opt, authors, paths)
  locate_alignments(opt, paths)
  map_to_exome(opt, paths)

  annotations <- if (is_base_edit_mode(opt$mode)) {
    run_base_edit(opt, paths)
    paths$file_base_edited
  } else {
    paths$file_exon_hits
  }

  master <- join_annotations(authors, annotations, paths)
  master <- finalise_output(master, opt)

  write_table(master, paths$file_out, sep = "\t")
  log_info("Wrote the exorcised library to ", paths$file_out, ".")
  invisible(paths$file_out)
}


# Bring the alignments and their exome annotations back onto the user's rows.
# Everything is read as character so that identifiers built by string
# concatenation upstream join reliably.
join_annotations <- function(authors, annotations_file, paths) {
  genome_hits <- read_table(paths$file_genomic_ranges, colClasses = "character")
  exome_hits <- read_table(annotations_file, colClasses = "character")

  # An alignment on a contig the exome says nothing about cannot be annotated
  # from it, so drop it rather than carry a blank symbol through the join. This is
  # what keeps alternative haplotypes and unplaced scaffolds out of the output.
  known <- unique(exome_hits$seqnames)
  dropped <- sum(!genome_hits$seqnames %in% known)
  if (dropped > 0) {
    log_info("Ignoring ", dropped, " alignments on ",
             length(setdiff(unique(genome_hits$seqnames), known)),
             " contigs that the exome does not cover.")
  }

  genome_hits <- genome_hits %>%
    as_tibble() %>%
    filter(seqnames %in% known) %>%
    transmute(exo_seq, exo_target, exo_cut)

  mappings <- left_join(
    genome_hits, as_tibble(exome_hits),
    by = c("exo_seq", "exo_cut", "exo_target"),
    relationship = "many-to-many"
  ) %>% unique()

  left_join(authors, mappings, by = "exo_seq", relationship = "many-to-many")
}
