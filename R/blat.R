#### Alignment against the genome ####


# BLAT and twoBitToFa are taken from the PATH. The Docker image points these
# environment variables at the binaries inside its conda environment, which is
# why this file no longer needs to be rewritten with sed at build time.
blat_binary <- function() Sys.getenv("EXORCISE_BLAT", unset = "blat")
twobit_binary <- function() Sys.getenv("EXORCISE_TWOBITTOFA", unset = "twoBitToFa")

# The 21 columns of a PSL record, in order. The first five lines of a PSL file
# are a two-row header and a rule, none of which parse as data.
PSL_COLUMNS <- c(
  "match", "mismatch", "rep_match", "n_count", "q_gap_count", "q_gap_bases",
  "t_gap_count", "t_gap_bases", "strand", "q_name", "q_size", "q_start", "q_end",
  "t_name", "t_size", "t_start", "t_end", "block_count", "block_sizes",
  "q_starts", "t_starts"
)
PSL_HEADER_LINES <- 5L


# Align the queries to the genome. Only perfect, full-length, ungapped hits are
# of interest, so the minimum score is the full length of the shortest query.
run_blat <- function(opt, authors, paths) {
  if (!needs_rebuild(paths$file_psl)) {
    log_info("Reusing the existing alignments in ", paths$file_psl, ".")
    return(invisible(paths$file_psl))
  }

  bin <- blat_binary()
  require_binary(bin, "EXORCISE_BLAT")

  # Ambiguous PAM positions are not sent to BLAT, so they do not count towards
  # the score either.
  pam_bases <- nchar(gsub("[^ACGT]", "", opt$pam))
  min_score <- min(nchar(authors$exo_seq)) + pam_bases

  args <- c(
    shQuote(resolve_path(paths$file_genome)),
    shQuote(paths$file_sgRNAs),
    shQuote(paths$file_psl),
    "-stepSize=4", "-tileSize=10", "-fine", "-repMatch=2000000",
    paste0("-minScore=", min_score), "-minIdentity=100"
  )
  dir.create(dirname(paths$file_psl), recursive = TRUE, showWarnings = FALSE)

  log_info("Running BLAT. Depending on the library and genome this can take ",
           "anywhere from minutes to days: ", bin, " ", paste(args, collapse = " "))
  status <- system2(bin, args)
  if (!identical(as.integer(status), 0L)) {
    abort("BLAT exited with status ", status, ".")
  }
  invisible(paths$file_psl)
}


read_psl <- function(path) {
  psl <- try(
    read_table(path, skip = PSL_HEADER_LINES, header = FALSE, fill = TRUE),
    silent = TRUE
  )
  if (inherits(psl, "try-error") || nrow(psl) == 0) {
    return(NULL)
  }
  if (ncol(psl) != length(PSL_COLUMNS)) {
    abort(path, " does not look like a PSL file: expected ", length(PSL_COLUMNS),
          " columns but found ", ncol(psl), ".")
  }
  psl <- setNames(tibble::as_tibble(psl), PSL_COLUMNS)
  # PSL stores block lists as trailing-comma strings; strip them so the numeric
  # columns coerce cleanly.
  psl %>%
    mutate(across(everything(), ~ sub(",$", "", as.character(.x)))) %>%
    mutate(across(all_of(c("q_size", "q_start", "q_end", "t_start", "t_end", "block_count")),
                  ~ suppressWarnings(as.integer(.x))))
}


# Turn alignments into genomic coordinates, then read the aligned sequence back
# out of the genome to confirm that each hit really is the query.
locate_alignments <- function(opt, paths) {
  if (!needs_rebuild(paths$file_genomic_ranges)) {
    log_info("Reusing the existing alignment coordinates in ", paths$file_genomic_ranges, ".")
    return(invisible(paths$file_genomic_ranges))
  }

  log_info("Reading alignments from ", paths$file_psl, "...")
  psl <- read_psl(paths$file_psl)
  if (is.null(psl)) {
    abort("BLAT found no alignments between ", opt$infile, " and ", opt$genome,
          ". Check that --seq points at the sequence column, that --pam is right, ",
          "and that --genome is the assembly you meant.")
  }

  total <- nrow(psl)
  psl <- psl %>%
    filter(!is.na(q_end), q_end == q_size, block_count == 1L, q_start == 0L) %>%
    filter(strand %in% c("+", "-"))
  log_info("Kept ", nrow(psl), " of ", total, " alignments that are perfect, ",
           "ungapped and full length.")
  if (nrow(psl) == 0) {
    abort("None of the ", total, " alignments in ", paths$file_psl,
          " is a perfect full-length match.")
  }

  # Cas9 cuts three to four bases 5' of the PAM. Ranges follow the convention
  # already established by the intermediate files: guideBegin is the PSL target
  # start and guideFinal the PSL target end.
  cut_offset <- 3L + nchar(opt$pam)
  pam_width <- nchar(opt$pam)

  hits <- psl %>%
    transmute(
      seqnames = t_name,
      guideBegin = t_start,
      guideFinal = t_end,
      strand = strand,
      assembly = paths$file_genome
    ) %>%
    mutate(
      start = if_else(strand == "+", guideFinal - cut_offset, guideBegin + cut_offset),
      end = start,
      seqSpec = if_else(
        strand == "+",
        paste0(seqnames, ":", guideBegin, "-", guideFinal - pam_width),
        paste0(seqnames, ":", guideBegin + pam_width, "-", guideFinal)
      ),
      exo_target = paste0(seqnames, ":", guideBegin, "-", guideFinal, "_", strand),
      exo_cut = paste0(seqnames, ":", start)
    )

  hits$exo_seq <- fetch_aligned_sequences(hits, paths)
  report_alignment_coverage(opt, hits, paths)

  write_table(hits, paths$file_genomic_ranges, sep = "\t")
  log_info("Wrote ", nrow(hits), " alignment coordinates to ", paths$file_genomic_ranges, ".")
  invisible(paths$file_genomic_ranges)
}


fetch_aligned_sequences <- function(hits, paths) {
  bin <- twobit_binary()
  require_binary(bin, "EXORCISE_TWOBITTOFA")

  dir.create(dirname(paths$file_genomic_seqSpecs), recursive = TRUE, showWarnings = FALSE)
  writeLines(hits$seqSpec, paths$file_genomic_seqSpecs)

  log_info("Reading ", nrow(hits), " aligned sequences back out of the genome...")
  status <- system2(bin, c(
    shQuote(resolve_path(paths$file_genome)),
    paste0("-seqList=", shQuote(paths$file_genomic_seqSpecs)),
    shQuote(paths$file_genomic_seqs)
  ))
  if (!identical(as.integer(status), 0L)) {
    abort("twoBitToFa exited with status ", status, ".")
  }

  seqs <- read_fasta(paths$file_genomic_seqs)
  if (length(seqs) != nrow(hits)) {
    abort("twoBitToFa returned ", length(seqs), " sequences for ", nrow(hits),
          " alignments. ", paths$file_genomic_seqs, " and ",
          paths$file_genomic_seqSpecs, " have got out of step.")
  }
  if_else(hits$strand == "+", toupper(seqs), revcomp(seqs))
}


# A mismatch here almost always means a checkpointed .psl was carried over from a
# run against a different input file.
report_alignment_coverage <- function(opt, hits, paths) {
  queries <- unique(read_fasta(paths$file_sgRNAs))
  queries <- substr(queries, 1L, nchar(queries) - nchar(opt$pam))
  aligned <- unique(hits$exo_seq)

  matched <- sum(aligned %in% queries)
  log_info("BLAT aligned ", matched, " of ", length(queries), " distinct query sequences.")

  unexpected <- length(aligned) - matched
  if (unexpected > 0) {
    log_warn(paths$file_psl, " contains ", unexpected, " sequence(s) that are not ",
             "in ", opt$infile, ". This is what a stale checkpoint looks like: ",
             "delete the intermediate files if the input has changed. Continuing.")
  }
}
