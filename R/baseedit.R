#### In silico base editing ####
#
# For each alignment, work out which bases a base editor could change, splice the
# affected transcripts in silico, and translate them to predict the protein-level
# consequence of every possible edit.


# All the editable bases in one window, one row per base.
find_editable_bases <- function(windows, opt) {
  windows <- windows %>%
    mutate(
      target_base = if_else(strand == "-", opt$be_from[[2]], opt$be_from[[1]]),
      new_base = if_else(strand == "-", opt$be_to[[2]], opt$be_to[[1]])
    )

  positions <- lapply(seq_len(nrow(windows)), function(i) {
    bases <- strsplit(windows$beWindowSeq[[i]], "", fixed = TRUE)[[1]]
    which(bases == windows$target_base[[i]])
  })

  editable <- windows[rep(seq_len(nrow(windows)), lengths(positions)), ]
  if (nrow(editable) == 0) {
    return(editable %>% mutate(be_position_in_window = integer(),
                               be_original = character(),
                               be_mutation = character(),
                               be_position = numeric()))
  }
  editable %>%
    mutate(
      be_position_in_window = unlist(positions, use.names = FALSE),
      be_original = substr(beWindowSeq, be_position_in_window, be_position_in_window),
      be_mutation = new_base,
      be_position = be_position_in_window + start - 1
    )
}


# Bystander edits: every combination of two or more bases a base editor could
# change at once within a single window. Each combination is reported as the span
# running from the first edited base to the last, so that the caller can treat it
# exactly like a single substitution.
find_bystander_edits <- function(sites) {
  n_sites <- nrow(sites)
  window <- sites$beWindowSeq[[1]]
  largest <- min(n_sites, nchar(window))
  if (largest < 2) {
    return(NULL)
  }

  combinations <- unlist(
    lapply(2:largest, function(k) utils::combn(n_sites, k, simplify = FALSE)),
    recursive = FALSE
  )

  bind_rows(lapply(combinations, function(index) {
    first <- sites$be_position_in_window[[min(index)]]
    last <- sites$be_position_in_window[[max(index)]]
    mutated <- window
    for (i in index) {
      position <- sites$be_position_in_window[[i]]
      substr(mutated, position, position) <- sites$be_mutation[[i]]
    }
    tibble(
      be_chr = sites$be_chr[[1]],
      exo_target = sites$exo_target[[1]],
      strand = sites$strand[[1]],
      beWindowSeq = window,
      be_position_in_window = first,
      be_original = substr(window, first, last),
      be_mutation = substr(mutated, first, last),
      be_position = first + sites$start[[1]] - 1
    )
  }))
}


# Read the base editing window of every alignment out of the genome.
read_edit_windows <- function(mapping, paths) {
  windows <- mapping %>%
    filter(!is.na(exo_target), !is.na(be_window_start), !is.na(be_window_end),
           be_window_end >= be_window_start)
  if (nrow(windows) == 0) {
    abort("No usable base editing window could be derived from the alignments.")
  }
  windows <- GRanges(
    seqnames = as.character(windows$seqnames),
    ranges = IRanges(start = windows$be_window_start, end = windows$be_window_end),
    strand = as.character(windows$strand),
    exo_target = as.character(windows$exo_target)
  )
  windows <- unique(windows)

  sequence <- as.character(import.2bit(resolve_path(paths$file_genome), which = windows))
  # import.2bit returns one base more than asked for; drop it to keep the
  # zero-based half-open convention used everywhere else.
  windows$be_window_seq <- sub("^.", "", sequence)
  windows
}


# Split every transcript's coding sequence out of the genome, restricted to the
# transcripts some guide could actually edit.
splice_affected_transcripts <- function(edit_windows, paths) {
  exome <- tibble::as_tibble(read_table(paths$file_exons, skip = "Starts"))
  names(exome) <- sub("^#", "", names(exome))

  column <- function(pattern, description) {
    name <- grep(pattern, names(exome), value = TRUE)
    if (length(name) == 0) {
      abort("Base editor mode needs a ", description, " column in ", paths$file_exons,
            " (looked for /", pattern, "/).")
    }
    name[[1]]
  }
  col <- list(
    chrom = column("chrom", "chromosome"),
    strand = column("strand", "strand"),
    exon_starts = column("exonStarts", "exon start"),
    exon_ends = column("exonEnds", "exon end"),
    exon_frames = column("exonFrame", "exon frame"),
    cds_start = column("cdsStart", "CDS start"),
    cds_end = column("cdsEnd", "CDS end"),
    transcript = column("name$", "transcript"),
    symbol = column("name2", "gene symbol")
  )

  exome <- exome %>%
    separate_longer_delim(
      all_of(c(col$exon_starts, col$exon_ends, col$exon_frames)),
      delim = ","
    ) %>%
    filter(nzchar(.data[[col$exon_frames]]))

  exons <- GRanges(
    seqnames = exome[[col$chrom]],
    ranges = IRanges(start = as.numeric(exome[[col$exon_starts]]),
                     end = as.numeric(exome[[col$exon_ends]])),
    strand = exome[[col$strand]],
    transcript = exome[[col$transcript]],
    name2 = exome[[col$symbol]],
    exonFrames = exome[[col$exon_frames]]
  )
  coding <- GRanges(
    seqnames = exome[[col$chrom]],
    ranges = IRanges(start = as.numeric(exome[[col$cds_start]]),
                     end = as.numeric(exome[[col$cds_end]])),
    strand = exome[[col$strand]]
  )

  # Intersecting each exon with its transcript's CDS bounds leaves only the
  # coding part of each exon. pintersect flags non-overlapping pairs, which are
  # wholly untranslated exons and of no use here.
  spliced <- pintersect(exons, coding)
  spliced <- spliced[mcols(spliced)$hit]

  affected <- unique(exons$transcript[to(findOverlaps(edit_windows, exons, ignore.strand = TRUE))])
  spliced <- spliced[spliced$transcript %in% affected]
  log_info(length(affected), " transcript(s) fall within a base editing window.")
  if (length(spliced) == 0) {
    abort("No coding exon overlaps any base editing window. Does ",
          paths$file_exons, " cover the right assembly?")
  }

  sequence <- as.character(import.2bit(resolve_path(paths$file_genome), which = spliced))
  spliced$cds <- sub("^.", "", sequence)
  spliced
}


run_base_edit <- function(opt, paths) {
  if (!needs_rebuild(paths$file_base_edited)) {
    log_info("Reusing the existing base edits in ", paths$file_base_edited, ".")
    return(invisible(paths$file_base_edited))
  }

  log_info("Predicting base edits: ", describe_mode(opt), "...")

  # Exon hits carry the annotation; the alignment coordinates carry the guide
  # bounds the editing window is measured from.
  mapping <- left_join(
    read_table(paths$file_exon_hits),
    read_table(paths$file_genomic_ranges),
    by = c("seqnames", "start", "end", "strand", "exo_seq", "exo_cut",
           "assembly", "exo_target")
  ) %>%
    as_tibble() %>%
    mutate(
      be_window_start = if_else(strand == "+",
                                guideBegin + opt$be_window_from,
                                guideFinal - opt$be_window_to),
      be_window_end = if_else(strand == "+",
                              guideBegin + opt$be_window_to,
                              guideFinal - opt$be_window_from)
    )
  if (nrow(mapping) == 0 || all(is.na(mapping$guideBegin))) {
    abort("Could not match the exon hits in ", paths$file_exon_hits,
          " to the alignment coordinates in ", paths$file_genomic_ranges,
          ". Delete both and let exorcise regenerate them.")
  }

  edit_windows <- read_edit_windows(mapping, paths)

  windows <- tibble(
    be_chr = as.character(seqnames(edit_windows)),
    exo_target = edit_windows$exo_target,
    strand = as.character(strand(edit_windows)),
    start = start(edit_windows),
    beWindowSeq = as.character(edit_windows$be_window_seq)
  ) %>% unique()

  single_edits <- find_editable_bases(windows, opt)
  log_info("Found ", nrow(single_edits), " single base edits across ",
           length(unique(single_edits$exo_target)), " alignment sites.")

  bystanders <- bind_rows(lapply(
    split(single_edits, single_edits$exo_target),
    find_bystander_edits
  ))
  log_info("Enumerated ", nrow(bystanders), " bystander edit combinations.")

  edits <- bind_rows(
    single_edits %>% select(-any_of(c("start", "target_base", "new_base"))),
    bystanders
  ) %>%
    filter(!is.na(be_position))
  if (nrow(edits) == 0) {
    abort("No editable base was found in any base editing window.")
  }

  edit_ranges <- GRanges(
    seqnames = edits$be_chr,
    ranges = IRanges(start = edits$be_position,
                     end = edits$be_position + nchar(edits$be_original)),
    strand = edits$strand,
    exo_target = edits$exo_target,
    be_original = edits$be_original,
    be_mutation = edits$be_mutation
  )

  spliced <- splice_affected_transcripts(edit_windows, paths)
  consequences <- predict_consequences(edit_ranges, spliced)

  result <- left_join(
    unique(mapping),
    consequences,
    by = "exo_target",
    relationship = "many-to-many"
  ) %>%
    mutate(exo_be = paste0(exo_symbol, ":", be_aa_mutation)) %>%
    unique()

  write_table(result, paths$file_base_edited, sep = "\t")
  log_info("Wrote ", nrow(result), " predicted edits to ", paths$file_base_edited, ".")
  invisible(paths$file_base_edited)
}


# Apply each edit to its transcript's coding sequence, translate, and describe
# what changed at the nucleotide and amino acid level.
predict_consequences <- function(edit_ranges, spliced) {
  spliced_table <- granges_tibble(spliced)

  overlaps <- findOverlaps(edit_ranges, spliced, ignore.strand = TRUE)
  query <- from(overlaps)
  subject <- to(overlaps)

  in_exon <- tibble(
    exo_target = edit_ranges$exo_target[query],
    be_original = edit_ranges$be_original[query],
    be_mutation = edit_ranges$be_mutation[query],
    be_chr = as.character(seqnames(edit_ranges))[query],
    be_strand = as.character(strand(edit_ranges))[query],
    be_position = start(edit_ranges)[query],
    transcript = spliced$transcript[subject],
    strand = as.character(strand(spliced))[subject],
    name2 = spliced$name2[subject],
    exonFrame = spliced$exonFrames[subject],
    exonStart = start(spliced)[subject],
    exonEnd = end(spliced)[subject],
    cds = as.character(spliced$cds)[subject]
  )

  # The exon sequence itself is only needed to confirm the exon is not empty; the
  # edit is applied to the assembled transcript further down, not here.
  in_exon <- in_exon %>%
    mutate(be_relative_position_in_exon = be_position - exonStart) %>%
    # Drop edits that findOverlaps admitted only because the edit range is one
    # base wider than the edit itself.
    filter(be_relative_position_in_exon >= 0,
           be_relative_position_in_exon < exonEnd - exonStart,
           nzchar(cds)) %>%
    select(-cds)
  if (nrow(in_exon) == 0) {
    abort("Every predicted edit fell outside the coding exons it was assigned to.")
  }

  # Position of each coding exon within its transcript's spliced CDS.
  coding_exons <- spliced_table %>%
    filter(exonFrames != -1) %>%
    group_by(transcript) %>%
    arrange(start, .by_group = TRUE) %>%
    mutate(widthBefore = pmax(0, c(0, cumsum(width - 1))[seq_len(dplyr::n())])) %>%
    ungroup()

  transcripts <- coding_exons %>%
    arrange(start) %>%
    summarise(flseq = paste0(cds, collapse = ""), .by = c(transcript, name2, strand))

  # unique() because a transcript duplicated in the exome would otherwise
  # multiply every edit assigned to it.
  exon_offsets <- coding_exons %>%
    transmute(exonStart = start, transcript, widthBefore) %>%
    unique()

  edited <- in_exon %>%
    left_join(exon_offsets, by = c("exonStart", "transcript"),
              relationship = "many-to-many") %>%
    mutate(be_relative_position_in_cds = be_relative_position_in_exon + widthBefore) %>%
    left_join(transcripts, by = c("transcript", "name2", "strand")) %>%
    filter(!is.na(flseq), !is.na(be_relative_position_in_cds),
           be_relative_position_in_cds < nchar(flseq))
  if (nrow(edited) == 0) {
    abort("No predicted edit could be placed within a spliced coding sequence.")
  }

  edited$be_flseq <- edited$flseq
  substr(edited$be_flseq,
         edited$be_relative_position_in_cds + 1,
         edited$be_relative_position_in_cds + nchar(edited$be_original)) <- edited$be_mutation

  # Transcripts on the minus strand are read in the other direction.
  transcripts <- transcripts %>%
    mutate(flseq = if_else(strand == "-", revcomp(flseq), flseq),
           translated = translate_cds(flseq))

  edited <- edited %>%
    mutate(
      be_flseq = if_else(strand == "-", revcomp(be_flseq), be_flseq),
      be_translated = translate_cds(be_flseq),
      # Amino acid positions are 1-based, unlike the nucleotide positions above.
      be_position_in_aa = if_else(
        strand == "-",
        floor((nchar(be_flseq) - 1 - be_relative_position_in_cds) / 3 + 1),
        floor(be_relative_position_in_cds / 3) + 1
      )
    ) %>%
    select(-flseq) %>%
    left_join(transcripts %>% select(transcript, name2, strand, translated),
              by = c("transcript", "name2", "strand"))

  # Comparing the translations directly is more reliable than predicting which
  # codon moved, because a multi-base bystander edit can shift several.
  edited <- bind_cols(edited, difference_range(edited$translated, edited$be_translated))
  edited %>%
    mutate(
      be_original_in_aa = if_else(
        !is.na(diff_start),
        substr(translated, diff_start, diff_start + diff_length - 1),
        substr(translated, be_position_in_aa,
               be_position_in_aa + floor(nchar(be_original) / 3))
      ),
      be_modified_in_aa = if_else(
        !is.na(diff_start),
        substr(be_translated, diff_start, diff_start + diff_length - 1),
        be_original_in_aa
      ),
      be_position_in_aa = if_else(!is.na(diff_start), as.numeric(diff_start),
                                  as.numeric(be_position_in_aa)),
      be_consequence = case_when(
        be_original_in_aa == be_modified_in_aa ~ "synonymous",
        be_modified_in_aa == "*" & be_original_in_aa != "*" ~ "stopgain",
        be_modified_in_aa != "*" & be_original_in_aa == "*" ~ "stoploss",
        be_modified_in_aa != be_original_in_aa ~ "missense"
      ),
      be_nt_mutation = if_else(
        strand == "-",
        paste0(revcomp(be_original), nchar(be_flseq) - be_relative_position_in_cds,
               revcomp(be_mutation)),
        paste0(be_original, be_relative_position_in_cds + 1, be_mutation)
      ),
      be_aa_mutation = paste0(be_original_in_aa, be_position_in_aa, be_modified_in_aa),
      be_position_in_genome = paste0(be_chr, ":", be_position, "-",
                                     be_position + nchar(be_original))
    ) %>%
    transmute(exo_target, transcript, transcript_strand = strand,
              be_position_in_genome, be_nt_mutation, be_aa_mutation, be_consequence)
}


# translate() warns for every sequence whose length is not a multiple of three,
# which is most of them once an exome has partial CDS records.
translate_cds <- function(x) {
  suppressWarnings(as.character(Biostrings::translate(
    Biostrings::DNAStringSet(x), if.fuzzy.codon = "solve"
  )))
}


# First and last positions at which two equal-length strings differ. Returns NA
# where the strings match or cannot be compared.
difference_range <- function(a, b) {
  n <- length(a)
  first <- rep(NA_integer_, n)
  last <- rep(NA_integer_, n)
  for (i in seq_len(n)) {
    if (is.na(a[[i]]) || is.na(b[[i]]) || nchar(a[[i]]) != nchar(b[[i]])) {
      next
    }
    differs <- which(strsplit(a[[i]], "", fixed = TRUE)[[1]] !=
                       strsplit(b[[i]], "", fixed = TRUE)[[1]])
    if (length(differs) == 0) {
      next
    }
    first[[i]] <- differs[[1]]
    last[[i]] <- differs[[length(differs)]]
  }
  tibble(diff_start = first, diff_end = last, diff_length = last - first + 1L)
}
