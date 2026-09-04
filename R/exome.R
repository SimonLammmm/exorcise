#### Mapping alignments onto the exome ####


# Exon overlap is deliberately strandless: a guide annealing to the antisense
# strand still disrupts the gene it cuts. Set this to the exome's strand column
# if you ever need sense-only matching.
EXOME_MATCH_STRAND <- "*"


# Locate the columns exorcise needs in an exome file. UCSC Table Browser output
# is the reference format, but the names vary enough between downloads that they
# have to be matched by pattern rather than by position.
exome_columns <- function(headers) {
  pick <- function(pattern, description) {
    index <- grep(pattern, headers)
    if (length(index) == 0) {
      abort("The exome has no ", description, " column (looked for /", pattern, "/).")
    }
    index[[1]]
  }
  list(
    chrom = pick("chr|seqname", "chromosome"),
    start = pick("exonStarts", "exon start"),
    end = pick("exonEnds", "exon end"),
    symbol = pick("symbol|name2", "gene symbol"),
    inherit = inherited_cols(headers)
  )
}


# UCSC packs every exon of a transcript into one comma-separated row. Split it so
# that each exon becomes its own interval.
melt_exons <- function(exome, cols) {
  start_col <- names(exome)[[cols$start]]
  end_col <- names(exome)[[cols$end]]
  exome %>%
    separate_longer_delim(all_of(c(start_col, end_col)), delim = ",") %>%
    filter(nzchar(.data[[start_col]]), nzchar(.data[[end_col]])) %>%
    unique()
}


import_exome <- function(opt, path) {
  exome <- tibble::as_tibble(read_table(path, skip = "Starts"))
  names(exome) <- sub("^#", "", names(exome))
  cols <- exome_columns(names(exome))

  if (any(grepl(",", exome[[cols$start]]))) {
    exome <- melt_exons(exome, cols)
  }

  exome <- tibble(
    seqnames = as.character(exome[[cols$chrom]]),
    start = suppressWarnings(as.numeric(exome[[cols$start]])),
    end = suppressWarnings(as.numeric(exome[[cols$end]])),
    strand = EXOME_MATCH_STRAND,
    exo_symbol = sub(",$", "", as.character(exome[[cols$symbol]])),
    exome[cols$inherit]
  ) %>%
    filter(!is.na(seqnames), !is.na(start), !is.na(end), end >= start)

  if (is_proximity_mode(opt$mode)) {
    # CRISPRi/a and arbitrary-distance chemistries do not cut, so collapse each
    # gene to a single interval and pad it. Padding is clamped at 1 because a
    # gene near the start of a contig would otherwise get a nonsensical range.
    flank <- proximity_flank(opt$mode)
    log_info("Collapsing exons to gene bodies padded by ", flank, " bases...")
    grouping <- c("seqnames", "exo_symbol", cols$inherit)
    exome <- exome %>%
      group_by(across(all_of(grouping))) %>%
      summarise(start = max(1, min(start) - flank), end = max(end) + flank,
                .groups = "drop") %>%
      mutate(strand = EXOME_MATCH_STRAND) %>%
      relocate(seqnames, start, end, strand, exo_symbol)
  }

  log_info("Loaded ", nrow(exome), " intervals covering ",
           length(unique(exome$exo_symbol)), " symbols from ", path, ".")

  metadata <- exome[c("exo_symbol", cols$inherit)]
  GRanges(
    seqnames = exome$seqnames,
    ranges = IRanges(start = exome$start, end = exome$end),
    strand = exome$strand,
    metadata
  )
}


import_hits <- function(path) {
  hits <- read_table(path)
  required <- c("seqnames", "start", "end", "strand", "exo_seq", "exo_target")
  missing <- setdiff(required, names(hits))
  if (length(missing) > 0) {
    abort(path, " is missing the column(s) ", paste(missing, collapse = ", "),
          ". Delete it and let exorcise regenerate it.")
  }
  hits <- hits %>%
    as_tibble() %>%
    filter(!is.na(seqnames), !is.na(start), !is.na(end), end >= start)
  # The alignment's own strand is kept, because downstream stages join the exon
  # hits back onto the alignment coordinates on it. Strandlessness is applied at
  # the exome end instead, via EXOME_MATCH_STRAND.
  GRanges(
    seqnames = as.character(hits$seqnames),
    ranges = IRanges(start = hits$start, end = hits$end),
    strand = as.character(hits$strand),
    exo_seq = as.character(hits$exo_seq),
    exo_target = as.character(hits$exo_target)
  )
}


# Two outputs: the exons each alignment lands in, and, for alignments that miss
# every exon, the distance to the nearest one. The latter is useful for
# diagnosing a library that was designed against a different assembly.
map_to_exome <- function(opt, paths) {
  if (!needs_rebuild(paths$file_exon_hits) && !needs_rebuild(paths$file_exon_distances)) {
    log_info("Reusing the existing exon hits in ", paths$file_exon_hits, ".")
    return(invisible(paths$file_exon_hits))
  }

  log_info("Intersecting alignments with ", paths$file_exons, " (mode: ",
           describe_mode(opt), ")...")
  exons <- import_exome(opt, paths$file_exons)
  hits <- import_hits(paths$file_genomic_ranges)

  exon_table <- granges_tibble(exons)
  inherit <- inherited_cols(exon_table)

  # Inherited exome columns are prefixed so that an arbitrary exome column can
  # never silently collide with one of exorcise's own.
  attach_inherited <- function(target, subject_index) {
    if (length(inherit) == 0) {
      return(target)
    }
    inherited <- exon_table[subject_index, inherit, drop = FALSE]
    names(inherited) <- paste0("inherit.", names(inherited))
    bind_cols(target, inherited)
  }

  overlaps <- suppressWarnings(findOverlaps(hits, exons))
  mapping <- granges_tibble(hits[from(overlaps)]) %>%
    mutate(exo_symbol = exons$exo_symbol[to(overlaps)]) %>%
    attach_inherited(to(overlaps)) %>%
    mutate(exo_cut = paste0(seqnames, ":", start),
           assembly = paths$file_genome,
           exome = paths$file_exons) %>%
    unique()

  nearest <- suppressWarnings(distanceToNearest(hits, exons, select = "arbitrary"))
  distances <- granges_tibble(hits[from(nearest)]) %>%
    mutate(nearestGene = exons$exo_symbol[to(nearest)],
           distance = mcols(nearest)$distance) %>%
    attach_inherited(to(nearest)) %>%
    mutate(exo_cut = paste0(seqnames, ":", start),
           assembly = paths$file_genome,
           exome = paths$file_exons) %>%
    unique()

  write_table(mapping, paths$file_exon_hits, sep = "\t")
  write_table(distances, paths$file_exon_distances, sep = "\t")
  log_info("Annotated ", length(unique(mapping$exo_target)), " of ",
           length(unique(hits$exo_target)), " alignment sites from the exome.")
  invisible(paths$file_exon_hits)
}
