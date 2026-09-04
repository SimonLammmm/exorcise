#### Shared helpers ####
#
# Small utilities used across the exorcise pipeline. Nothing here knows about
# the pipeline itself, so this file is safe to source first.


# Columns in an exome file that exorcise consumes directly. Everything else is
# inherited verbatim onto the output. This pattern used to be copy-pasted in
# three places, which meant a change in one only half-worked.
EXOME_RESERVED_PATTERN <- paste0(
  "chr|seqname|exonStarts|exonEnds|str|symbol|name2|name$|cdsStart|cdsEnd|",
  "exonFrame|^ranges$|^seqlevels$|^seqlengths$|^isCircular$|^start$|^end$|",
  "^width$|^element$|^hit$"
)

# Names of the columns that should be passed through untouched.
inherited_cols <- function(x) {
  headers <- if (is.character(x)) x else names(x)
  grep(EXOME_RESERVED_PATTERN, headers, value = TRUE, invert = TRUE)
}


#### Paths and I/O ####

# Historically exorcise was called with paths that omitted the compression
# suffix (`foo.tsv` for `foo.tsv.gz`), so resolve those before reading.
resolve_path <- function(path) {
  if (file.exists(path)) {
    return(path)
  }
  candidates <- Sys.glob(paste0(path, "*"))
  if (length(candidates) == 0) {
    return(path) # let the caller report the missing file
  }
  candidates[[1]]
}

is_readable_file <- function(path) {
  path <- resolve_path(path)
  file.exists(path) && !dir.exists(path)
}

# data.table::fread hands gzipped files to R.utils, which it neither declares nor
# loads. Most exorcise inputs are gzipped, so check once, up front, rather than
# letting fread fail obscurely part way through a run.
require_gzip_support <- function() {
  if (requireNamespace("R.utils", quietly = TRUE)) {
    return(invisible(TRUE))
  }
  stop(
    "The R.utils package is not installed. data.table::fread needs it to read ",
    "gzipped files, and most exorcise inputs are gzipped. Install it with ",
    'install.packages("R.utils"), or recreate the conda environment from ',
    "env/exorcise.yaml.",
    call. = FALSE
  )
}

read_table <- function(path, ...) {
  data.table::fread(file = resolve_path(path), ...)
}

write_table <- function(x, path, ...) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, ...)
  invisible(path)
}

# TRUE when a checkpoint is missing or empty and so has to be regenerated.
needs_rebuild <- function(path) {
  path <- resolve_path(path)
  !file.exists(path) || isTRUE(file.size(path) == 0)
}

# Read a FASTA into a plain character vector of sequences, joining records that
# were wrapped across several lines.
read_fasta <- function(path) {
  lines <- readLines(resolve_path(path), warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) {
    return(character())
  }
  is_header <- startsWith(lines, ">")
  record <- cumsum(is_header)[!is_header]
  body <- lines[!is_header]
  if (length(body) == 0) {
    return(character())
  }
  parts <- split(body, factor(record, levels = sort(unique(record))))
  unname(vapply(parts, paste0, character(1), collapse = ""))
}

write_fasta <- function(names, seqs, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste0(">", names, "\n", seqs), con)
  invisible(path)
}


#### Sequences and ranges ####

# Reverse complement. Biostrings does this in C and, unlike the chartr() version
# this replaces, handles IUPAC ambiguity codes correctly.
revcomp <- function(x) {
  if (length(x) == 0) {
    return(character())
  }
  as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(toupper(x))))
}

# GRanges to tibble, without relying on an as_tibble() method being registered
# for S4 objects by whichever package happened to be attached.
#
# as.data.frame() returns seqnames and strand as factors, which then refuse to
# join against the character vectors that come out of as.character(seqnames(x))
# elsewhere. Coercing them here keeps every table in the pipeline comparable.
granges_tibble <- function(x) {
  out <- tibble::as_tibble(as.data.frame(x))
  for (column in intersect(c("seqnames", "strand"), names(out))) {
    out[[column]] <- as.character(out[[column]])
  }
  out
}


#### Errors ####

# Log and abort. Every fatal path in exorcise goes through here so that the
# message reaches both the console and the logfile.
abort <- function(...) {
  message <- paste0(...)
  log_error(message)
  stop(paste0("FATAL: ", message), call. = FALSE)
}

require_binary <- function(bin, env_var) {
  if (nzchar(Sys.which(bin)) || is_readable_file(bin)) {
    return(invisible(TRUE))
  }
  abort(
    "Could not find `", bin, "` on the PATH. Install it, or point ", env_var,
    " at the executable."
  )
}
