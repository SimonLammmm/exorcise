#### Command-line options ####
#
# Parsing, validation and normalisation of exorcise's arguments. Validation
# collects every problem it can find before reporting, so a user with three
# mistakes in their command learns about all three in one run.


#### Base editor chemistry ####

# from/to are given as (sense, antisense) pairs so that a guide on either strand
# can be edited without re-deriving the complement at each use site.
BE_PRESETS <- list(
  cbe = list(from = c("C", "G"), to = c("T", "A"), window = c(2L, 8L)),
  abe = list(from = c("A", "T"), to = c("G", "C"), window = c(4L, 9L))
)

# Custom base editor specification: be<FROM><TO><START><END>, where START and
# END are two digits each, e.g. "beCT0208" is equivalent to "cbe".
parse_be_spec <- function(spec) {
  match <- regmatches(
    toupper(spec),
    regexec("^BE([ACGT])([ACGT])([0-9]{2})([0-9]{2})$", toupper(spec))
  )[[1]]
  if (length(match) != 5) {
    return(NULL)
  }
  window <- as.integer(match[4:5])
  if (window[[1]] > window[[2]] || window[[1]] < 1) {
    return(NULL)
  }
  list(
    from = c(match[[2]], chartr("ACGT", "TGCA", match[[2]])),
    to = c(match[[3]], chartr("ACGT", "TGCA", match[[3]])),
    window = window
  )
}

is_base_edit_mode <- function(mode) {
  is.character(mode) && (mode %in% names(BE_PRESETS) || grepl("^be", mode))
}

# CRISPRi/a and arbitrary-distance modes annotate a guide with any feature within
# this many bases of the target site, rather than only features it cuts inside.
is_proximity_mode <- function(mode) {
  is.numeric(mode) || (is.character(mode) && mode %in% c("a", "i"))
}

proximity_flank <- function(mode) {
  if (is.numeric(mode)) as.numeric(mode) else 500
}


#### Parsing ####

exorcise_option_list <- function() {
  list(
    make_option(c("-i", "--infile"), type = "character", default = NULL,
                metavar = "FILE",
                help = "Required. File containing the sequences to exorcise."),
    make_option(c("-o", "--outdir"), type = "character", default = NULL,
                metavar = "DIR",
                help = "Required. Directory to write results into. Created if absent."),
    make_option(c("-g", "--seq"), type = "character", default = NULL,
                metavar = "N",
                help = "Required. 1-based column number in --infile holding the sequences."),
    make_option(c("-v", "--genome"), type = "character", default = NULL,
                metavar = "FILE",
                help = "Required. Genome assembly in 2bit format."),
    make_option(c("-w", "--exome"), type = "character", default = NULL,
                metavar = "FILE",
                help = "Required. Exome annotation in UCSC Table Browser format."),
    make_option(c("-z", "--pam"), type = "character", default = NULL,
                metavar = "SEQ",
                help = "PAM sequence appended to each guide before alignment, e.g. NGG. [ACGTN]."),
    make_option(c("-q", "--mode"), type = "character", default = "ko",
                metavar = "MODE",
                help = paste("CRISPR chemistry: ko (default), a, i, cbe, abe,",
                             "beXYnnmm for a custom base editor, or an integer",
                             "distance in bases.")),
    make_option(c("-n", "--harm"), type = "character", default = NULL,
                metavar = "N",
                help = "1-based column number in --infile holding the existing annotations."),
    make_option(c("-c", "--control"), type = "character", default = NULL,
                metavar = "PATTERNS",
                help = paste("Comma-separated regular expressions matching control",
                             "guides in the --harm column.")),
    make_option(c("-d", "--control_type"), type = "character", default = NULL,
                metavar = "NAMES",
                help = "Comma-separated control names, one per --control pattern."),
    make_option(c("-x", "--expression"), type = "character", default = NULL,
                metavar = "FILE",
                help = paste("Two-column file of gene symbol and expression value.",
                             "Guides against genes below --exprcutoff are dropped.")),
    make_option(c("-k", "--exprcutoff"), type = "character", default = "10",
                metavar = "VALUE",
                help = "Expression value below which a gene counts as not expressed. [default %default]"),
    make_option("--ref", action = "store_true", default = FALSE,
                help = "Print the citation and exit."),
    make_option("--version", action = "store_true", default = FALSE,
                help = "Print the version and exit.")
  )
}

parse_exorcise_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  parse_args(OptionParser(option_list = exorcise_option_list()), args = args)
}

commasplit <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  values <- trimws(unlist(strsplit(as.character(x), ",")))
  values[nzchar(values)]
}


#### Validation ####

# Options that only ever take one value. Historically every one of these had its
# own near-identical block of "if more than one, take the first and warn".
SCALAR_OPTIONS <- c(
  "infile", "outdir", "seq", "pam", "mode", "genome", "exome", "harm",
  "expression", "exprcutoff"
)

EXOME_REQUIRED_COLUMNS <- c("chrom", "strand", "exonStarts", "exonEnds", "name2")
EXOME_BE_COLUMNS <- c("name$", "cdsStart", "cdsEnd", "exonFrame")

validate_options <- function(opt) {
  # Environments have reference semantics, so err() and warn() can accumulate
  # messages without a return value being threaded through every check. assign()
  # rather than report$x <- ... because this is load-bearing: if the messages were
  # silently dropped, a bad command would run instead of being rejected.
  report <- new.env(parent = emptyenv())
  assign("errors", character(), envir = report)
  assign("warnings", character(), envir = report)
  err <- function(...) assign("errors", c(report$errors, paste0(...)), envir = report)
  warn <- function(...) assign("warnings", c(report$warnings, paste0(...)), envir = report)

  for (flag in c(SCALAR_OPTIONS, "control", "control_type")) {
    opt[[flag]] <- commasplit(opt[[flag]])
  }
  for (flag in SCALAR_OPTIONS) {
    if (length(opt[[flag]]) > 1) {
      warn("--", flag, " takes one value but received ", length(opt[[flag]]),
           ". Using the first only: ", opt[[flag]][[1]], ".")
      opt[[flag]] <- opt[[flag]][[1]]
    }
  }

  # --infile
  infile_headers <- character()
  if (length(opt$infile) == 0) {
    err("--infile is required but was not given.")
  } else if (!is_readable_file(opt$infile)) {
    err("--infile ", opt$infile, " is not a readable file.")
    opt$infile <- NULL
  } else {
    infile_headers <- names(read_table(opt$infile, nrows = 1L))
    if (any(startsWith(infile_headers, "exo_"))) {
      warn("--infile ", opt$infile, " already has exo_ columns. ",
           "These will be overwritten in the output.")
    }
  }

  # --outdir
  if (length(opt$outdir) == 0) {
    err("--outdir is required but was not given.")
  } else if (!dir.exists(opt$outdir)) {
    dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(opt$outdir)) {
      err("--outdir ", opt$outdir, " could not be created.")
    }
  }

  # --seq and --harm are 1-based column numbers into --infile
  n_cols <- if (length(infile_headers) > 0) length(infile_headers) else NULL
  column_index <- function(flag) {
    value <- opt[[flag]]
    if (!grepl("^[0-9]+$", value)) {
      err("--", flag, " must be a positive integer, got ", value, ".")
      return(NULL)
    }
    value <- as.integer(value)
    if (value < 1) {
      err("--", flag, " is 1-based, got ", value, ".")
      return(NULL)
    }
    if (!is.null(n_cols) && value > n_cols) {
      err("--", flag, " is out of bounds: got ", value, " but --infile has ",
          n_cols, " columns.")
      return(NULL)
    }
    value
  }

  if (length(opt$seq) == 0) {
    err("--seq is required but was not given.")
  } else {
    opt$seq <- column_index("seq")
  }
  if (length(opt$harm) > 0) {
    opt$harm <- column_index("harm")
  }
  if (!is.null(opt$seq) && !is.null(opt$harm) && identical(opt$seq, opt$harm)) {
    err("--seq and --harm both point at column ", opt$seq,
        ". They must be different columns.")
  }

  # --pam
  if (length(opt$pam) == 0) {
    opt$pam <- ""
    warn("--pam was not given. Aligning the sequences as supplied, with no PAM.")
  } else if (grepl("[^ACGTNacgtn]", opt$pam)) {
    err("--pam ", opt$pam, " contains non-nucleotide characters. Only [ACGTN] are accepted.")
  } else {
    opt$pam <- toupper(opt$pam)
  }

  # --mode
  mode_result <- resolve_mode(opt$mode, warn)
  opt <- utils::modifyList(opt, mode_result)

  # --genome
  if (length(opt$genome) == 0) {
    err("--genome is required but was not given.")
  } else if (!is_readable_file(opt$genome)) {
    err("--genome ", opt$genome, " is not a readable file.")
  }

  # --exome
  if (length(opt$exome) == 0) {
    err("--exome is required but was not given.")
  } else if (!is_readable_file(opt$exome)) {
    err("--exome ", opt$exome, " is not a readable file.")
  } else {
    headers <- names(read_table(opt$exome, nrows = 0L, skip = "Starts"))
    required <- EXOME_REQUIRED_COLUMNS
    if (is_base_edit_mode(opt$mode)) {
      required <- c(required, EXOME_BE_COLUMNS)
    }
    missing <- required[!vapply(required, function(p) any(grepl(p, headers)), logical(1))]
    if (length(missing) > 0) {
      err("--exome ", opt$exome, " is missing required column(s): ",
          paste(sub("\\$$", "", missing), collapse = ", "),
          ". Exome files must be in UCSC Table Browser format",
          if (is_base_edit_mode(opt$mode)) " and, in base editor mode, include CDS coordinates." else ".")
    }
  }

  # --control and --control_type
  if (length(opt$control) > 0 && length(opt$harm) == 0) {
    warn("--control needs --harm to know which column to search. Ignoring --control.")
    opt$control <- NULL
    opt$control_type <- NULL
  }
  if (length(opt$control) == 0 && length(opt$control_type) > 0) {
    warn("--control_type was given without --control. Ignoring --control_type.")
    opt$control_type <- NULL
  }
  if (length(opt$control) > 0) {
    if (length(opt$control_type) == 0) {
      opt$control_type <- paste0("Non-targeting", seq_along(opt$control))
      warn("--control_type was not given. Naming the control types ",
           paste(opt$control_type, collapse = ", "), ".")
    } else if (length(opt$control_type) == 1 && length(opt$control) > 1) {
      opt$control_type <- rep(opt$control_type, length(opt$control))
    } else if (length(opt$control_type) != length(opt$control)) {
      err("--control has ", length(opt$control), " pattern(s) but --control_type has ",
          length(opt$control_type), " name(s). They must match.")
    }
    if (NON_TARGETING_PREFIX %in% opt$control_type) {
      err("--control_type may not be \"", NON_TARGETING_PREFIX,
          "\"; exorcise uses that name for guides it could not annotate.")
    }
    if (any(duplicated(opt$control_type))) {
      err("--control_type contains duplicated names: ",
          paste(unique(opt$control_type[duplicated(opt$control_type)]), collapse = ", "), ".")
    }
    invalid <- opt$control[!vapply(opt$control, is_valid_regex, logical(1))]
    if (length(invalid) > 0) {
      err("--control contains invalid regular expression(s): ",
          paste(invalid, collapse = ", "), ".")
    }
  }

  # --expression and --exprcutoff
  cutoff <- suppressWarnings(as.numeric(opt$exprcutoff))
  if (is.na(cutoff)) {
    warn("--exprcutoff ", opt$exprcutoff, " is not a number. Using the default of 10.")
    cutoff <- 10
  }
  opt$exprcutoff <- cutoff

  if (length(opt$expression) > 0) {
    if (!is_readable_file(opt$expression)) {
      warn("--expression ", opt$expression, " is not a readable file. Ignoring it.")
      opt$expression <- NULL
    } else if (length(opt$harm) == 0) {
      warn("--expression needs --harm to know which symbols to match. Ignoring --expression.")
      opt$expression <- NULL
    }
  } else {
    opt$expression <- NULL
  }

  for (message in report$warnings) log_warn(message)
  if (length(report$errors) > 0) {
    for (message in report$errors) log_error(message)
    stop("FATAL: ", length(report$errors), " problem(s) with the arguments. Quitting.",
         call. = FALSE)
  }

  opt
}

is_valid_regex <- function(pattern) {
  !inherits(try(grepl(pattern, "", perl = FALSE), silent = TRUE), "try-error")
}

# Turn the --mode string into the mode itself plus, for base editors, the window
# and substitution it implies.
resolve_mode <- function(mode, warn) {
  if (length(mode) == 0 || !nzchar(mode)) {
    return(list(mode = "ko"))
  }

  # An integer is an arbitrary proximity distance in bases.
  if (grepl("^[0-9]+$", mode)) {
    return(list(mode = as.numeric(mode)))
  }

  mode <- tolower(mode)
  if (mode %in% c("ko", "a", "i")) {
    return(list(mode = mode))
  }

  if (mode %in% names(BE_PRESETS)) {
    preset <- BE_PRESETS[[mode]]
    return(list(mode = mode, be_from = preset$from, be_to = preset$to,
                be_window_from = preset$window[[1]], be_window_to = preset$window[[2]]))
  }

  if (startsWith(mode, "be")) {
    spec <- parse_be_spec(mode)
    if (is.null(spec)) {
      warn("--mode ", mode, " is not a valid custom base editor specification ",
           "(expected beXYnnmm, e.g. beCT0208). Falling back to cbe.")
      spec <- BE_PRESETS$cbe
      mode <- "cbe"
    }
    return(list(mode = mode, be_from = spec$from, be_to = spec$to,
                be_window_from = spec$window[[1]], be_window_to = spec$window[[2]]))
  }

  warn("--mode ", mode, " is not a recognised chemistry. Falling back to ko.")
  list(mode = "ko")
}

describe_mode <- function(opt) {
  if (is.numeric(opt$mode)) {
    return(paste0("proximity, within ", opt$mode, " bases of the target site"))
  }
  switch(opt$mode,
    ko = "CRISPR knockout, annotating the cut site",
    a = "CRISPR activation, within 500 bases of the target site",
    i = "CRISPR interference, within 500 bases of the target site",
    paste0("base editor, ", opt$be_from[[1]], " to ", opt$be_to[[1]],
           " in window [", opt$be_window_from, ", ", opt$be_window_to, "]")
  )
}
