#!/usr/bin/env Rscript

# exorcise-trace: plot the proportion of each nucleotide called at each
# sequencing cycle in a FASTQ file. Useful for spotting the fixed-sequence
# regions of a CRISPR amplicon and confirming where the variable guide region
# starts.
#
# Was `ntByCycle`, which still works but warns.
#
# Accepts plain and gzipped FASTQ.
#
# Author:  Dr Simon Lam, University of Cambridge <sl681@cam.ac.uk>
# Source:  https://github.com/SimonLammmm/exorcise
# Licence: Creative Commons Zero v1.0 Universal

TRACE_VERSION <- "1.0.0"

suppressWarnings(suppressMessages({
  library(optparse)
  library(logger)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(plotly)
}))

FASTQ_PATTERN <- "\\.(fastq|fq)(\\.gz)?$"


#### Reading ####

# A FASTQ record is four lines and the sequence is the second, so read the file
# as lines and take every fourth. Earlier versions guessed at which lines were
# sequences by testing for [ATCGN] only, which can also match a quality line.
read_fastq_sequences <- function(path, max_lines = Inf) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, n = if (is.finite(max_lines)) as.integer(max_lines) else -1L,
                     warn = FALSE)
  if (length(lines) < 2) {
    return(character())
  }
  toupper(lines[seq(2L, length(lines), by = 4L)])
}


#### Counting ####

# Proportion of each base called at each cycle. substring() is vectorised over
# the reads, so one pass per cycle is enough; the previous implementation built
# the table row by row and then patched in the zero counts with a join.
nucleotide_frequencies <- function(reads, first_cycle, last_cycle) {
  if (length(reads) == 0) {
    return(NULL)
  }
  longest <- max(nchar(reads))
  last_cycle <- min(longest, last_cycle)
  if (first_cycle > last_cycle) {
    log_warn("Requested cycles ", first_cycle, " to ", last_cycle,
             " but the longest read is only ", longest, " bases. Skipping.")
    return(NULL)
  }

  cycles <- first_cycle:last_cycle
  counts <- bind_rows(lapply(cycles, function(cycle) {
    called <- substring(reads, cycle, cycle)
    tallied <- table(called[nzchar(called)])
    tibble(Cycle = cycle, Base = names(tallied), Freq = as.numeric(tallied))
  }))
  if (nrow(counts) == 0) {
    return(NULL)
  }

  counts %>%
    mutate(Pct = Freq / length(reads)) %>%
    complete(Cycle = cycles, Base, fill = list(Freq = 0, Pct = 0)) %>%
    arrange(Cycle, Base)
}


#### Plotting ####

nucleotide_plot <- function(frequencies, title) {
  ggplot(frequencies, aes(x = Cycle, y = Pct, colour = Base)) +
    geom_line() +
    scale_x_continuous(minor_breaks = seq(1, max(frequencies$Cycle), 1)) +
    labs(title = title, x = "Cycle", y = "Proportion of reads") +
    theme_classic() +
    theme(
      panel.grid.major = element_line(colour = "black", linewidth = 0.1),
      panel.grid.minor = element_line(colour = "grey", linewidth = 0.1)
    )
}


analyse_fastq <- function(path, outdir, max_lines, first_cycle, last_cycle) {
  label <- basename(path)
  log_info("Reading ", path, "...")
  reads <- read_fastq_sequences(path, max_lines)
  if (length(reads) == 0) {
    log_warn("No sequences found in ", path, ". Skipping.")
    return(invisible(NULL))
  }
  log_info("Read ", length(reads), " sequences. Counting nucleotides by cycle...")

  frequencies <- nucleotide_frequencies(reads, first_cycle, last_cycle)
  if (is.null(frequencies)) {
    log_warn("No nucleotide calls to plot for ", path, ". Skipping.")
    return(invisible(NULL))
  }

  plot <- nucleotide_plot(frequencies, label)

  table_file <- file.path(outdir, paste0(label, ".tsv"))
  utils::write.table(frequencies, table_file, sep = "\t", quote = FALSE, row.names = FALSE)
  log_info("Wrote counts to ", table_file, ".")

  static_file <- file.path(outdir, paste0(label, ".pdf"))
  pdf(file = static_file, width = 7, height = 7)
  print(plot)
  dev.off()
  log_info("Wrote the static plot to ", static_file, ".")

  write_interactive_plot(plot, file.path(outdir, paste0(label, ".html")))

  invisible(frequencies)
}


# A self-contained widget is one file you can email; producing one needs pandoc,
# which htmlwidgets shells out to. pandoc is a declared dependency, so needing
# the fallback means something is wrong with the installation, but a plot with a
# sidecar directory beats no plot at all.
write_interactive_plot <- function(plot, path) {
  widget <- config(ggplotly(plot), showLink = TRUE)

  tryCatch({
    htmlwidgets::saveWidget(widget, path, selfcontained = TRUE)
    log_info("Wrote the interactive plot to ", path, ".")
  }, error = function(e) {
    log_warn("Could not write a self-contained interactive plot: ",
             conditionMessage(e))
    tryCatch({
      htmlwidgets::saveWidget(widget, path, selfcontained = FALSE)
      log_warn("Wrote the interactive plot to ", path, " instead, with its ",
               "JavaScript in a sibling directory. Install pandoc for a single ",
               "self-contained file: it ships in env/exorcise.yaml.")
    }, error = function(e2) {
      log_warn("Could not write the interactive plot at all: ",
               conditionMessage(e2), ". The static plot and the counts table ",
               "were still written.")
    })
  })
}


#### Arguments ####

trace_option_list <- function() {
  list(
    make_option(c("-f", "--file"), type = "character", default = NULL,
                metavar = "PATHS",
                help = paste("Required. Comma-separated FASTQ files or directories.",
                             "Directories are searched recursively for .fastq,",
                             ".fq and their gzipped forms.")),
    make_option(c("-o", "--out"), type = "character", default = ".",
                metavar = "DIR",
                help = "Directory to write results into. [default %default]"),
    make_option(c("-n", "--nrows"), type = "numeric", default = Inf,
                metavar = "N",
                help = paste("Read only the first N lines of each file. Four lines",
                             "make one read. [default all]")),
    make_option(c("-s", "--start"), type = "numeric", default = 1,
                metavar = "N",
                help = "First cycle to plot. [default %default]"),
    make_option(c("-e", "--end"), type = "numeric", default = Inf,
                metavar = "N",
                help = "Last cycle to plot. [default the longest read]"),
    make_option("--version", action = "store_true", default = FALSE,
                help = "Print the version and exit.")
  )
}


# Expand the --file argument into a list of FASTQ files.
collect_fastq_files <- function(inputs) {
  found <- character()
  for (input in inputs) {
    if (dir.exists(input)) {
      log_info("Searching ", input, " for FASTQ files...")
      found <- c(found, list.files(input, pattern = FASTQ_PATTERN,
                                   recursive = TRUE, full.names = TRUE))
    } else if (!file.exists(input)) {
      log_warn(input, " does not exist. Ignoring.")
    } else if (grepl(FASTQ_PATTERN, input)) {
      found <- c(found, input)
    } else {
      log_warn(input, " is not a FASTQ file. Ignoring.")
    }
  }
  unique(found)
}


main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_args(OptionParser(option_list = trace_option_list()), args = args)

  if (isTRUE(opt$version)) {
    cat(TRACE_VERSION, "\n", sep = "")
    return(invisible(0L))
  }

  outdir <- opt$out
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(outdir)) {
    stop("Could not create the output directory ", outdir,
         ". Check the path and permissions.", call. = FALSE)
  }

  logfile <- file.path(outdir, paste0("logfile_exorcise-trace_",
                                      format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log"))
  log_appender(appender_tee(logfile))
  log_info("exorcise-trace ", TRACE_VERSION, " starting.")

  if (is.null(opt$file)) {
    stop("--file is required. Run exorcise-trace --help for the syntax.", call. = FALSE)
  }
  files <- collect_fastq_files(trimws(unlist(strsplit(opt$file, ","))))
  if (length(files) == 0) {
    stop("No FASTQ files were found in --file.", call. = FALSE)
  }
  log_info("Analysing ", length(files), " file(s):\n", paste(files, collapse = "\n"))

  started <- proc.time()
  for (file in files) {
    analyse_fastq(file, outdir, opt$nrows, opt$start, opt$end)
  }
  log_info("exorcise-trace finished in ", round((proc.time() - started)[["elapsed"]], 1),
           " seconds.")
  invisible(0L)
}


if (!interactive()) {
  main()
}
