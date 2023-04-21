#!/usr/bin/env Rscript
#
# ntByCycle.R
#
# Visualise % nucleotide by cycle in sequencing FASTQ files
# Accepts FASTQ and gzipped FASTQ files
#
# Author: Dr Simon Lam
# Affiliation: University of Cambridge
# Contact: sl681@cam.ac.uk
#
# Version history:
# SL    25-Jan-2023     0.1.0     Initial release
# SL    05-Apr-2023     0.1.1     Fixed multi-file checker, added nrows functionality, fixed for "multi-column" FASTQ files
#
#

#### Init ####
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(plotly)
  library(optparse)
  library(logger)
  library(foreach)
})

ver <- "0.1.0"

#### Main loop ####

ntByCycle <- function(infiles, outdir, nrows) {
  
  # Analyse one infile at a time
  for (f in 1:length(infiles)) {
    
    # Obtain friendly infile name
    infile_friendly <- sub(".+\\/(.+)$", "\\1", infiles[f])
    
    # Read infile
    log_info("Reading file ", f, " of ", length(infiles), ": ", infiles[f], "...")
    fastq <- fread(infiles[f], header = F, nrows = nrows, fill = TRUE)                                      # Read infile
    fastq <- fastq %>% filter(!grepl("[^ATCGN]+", V1)) %>% select(V1)                          # Keep only sequence lines
    
    # Determine number of reads
    reads <- length(fastq$V1)
    
    # Determine number of cycles
    cycles <- max(sapply(fastq, function(x) nchar(x)))
    
    # Populate nucleotide frequencies by cycle number
    log_info("Determining nucleotide frequencies...")
    freq <- data.frame()
    for (c in 1:cycles) {
      freq <- rbind(freq, data.frame(Cycle = c, table(sapply(fastq, function(x) substr(x, c,c)))))
    }
    
    # Determine nucleotide percentages by cycle number
    freq <- freq %>% transmute(Cycle, Base = Var1, Freq, Pct = Freq / reads)    # We assume that exactly one call was made per cycle in each read
    
    # Plot
    log_info("Plotting...")
    p <- ggplot(freq, aes(x = Cycle, y = Pct, colour = Base)) +
      geom_line() +
      ggtitle(infile_friendly)
    
    #q <- ggplotly(p)
    
    # Export
    outfile <- paste0(outdir, "/", infile_friendly, ".pdf")
    pdf(file = outfile, width = 7, height = 7)
    print(p)
    dev.off()
    log_info("Exported static plot to ", outfile, ".")
    
    #outfile2 <- paste0(outdir, "/", infile_friendly, ".html")
    #htmlwidgets::saveWidget(config(q, showLink = T), outfile2)
    #log_info("Exported interactive plot to ", outfile2, ".")
    
  }
  
  
}

#### Arguments ####

if (!interactive()) {                                                           # Take command-line inputs if we're running on the command-line
  option_list <- list(
    make_option(opt_str = c("-f", "--file"), type = "character", default = NULL,
                help = "Path to FASTQ file(s). Accepts gzipped FASTQ files. If a directory is supplied, then all files ending with .fastq and .fastq.gz will be accepted.", metavar = "character"),
    make_option(opt_str = c("-o", "--out"), type = "character", default = ".",
                help = "Optional. Path to directory to output the result. Defaults to the current working directory.", metavar = "character"),
    make_option(opt_str = c("-n", "--nrows"), type = "numeric", default = Inf,
                help = "Optional. How many rows to read from the top of the FASTQ file(s). Default Inf.", metavar = "numeric")
  )
  
  opt_parser = OptionParser(option_list = option_list)
  opt = parse_args(opt_parser)
  
} else {                                                                        # If we're running in an IDE, take inputs from this block instead
  opt <- NULL
  opt$file <- "/Volumes/files/research/sjlab/archive/fs01/Sequencing/NovaSeq FASTQ/NVS044-329912606/FASTQ_Generation_2022-02-15_13_02_45Z-528623110/SPJms2050_d9_1_L001-ds.42b4dfdc570040c383437762aeb36e89/SPJms2050-d9-1_S1_L001_R1_001.fastq.gz"
  opt$out <- "~/bio/Sandbox/ntByCycle"
  opt$nrows <- 10
}

#### Error checking ####

panic <- function() stop("There were errors while running the script. Check the logfile for details.", call. = F)

# Try to create the output folder
if(!file.exists(opt$out)) { dir.create(opt$out, recursive = T, showWarnings = F) }

# Stop if output folder couldn't be created
if(!file.exists(opt$out)) {
  stop("Invalid argument to --out: ", opt$out, ". Could not create directory. Maybe a permissions issue?")
}

outdir <- opt$out

# Start logging
logfile <- paste0(opt$out, "/logfile_ntByCycle_", format(Sys.time(), "%Y-%m-%dT%H-%M-%S%Z"), ".log")
log_appender(appender_tee(logfile))
start_time <- proc.time()
log_info("Welcome to ntByCycle, version ", ver, ".")

# Check if the --file is valid
if(is.null(opt$file)) {
  log_error("Nothing supplied to --file. Run ntByCycle.R --help for syntax. Exiting.")
  panic()
}

opt$file <- strsplit(opt$file, ",")[[1]]


# Check if the --file exists
if (any(!file.exists(opt$file))) {
  fileNotFound <- opt$file[which(!file.exists(opt$file))]
  log_error("Invalid argument to --file: ", as.character(fileNotFound), ". File/directory not found.")
  panic()
}

# Populate infiles by going through --file inputs
infiles <- NULL

for (i in 1:length(opt$file)) {
  # If --file exists, then check if it is a directory
  if (file.exists(opt$file[i]) & dir.exists(opt$file[i])) {
    # If --file is a directory, find all .fastq and .fastz.gz files recursively and add them to the infiles list
    log_info("Accepted --file as a directory: ", opt$file[i], ".")
    infiles <- c(infiles, list.files(opt$file[i], pattern = "\\.fastq$|\\.fastq\\.gz", recursive = T, full.names = T))
  } else {
    # If --file is a file, then check if it's a fastq file. If so, add it; if not, ignore it and warn
    if (grepl("\\.fastq$|\\.fastq\\.gz$", opt$file[i])) {
      log_info("Accepted --file as a file: ", opt$file[i], ".")
      infiles <- c(infiles, opt$file[i])
    }
    else {
      log_warn("Supplied --file in not a .fastq or .fastq.gz file: ", opt$file[i], ". Ignoring.")
    }
  }
}

# Check if zero infiles were found
if (length(infiles) == 0) {
  log_error("No fastq files were found within --files argument. Exiting.")
  panic()
}

# Ask the user if it's OK to go
log_info("Found these fastq files:\n", paste(infiles, collapse = "\n"))

# Warn the user that it'll take a long time if more than one infile was selected
if (length(infiles) > 1) {
  log_warn("You have selected ", length(infiles), " fastq files to analyse with ntByCycle. It is recommended to analyse one file at a time. Analysing too many files will take a long time. Are you sure you want to continue?")

# Wait for user to authorise
  Sys.sleep(1)
  invisible(readline(prompt = "Press [enter] to continue."))
}

nrows <- opt$nrows

#### Execution ####

ntByCycle(infiles, outdir, nrows)
end_time <- proc.time()
log_info("crispieR-cts process completed in ", (end_time - start_time)[[3]], " seconds.")
