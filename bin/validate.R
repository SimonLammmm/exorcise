#### Validate inputs ####

fixOpts <- function(opt) {
  
  # fix comma-separated argument inputs
  if(!is.null(opt$infile))       opt$infile       <- commasplit(opt$infile)
  if(!is.null(opt$outdir))       opt$outdir       <- commasplit(opt$outdir)
  if(!is.null(opt$seq))          opt$seq          <- commasplit(opt$seq)
  if(!is.null(opt$pam))          opt$pam          <- commasplit(opt$pam)
  #if(!is.null(opt$mode))         opt$mode         <- commasplit(opt$mode)
  #if(!is.null(opt$library))      opt$library      <- commasplit(opt$library)
  if(!is.null(opt$genome))       opt$genome       <- commasplit(opt$genome)
  if(!is.null(opt$exome))        opt$exome        <- commasplit(opt$exome)
  #if(!is.null(opt$priorities))   opt$priorities   <- commasplit(opt$priorities)
  #if(!is.null(opt$id))           opt$id           <- commasplit(opt$id)
  if(!is.null(opt$harm))         opt$harm         <- commasplit(opt$harm)
  if(!is.null(opt$control))      opt$control      <- commasplit(opt$control)
  if(!is.null(opt$control_type)) opt$control_type <- commasplit(opt$control_type)
  
  opt$adhoc <- F
  
  # check arguments
  errors <- list()
  warnings <- list()
  
  # check if number of inputs are correct
  lengthChecker <- function(got, need, label) {
    got <- length(got)
    need <- length(need)
    if(got != 0 & got != 1 & got != need) {
      return(paste0("Invalid number of inputs to --", label, ". Should be 1 or ", need, ", but got ", got, ". Using first value only."))
    } else {
      return(NULL)
    }
  }
  
  # check infile
  if(length(opt$infile) > 0) {
    if(length(opt$infile) > 1) {
      opt$infile <- opt$infile[1]
      warnings <- c(warnings, paste0("Warning: --infile received more than one argument. Using first value only: ", opt$infile, "."))
    }
    if(!file.exists(opt$infile)) {
      opt$infile_glob <- Sys.glob(paste0(opt$infile, "*"))
      if(length(opt$infile_glob) == 0) {
        errors <- c(errors, paste0("Error: --infile ", opt$infile, " not found."))
        opt$infile <- NULL
      } else {
        opt$infile <- opt$infile_glob
        warnings <- c(warnings, paste0("Warning: --infile accepted as globbed argument ", opt$infile, "."))
      }
    }
    if(!isFile(opt$infile)) {
      errors <- c(errors, paste0("Error: --infile ", opt$infile, " is not a file."))
      opt$infile <- NULL
    }
  } else {
    errors <- c(errors, paste0("Error: --infile not specified."))
  }
  
  # check infile headers
  if(length(opt$infile) > 0) {
    opt$infile_headers <- names(fread(opt$infile, nrows = 1))
    if (any(grepl("^exo_", opt$infile_headers))) {
      warnings <- c(warnings, paste0("Warning: --infile ", opt$infile, " contains exorcise-like column names (\"exo_\"). These might be clobbered."))
    }
  }
  
  # check outdir
  if(length(opt$outdir) > 0) {
    if(length(opt$outdir) > 1) {
      opt$outdir <- opt$outdir[1]
      warnings <- c(warnings, paste0("Warning: --outdir received more than one argument. Using first value only: ", opt$outdir, "."))
    }
    if(!dir.exists(opt$outdir)) {
      dir.create(opt$outdir, showWarnings = F, recursive = T)
      warnings <- c(warnings, paste0("Info: --outdir ", opt$outdir, " being recursively created."))
    }
  } else {
    errors <- c(errors, paste0("Error: --outdir not specified."))
  }
  
  # check seq
  if(length(opt$seq) > 0) {
    if(length(opt$seq) > 1) {
      opt$seq <- opt$seq[1]
      warnings <- c(warnings, paste0("Warning: --seq received more than one argument. Using first value only: ", opt$seq, "."))
    }
    if(!grepl("^\\d+$", opt$seq)) {
      errors <- c(errors, paste0("Error: --seq must be an integer, got ", opt$seq, "."))
    } else {
      opt$seq <- as.numeric(opt$seq)
      if (length(opt$infile) > 0) {
        if (opt$seq > length(opt$infile_headers)) {
          errors <- c(errors, paste0("Error: --seq is out of bounds, got ", opt$seq, " but --infile only contains ", length(opt$infile_headers), " columns."))
        }
      }
    }
  } else {
    errors <- c(errors, paste0("Error: --seq not specified."))
  }
  
  # check library
  if(length(opt$library) > 0) {
    warnings <- c(warnings, paste0("Warning: deprecated option --library specified. Ignoring."))
    # if(length(opt$library) > 1) {
    #   opt$library <- opt$library[1]
    #   warnings <- c(warnings, paste0("Warning: --library received more than one argument. Using first value only: ", opt$library, "."))
    # }
    # if(!file.exists(opt$library)) {
    #   opt$library_glob <- Sys.glob(paste0(opt$library, "*"))
    #   if(length(opt$library_glob) == 0) {
    #     errors <- c(errors, paste0("Error: --library ", opt$library, " not found."))
    #   } else {
    #     opt$library <- opt$library_glob
    #     warnings <- c(warnings, paste0("Warning: --library accepted as globbed argument ", opt$library))
    #   }
    # }
    # if(file.exists(opt$library)) {
    #   opt$library_headers <- names(fread(opt$library, nrows = 0))
    #   if(!("exo_id" %in% opt$library_headers)) {
    #     errors <- c(errors, paste0("Error: --library ", opt$library, " doesn't look like an exorcised library. Expected an `exo_id` column."))
    #   }
    #   if(!("exo_seq" %in% opt$library_headers)) {
    #     errors <- c(errors, paste0("Error: --library ", opt$library, " doesn't look like an exorcised library. Expected an `exo_seq` column."))
    #   }
    #   if(!("exo_symbol" %in% opt$library_headers)) {
    #     errors <- c(errors, paste0("Error: --library ", opt$library, " doesn't look like an exorcised library. Expected an `exo_symbol` column."))
    #   }
    # }
  } #else {
    opt$adhoc <- T
    #warnings <- c(warnings, paste0("Info: --library not specified. Using ad-hoc mode."))
  #}
  
  # check pam
  if(length(opt$pam) > 0) {
    if(opt$adhoc) {
      if(length(opt$pam) > 1) {
        opt$pam <- opt$pam[1]
        warnings <- c(warnings, paste0("Warning: --pam received more than one argument. Using first value only: ", opt$pam, "."))
      }
      if(grepl("[^ACTGNactgn]", opt$pam)) {
        errors <- c(errors, paste0("Error: --pam ", opt$pam," contains non-nucleotide letters. Only [ACTGN] are accepted."))
      }
    } else {
      warnings <- c(warnings, paste0("Warning: --pam specified outside of ad-hoc mode. Ignoring"))
    }
  } else if(opt$adhoc) {
    opt$pam <- ""
    warnings <- c(warnings, paste0("Warning: --pam not specified. Not appending any PAM."))
  }
  
  # check mode
  if(length(opt$mode) > 0) {
    if(length(opt$mode) > 1) {
      opt$mode <- opt$mode[1]
      warnings <- c(warnings, paste0("Warning: --mode received more than one argument. Using first value only: ", opt$mode, "."))
    }
    opt$mode_tolower <- tolower(opt$mode)
    if(grepl("^\\d+$", opt$mode)) {  # Arbitrary chemistry mode, accept hits this many positions from first/last exon bounds
      opt$mode <- as.numeric(opt$mode)
    } else if(!(opt$mode_tolower %in% c("ko", "a", "i", "cbe", "abe")) & !(grepl("^be", opt$mode_tolower))) { # Standard chemistry modes
      warnings <- c(warnings, paste0("Warning: --mode received an invalid value: ", opt$mode, ". Falling back to KO."))
      opt$mode <- "ko"
    } else {
      opt$mode <- opt$mode_tolower
      if (opt$mode == "cbe") { # Apply cytosine base editor settings
        opt$be_window_from = 2
        opt$be_window_to = 8
        opt$be_from = c("C", "G") # C on the forward strand, G on the reverse strand
        opt$be_to = c("T", "A")
      } else if (opt$mode == "abe") { # Apply adenine base editor settings
        opt$be_window_from = 4
        opt$be_window_to = 9
        opt$be_from = c("A", "T") # A on the forward strand, T on the reverse strand
        opt$be_to = c("G", "C")
      } else if (grepl("^be", opt$mode)) { # Check and apply custom base editor settings
        if (toupper(substr(opt$mode, 3, 3)) %in% c("A", "C", "G", "T") &
            toupper(substr(opt$mode, 4, 4)) %in% c("A", "C", "G", "T") &
            grepl("\\d{2}", substr(opt$mode, 5, 6)) &
            grepl("\\d{2}", substr(opt$mode, 7, 8))
        ) {
          opt$be_from = toupper(substr(opt$mode, 3, 3))
          opt$be_to = toupper(substr(opt$mode, 4, 4))
          opt$be_from = c(opt$be_from, chartr("ACGT", "TGCA", opt$be_from))
          opt$be_to = c(opt$be_to, chartr("ACGT", "TGCA", opt$be_to))
          opt$be_window_from = as.numeric(substr(opt$mode, 5, 6))
          opt$be_window_to = as.numeric(substr(opt$mode, 7, 8))
          log_info("Applying custom BE settings: original ", opt$be_from[1], "; edit ", opt$be_to[1], "; window start ", opt$be_window_from, "; window end ", opt$be_window_end, ".")
        } else { # Fall back to CBE if fail
          warnings <- c(warnings, paste0("Warning: --mode received an invalid value: ", opt$mode, ". Falling back to CBE"))
          opt$be_window_from = 2
          opt$be_window_to = 8
          opt$be_from = c("C", "G") # C on the forward strand, G on the reverse strand
          opt$be_to = c("T", "A")
        }
      }
    }
  } else {
    opt$mode <- "ko"
    warnings <- c(warnings, paste0("Info: --mode not specified. Using CRISPRko chemistry."))
  }
  
  # check genome
  if(length(opt$genome) > 0) {
    if(length(opt$genome) > 1) {
      opt$genome <- opt$genome[1]
      warnings <- c(warnings, paste0("Warning: --genome received more than one argument. Using first value only: ", opt$genome, "."))
    }
    if(!file.exists(opt$genome)) {
      opt$genome_glob <- Sys.glob(paste0(opt$genome, "*"))
      if(length(opt$genome_glob) == 0) {
        errors <- c(errors, paste0("Error: --genome ", opt$genome, " not found."))
      } else {
        opt$genome <- opt$genome_glob
        warnings <- c(warnings, paste0("Warning: --genome accepted as globbed argument ", opt$genome))
      }
    }
  # } else if(opt$adhoc) {
  #   errors <- c(errors, paste0("Error: --genome not passed while in ad-hoc mode."))
  } else {
    errors <- c(errors, paste0("Error: --genome not specified."))
  }
  
  # check exome
  if(length(opt$exome) > 0) {
    if(length(opt$exome) > 1) {
      opt$exome <- opt$exome[1]
      warnings <- c(warnings, paste0("Warning: --exome received more than one argument. Using first value only: ", opt$exome, "."))
    }
    if(!file.exists(opt$exome)) {
      opt$exome_glob <- Sys.glob(paste0(opt$exome, "*"))
      if(length(opt$exome_glob) == 0) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " not found."))
      } else {
        opt$exome <- opt$exome_glob
        warnings <- c(warnings, paste0("Warning: --exome accepted as globbed argument ", opt$exome))
      }
    }
    if(file.exists(opt$exome)) {
      opt$exome_headers <- names(fread(opt$exome, nrows = 0, skip = "Starts"))
      if(!(any(grepl("chrom", opt$exome_headers)))) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected a `chrom` column."))
      }
      if(!(any(grepl("strand", opt$exome_headers)))) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected a `strand` column."))
      }
      if(!(any(grepl("exonStarts", opt$exome_headers)))) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected an `exonStarts` column."))
      }
      if(!(any(grepl("exonEnds", opt$exome_headers)))) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected an `exonEnds` column."))
      }
      if(!(any(grepl("name2", opt$exome_headers)))) {
        errors <- c(errors, paste0("Error: --exome ", opt$exome, " doesn't look like an exome. Expected a `name2` column."))
      }
      if(opt$mode %in% c("cbe", "abe") | grepl("^be", opt$mode)) {
        if(!(any(grepl("name$", opt$exome_headers)))) {
          errors <- c(errors, paste0("Error: --exome ", opt$exome, " is invalid while in ,", opt$mode, " mode. Expected a `name` column."))
        }
        if(!(any(grepl("cdsStart", opt$exome_headers)))) {
          errors <- c(errors, paste0("Error: --exome ", opt$exome, " is invalid while in ,", opt$mode, " mode. Expected a `cdsStart` column."))
        }
        if(!(any(grepl("cdsEnd", opt$exome_headers)))) {
          errors <- c(errors, paste0("Error: --exome ", opt$exome, " is invalid while in ,", opt$mode, " mode. Expected a `cdsEnd` column."))
        }
      }
    }
  # } else if(opt$adhoc) {
  #   errors <- c(errors, paste0("Error: --exome not passed while in ad-hoc mode."))
  } else {
    errors <- c(errors, paste0("Error: --exome not specified."))
  }
  
  # check id
  # if(length(opt$id) > 0) {
  #   if(opt$adhoc) {
  #     if(length(opt$id) > 1) {
  #       opt$id <- opt$id[1]
  #       warnings <- c(warnings, paste0("Warning: --id received more than one argument. Using first value only: ", opt$id, "."))
  #     }
  #     if(!grepl("^\\d+$", opt$id)) {
  #       errors <- c(errors, paste0("Error: --id must be an integer, got ", opt$id, "."))
  #     } else {
  #       opt$id <- as.numeric(opt$id)
  #       if (length(opt$infile) > 0) {
  #         if (opt$id > length(opt$infile_headers)) {
  #           errors <- c(errors, paste0("Error: --id is out of bounds, got ", opt$id, " but --infile only contains ", length(opt$infile_headers), " columns."))
  #         }
  #       }
  #     }
  #   } else {
  #     warnings <- c(warnings, paste0("Warning: --id specified outside of ad-hoc mode. Ignoring."))
  #   }
  # }
  
  # check harm
  if(length(opt$harm) > 0) {
    if(length(opt$harm) > 1) {
      opt$harm <- opt$harm[1]
      warnings <- c(warnings, paste0("Warning: --harm received more than one argument. Using first value only: ", opt$harm, "."))
    }
    if(!grepl("^\\d+$", opt$harm)) {
      errors <- c(errors, paste0("Error: --harm must be an integer, got ", opt$harm, "."))
    } else {
      opt$harm <- as.numeric(opt$harm)
      if (length(opt$infile) > 0) {
        if (opt$harm > length(opt$infile_headers)) {
          errors <- c(errors, paste0("Error: --harm is out of bounds, got ", opt$harm, " but --infile only contains ", length(opt$infile_headers), " columns."))
        }
      }
    }
  }
  
  # check priorities
  if(length(opt$priorities) > 0) {
    warnings <- c(warnings, paste0("Warning: Deprecated option --priorities specified. Ignoring."))
  #   if(length(opt$priorities) > 1) {
  #     opt$priorities <- opt$priorities[1]
  #     warnings <- c(warnings, paste0("Warning: --priorities received more than one argument. Using first value only: ", opt$priorities, "."))
  #   }
  #   if(!file.exists(opt$priorities)) {
  #     opt$priorities_glob <- Sys.glob(paste0(opt$priorities, "*"))
  #     if(length(opt$priorities_glob) == 0) {
  #       errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " not found."))
  #     } else {
  #       opt$priorities <- opt$priorities_glob
  #       warnings <- c(warnings, paste0("Warning: --priorities accepted as globbed argument ", opt$priorities))
  #     }
  #   }
  #   if(file.exists(opt$priorities)) {
  #     opt$priorities_headers <- names(fread(opt$priorities, nrows = 0))
  #     if(!("Symbol" %in% opt$priorities_headers)) {
  #       errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " doesn't look like a feature priorities file. Expected a `Symbol` column."))
  #     }
  #     if(!("Gene Type" %in% opt$priorities_headers)) {
  #       errors <- c(errors, paste0("Error: --priorities ", opt$priorities, " doesn't look like a feature priorities file. Expected a `Gene Type` column."))
  #     }
  #   }
  #   if(length(opt$harm) == 0) {
  #     warnings <- c(warnings, paste0("Warning: --priorities passed without --harm. Ignoring."))
  #   }
  # } else if(length(opt$harm) > 0) {
  #   if(opt$harm != 0) {
  #     errors <- c(errors, paste0("Error: --priorities not passed while --harm passed."))
  #   }
  }
  
  # check control
  if(length(opt$control) > 0) {
    if(length(opt$harm) == 0) {
      opt$control <- NULL
      warnings <- c(warnings, paste0("Warning: --control passed without --harm. Ignoring."))
    }
  }
  
  # check control_type
  if(length(opt$control_type) > 0) {
    if(length(opt$control) == 0) {
      opt$control_type <- NULL
      warnings <- c(warnings, paste0("Warning: --control_types passed without --control. Ignoring."))
    } else if(length(opt$control) == 1) {
      if(length(opt$control_type) > 1) {
        opt$control_type <- opt$control_type[1]
        warnings <- c(warnings, paste0("Warning: --control_type received more than one argument. Using first value only: ", opt$control_type, "."))
      }
    } else if(length(opt$control) > 1) {
      if(length(opt$control_type) == 1) {
        opt$control_type <- rep(opt$control_type, length(opt$control))
      } else if(length(opt$control_type) != length(opt$control)) {
        errors <- c(errors, paste0("Error: --control and --control_type lengths are not equal."))
      }
    }
    if("exo_Non-targeting" %in% opt$control_type) {
      errors <- c(errors, paste0("Error: --control_type contains the disallowed value \"exo_Non-targeting\"."))
    }
    if(any(duplicated(opt$control_type))) {
      errors <- c(errors, paste0("Error: --control_type contains duplicated values."))
    }
  } else if(length(opt$control) > 0) {
    opt$control_type <- paste0("Non-targeting", 1:length(opt$control))
    warnings <- c(warnings, paste0("Warning: --control passed without --control_type. Assuming --control_type for ", opt$control, " is ", opt$control_type, "."))
  }
  
    # check expression
    if(length(opt$expression) > 0) {
      if(length(opt$expression) > 1) {
        opt$expression <- opt$expression[1]
        warnings <- c(warnings, paste0("Warning: --expression received more than one argument. Using first value only: ", opt$expression, "."))
      }
      if(!file.exists(opt$expression)) {
        opt$expression_glob <- Sys.glob(paste0(opt$expression, "*"))
        if(length(opt$expression_glob) == 0) {
          errors <- c(errors, paste0("Error: --expression ", opt$expression, " not found."))
          opt$expression <- NULL
        } else {
          opt$expression <- opt$expression_glob
          warnings <- c(warnings, paste0("Warning: --expression accepted as globbed argument ", opt$expression, "."))
        }
      }
      if(!isFile(opt$expression)) {
        errors <- c(errors, paste0("Error: --expression ", opt$expression, " is not a file."))
        opt$expression <- NULL
      }
    }
    
  # check ref
  if(length(opt$ref) > 0) {
    ref = paste0("
 exorcise version ", ver, " was developed by Dr Simon Lam, University of Cambridge
 https://github.com/SimonLammmm/exorcise

 If you found exorcise useful in your work, please cite:
 * Lam S et al, 2024, Genome-aware annotation of CRISPR guides validates targets in variant cell lines and enhances discovery in screens, Genome Med 16(139), doi: 10.1186/s13073-024-01414-4
 * Kent WJ, 2002, BLAT--the BLAST-like alignment tool, Genome Res 12(4): 656-664, doi: 10.1101/gr.229202

=======================================================================================================")
    cat(ref, "\n\n")
    options(show.error.messages = FALSE)
    stop()
  }
  
  # remove temporary options
  opt$infile_glob <- NULL
  opt$infile_headers <- NULL
  opt$library_glob <- NULL
  opt$library_headers <- NULL
  opt$genome_glob <- NULL
  opt$exome_glob <- NULL
  opt$exome_headers <- NULL
  opt$mode_tolower <- NULL
  opt$priorities_glob <- NULL
  opt$priorities_headers <- NULL
  
  # print warnings and errors
  if(length(warnings) > 0) {
    for (i in 1:length(warnings)) {
      log_warn(warnings[[i]])
    }
  }
  
  if(length(errors) > 0) {
    for (i in 1:length(errors)) {
      log_error(errors[[i]])
    }
    if(!is.null(opt$outdir)) {
      stop("FATAL: There were errors in the input. Quitting. Please check the logfile: ", logfile, ".")
    } else {
      stop("FATAL: There were errors in the input. Quitting.", call. = F)
    }
  }
  
  # if(length(warnings) > 0 & opt$just_go == F) {
  #   readline(prompt = "Accept these settings? Press [enter] to continue.")
  # }
  
  return(opt)
}