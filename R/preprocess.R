#### Input handling ####


# Read the user's file and stage it for exorcism: de-duplicate, give every row a
# primary key, and move the sequence column (plus the prior annotation, if given)
# to the front under exorcise names. Columns exorcise does not understand are
# left alone and trail behind in the output.
import_input <- function(opt) {
  log_info("Reading ", opt$infile, "...")
  input <- tibble::as_tibble(read_table(opt$infile))
  if (nrow(input) == 0) {
    abort("--infile ", opt$infile, " has no rows.")
  }

  headers <- names(input)
  rename_map <- c(exo_seq = headers[[opt$seq]])
  if (!is.null(opt$harm)) {
    rename_map <- c(rename_map, c(exo_orig = headers[[opt$harm]]))
  }
  log_info("Taking sequences from column ", opt$seq, " (", headers[[opt$seq]], ")",
           if (!is.null(opt$harm)) {
             paste0(" and prior annotations from column ", opt$harm, " (", headers[[opt$harm]], ")")
           } else "", ".")

  authors <- input %>%
    unique() %>%
    mutate(exo_id = paste0("exorcise_", dplyr::row_number())) %>%
    dplyr::rename(all_of(rename_map)) %>%
    relocate(exo_id, exo_seq, any_of("exo_orig")) %>%
    mutate(exo_seq = toupper(as.character(exo_seq)))

  blank <- !nzchar(authors$exo_seq) | is.na(authors$exo_seq)
  if (any(blank)) {
    log_warn("Dropping ", sum(blank), " rows with an empty sequence.")
    authors <- authors[!blank, ]
  }
  invalid <- grepl("[^ACGTN]", authors$exo_seq)
  if (any(invalid)) {
    log_warn("Dropping ", sum(invalid), " rows whose sequence contains characters ",
             "other than [ACGTN]. Is column ", opt$seq, " really the sequence column?")
    authors <- authors[!invalid, ]
  }
  if (nrow(authors) == 0) {
    abort("No usable sequences left in ", opt$infile, " after filtering.")
  }

  log_info("Staged ", nrow(authors), " unique rows carrying ",
           length(unique(authors$exo_seq)), " distinct sequences.")
  authors
}


# Drop guides whose prior annotation names a gene the sample barely expresses.
# Only meaningful alongside --harm, which validation already enforces.
mask_low_expression <- function(opt, authors) {
  log_info("Masking genes expressed below ", opt$exprcutoff, " according to ",
           opt$expression, "...")

  expression <- read_table(opt$expression, select = 1:2) %>%
    setNames(c("exo_orig", "expression")) %>%
    mutate(expression = suppressWarnings(as.numeric(expression)))

  unparsed <- sum(is.na(expression$expression))
  if (unparsed > 0) {
    log_warn(unparsed, " row(s) in ", opt$expression, " have a non-numeric ",
             "expression value and were ignored. Does the file have a header row?")
  }

  low <- unique(expression$exo_orig[!is.na(expression$expression) &
                                      expression$expression < opt$exprcutoff])
  before <- nrow(authors)
  authors <- authors %>% filter(!(exo_orig %in% low))
  log_info("Expression mask removed ", before - nrow(authors), " of ", before,
           " rows across ", length(low), " low-expressing genes.")

  if (nrow(authors) == 0) {
    abort("The expression mask removed every guide. Is --exprcutoff (",
          opt$exprcutoff, ") too high, or do the symbols in ", opt$expression,
          " not match column ", opt$harm, " of --infile?")
  }
  authors
}


# Write the queries for BLAT. Deliberately not checkpointed: the FASTA is cheap
# to rebuild and must always reflect the current --infile, otherwise a stale
# alignment would go unnoticed.
write_guide_fasta <- function(opt, authors, paths) {
  log_info("Writing ", nrow(authors), " queries to ", paths$file_sgRNAs, "...")
  write_fasta(authors$exo_id, paste0(authors$exo_seq, opt$pam), paths$file_sgRNAs)
}
