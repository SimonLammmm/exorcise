#### Assembling the output ####


# Reserved name for guides exorcise could not annotate. --control_type may not
# use it, so a reader can always tell a user-declared control from a guide that
# simply failed to map into an exon.
NON_TARGETING_PREFIX <- "exo_Non-targeting"

# Columns exorcise itself produces, in the order they should appear. Anything
# else the user brought with them trails behind, untouched.
EXORCISE_COLUMN_ORDER <- c(
  "exo_id", "exo_seq", "exo_symbol", "exo_orig", "exo_target", "exo_cut", "exo_be"
)


# Guides with no exonic hit get a name rather than an empty cell, so that
# downstream counting tools do not silently drop them. A guide whose prior
# annotation matches a --control pattern is named after that control; the rest
# are numbered as non-targeting.
label_unannotated <- function(master, opt) {
  if (!"exo_symbol" %in% names(master)) {
    master$exo_symbol <- NA_character_
  }
  master$exo_symbol <- as.character(master$exo_symbol)

  unannotated <- is.na(master$exo_symbol) |
    !nzchar(master$exo_symbol) |
    startsWith(master$exo_symbol, paste0(NON_TARGETING_PREFIX, "_"))

  if (length(opt$control) > 0) {
    prior <- if ("exo_orig" %in% names(master)) {
      as.character(master$exo_orig)
    } else {
      rep(NA_character_, nrow(master))
    }
    for (i in seq_along(opt$control)) {
      matched <- unannotated & !is.na(prior) & grepl(opt$control[[i]], prior)
      if (any(matched)) {
        log_info("Naming ", sum(matched), " unannotated guides matching /",
                 opt$control[[i]], "/ as ", opt$control_type[[i]], ".")
        master$exo_symbol[matched] <- paste0(opt$control_type[[i]], "_",
                                             seq_len(sum(matched)))
      } else {
        log_warn("No unannotated guide matched --control pattern /",
                 opt$control[[i]], "/.")
      }
      unannotated <- unannotated & !matched
    }
  }

  remaining <- sum(unannotated)
  if (remaining > 0) {
    log_info("Naming ", remaining, " remaining unannotated guides as ",
             NON_TARGETING_PREFIX, ".")
    master$exo_symbol[unannotated] <- paste0(NON_TARGETING_PREFIX, "_",
                                             seq_len(remaining))
  }
  master
}


finalise_output <- function(master, opt) {
  master <- tibble::as_tibble(master) %>% label_unannotated(opt)

  # exo_id identifies a reannotation, not an input row: the same sequence hitting
  # two loci is two reannotations, and in base editor mode each predicted protein
  # change is its own reannotation too.
  master <- master %>%
    mutate(exo_id = paste0("exorcise_", exo_seq, "_", exo_target, "_", exo_symbol))
  if (is_base_edit_mode(opt$mode)) {
    master <- master %>% mutate(exo_id = paste0(exo_id, ":", be_aa_mutation))
  }

  master <- master %>%
    relocate(any_of(EXORCISE_COLUMN_ORDER)) %>%
    unique()

  log_info("Assembled ", nrow(master), " reannotations covering ",
           length(unique(master$exo_seq)), " sequences and ",
           length(unique(master$exo_symbol)), " symbols.")
  unannotated <- sum(startsWith(master$exo_symbol, NON_TARGETING_PREFIX))
  if (unannotated > 0) {
    log_info(unannotated, " of those reannotations are non-targeting.")
  }
  master
}
