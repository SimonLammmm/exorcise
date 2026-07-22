### Process inputs ####

# Import authors' original library files, extracting original symbols, detecting/creating a primary key, and removing adapters where necessary
importAuthorsLib <- function(opt) {
  log_info("Opening file ", opt$infile, "...")
  authors <- fread(opt$infile)
  authors <- authors %>% mutate(exo_id = paste0("exorcise_", 1:dplyr::n()))
  authors <- authors %>%                               # read authors' file
    as_tibble() %>%
    relocate(exo_id, exo_seq = opt$seq, exo_orig = opt$harm) %>%     # select appropriate columns
    mutate(exo_seq = toupper(exo_seq)) %>%
    unique()
  return(authors)
}

# Mask low expression
# Assume the file has symbols in column 1 and expression values in column 2
maskLowExpression <- function(opt, authors) {
  log_info("Using ", opt$expression, " to mask out low expressing genes...")
  expression <- fread(opt$expression) %>% setNames(c("exo_orig", "expression"))
  low_genes <- expression$exo_orig[expression$expression < opt$expression_cutoff]
  authors <- authors %>%
    filter(!(exo_orig %in% low_genes))
  return(authors)
}

# Use authors' original library to extract guide RNA sequences for alignment
extractGuides <- function(opt, authors, blats) {
  log_info("Extracting guide sequences...")
  sgRNAs <- paste0(">", authors$exo_id, "\n", authors$exo_seq, opt$pam)
  dirwrite(list(sgRNAs), blats$file_sgRNAs, sep = "\n", quote = F)
  return()
}