#### Output loop ####

# Write master mapping file
exorcisemaster <- function(premaster, opt) {
  
  
  for (i in 1:length(premaster)) {
    premaster[[i]][which(is.na(premaster[[i]]) | premaster[[i]] == "" | grepl("^exo_Non-targeting_", premaster[[i]]))] <- "X"
  }
  
  if(is.null(opt$control)) {
    ncontrols <- length(premaster$exo_symbol[which(premaster$exo_symbol == "X")])
    premaster$exo_symbol[which(premaster$exo_symbol == "X")] <- paste0("exo_Non-targeting", 1:ncontrols)
  }
  
  # Fix controls
  if (!is.null(opt$control)) {
    for (s in 1:length(opt$control)) {
      log_info("Finding control sequences: replacing ", opt$control[s], " with ", opt$control_type[s], ".")
      control_guides <- rep(F, nrow(premaster))
      if ("exo_orig" %in% names(premaster)) control_guides <- control_guides | grepl(opt$control[s], premaster$exo_orig)  # search for the control string in authors' symbols and IDs if those columns exist
      
      nThisControl <- length(premaster$exo_symbol[control_guides & premaster$exo_symbol == "X"])
      premaster$exo_symbol[control_guides & premaster$exo_symbol == "X"] <- paste0(opt$control_type[s], "_", 1:nThisControl)   # Map control guides to control annotations unless there is an approved Symbol column
    }
    
    ncontrols <- length(premaster$exo_symbol[which(premaster$exo_symbol == "X")])
    premaster$exo_symbol[which(premaster$exo_symbol == "X")] <- paste0("exo_Non-targeting_", 1:ncontrols) # Catch the remaining non-targeting guides
  }
  
  master <- premaster %>% unique()
  
  # Harmonise
  
  
  # if ("exo_orig" %in% names(master)) {
  #   log_info("Harmonising...")
  #   simple <- exorciseinferTargets(master, opt)
  #   master <- left_join(master, simple, by = "exo_orig") %>%
  #     mutate(exo_harm = case_when(is.na(exo_harm) ~ exo_symbol, # Accept harmonised symbol if found
  #                                 T ~ exo_harm))                # Accept exorcised symbol if no harmonised symbol
  #   
  #   if (length(opt$control) > 0) { # For control guides, retain control types in the harmonised column, ie. without exorcism
  #     for (i in 1:length(opt$control)) {
  #       thisControlPattern <- opt$control[i]
  #       thisControlType <- opt$control_type[i]
  #       master <- master %>%
  #         mutate(exo_harm = case_when(grepl(thisControlPattern, exo_orig) ~ paste0(thisControlType, "_", 1:dplyr::n()),
  #                                     T ~ exo_harm))
  #     }
  #   }
  #   
  #   # Create a second ID column uniquely identifying harmonised sequences
  #   master <- master %>% mutate(exo_id_harm = paste0("exorcise_", exo_seq, "_", exo_harm))
  #   
  #   if(opt$mode %in% c("cbe", "abe") | grepl("^be", opt$mode)) {
  #     master <- master %>%
  #       mutate(exo_id_harm = paste0(exo_id_harm, ":", be_aa_mutation),
  #              exo_be_harm = paste0(exo_harm, ":", be_aa_mutation)) # Add harmonised mutation summary column if using BE mode
  #   }
  #   
  #   master <- master %>% relocate(exo_id, exo_id_harm, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut)
  # } else {
  #   
  #   # If not harmonise
    master <- master %>% relocate(exo_id, exo_seq, exo_symbol, exo_target, exo_cut)
  # }
  
  # Finalise
  master <- master %>% mutate(exo_id = paste0("exorcise_", exo_seq, "_", exo_target, "_", exo_symbol))
  
  if(opt$mode %in% c("cbe", "abe") | grepl("^be", opt$mode)) {
    master <- master %>% mutate(exo_id = paste0(exo_id, ":", be_aa_mutation)) # Add mutation if using BE mode
  }
  
  # Fix where values are "X"
  master <- lapply(master, function(x) { x[x == "X"] <- NA; x } )
  
  return(master)
}

# # Infer intended target per authors' original symbols and write
# exorciseinferTargets <- function(master, opt) {
#   
#   simple <- master %>%
#     dplyr::select(exo_symbol, exo_orig, exo_seq, matches("^inherit")) %>% unique() %>% dplyr::select(-exo_seq)
#   
#   # Load in gene symbol annotations for ranking. This file contains all possible approved symbols and gene classes
#   annot <- fread(opt$priorities)
#   names(annot) <- c("Symbol", "Gene Type")
#   annot$`Gene Type` <- factor(annot$`Gene Type`, ordered = T, levels = c("PROTEIN_CODING", "ncRNA", "PSEUDO", "rRNA", "tRNA", "snRNA", "snoRNA", "scRNA", "BIOLOGICAL_REGION", "OTHER", "")) 
#   
#   # Infer
#   
#   simple_dup <- simple %>%
#     left_join(annot, by = c("exo_symbol" = "Symbol"), relationship = "many-to-many")                  # Join the multi-mapped original symbols to annotations for their mapping target
#   simple_dup$`Gene Type`[which(is.na(simple_dup$`Gene Type`))] <- "OTHER"
#   
#   simple_dup2 <- simple_dup %>% filter(!grepl(paste0(paste0("(^", c("X", opt$control_type), "\\d+$)"), collapse = "|"), exo_symbol)) # Remove controls and unmapped guides when considering intended target
#   
#   # Harmonise
#   simple <- simple_dup2 %>%
#     group_by(exo_symbol, exo_orig, `Gene Type`) %>% # For each author's symbol
#     reframe(n = dplyr::n()) %>%
#     group_by(exo_orig) %>%
#     filter(n == max(n)) %>% # Count candidate frequency and accept the most frequent
#     filter(as.numeric(`Gene Type`) == min(as.numeric(`Gene Type`))) %>% # Of those remaining, accept the highest priority gene
#     mutate(m = exo_symbol == exo_orig) %>%
#     filter(m == max(m)) %>% # Of those remaining, accept matching gene symbol with the authors', if any
#     arrange(exo_symbol) %>%
#     mutate(t = 1:dplyr::n()) %>%
#     filter(t == min(t)) %>% # Of those remaining, accept the first lexicographical symbol
#     dplyr::select(exo_harm = exo_symbol, exo_orig) %>%
#     ungroup()
#   
#   if(any(grepl("^inherit", names(simple_dup2)))) {
#     simple <- simple %>%
#       left_join(simple_dup2 %>% select(exo_symbol, matches("^inherit")), by = c("exo_harm" = "exo_symbol"), relationship = "many-to-many") %>% unique()
#     names(simple)[grep("^inherit", names(simple))] <- paste0("exo_harm.", names(simple)[grep("^inherit", names(simple))])
#   }
#   
#   return (simple)
#   
# }

# reannotateExisting <- function(opt) {
#   
#   input <- fread(opt$infile) %>%
#     relocate(exo_seq = opt$seq) %>%
#     unique()
#   
#   library <- fread(opt$library)
#   
#   if(!is.null(opt$harm)) {
#     input <- input %>% relocate(exo_orig = opt$harm, exo_seq)
#     library <- library %>% select(exo_id, exo_seq, exo_symbol, exo_target, exo_cut) %>% unique()                         # harm specified: harmonise
#     exorcised <- left_join(input, library, by = "exo_seq")
#     exorcised <- exorcisemaster(exorcised, opt)
#   } else {
#     if("exo_harm" %in% names(library)) {
#       log_info("Inheriting existing harmonisations from the library...")
#       library <- library %>% select(exo_id, exo_id_harm, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut) %>% unique()   # harm not specified: inherit harmonisations if exist
#     } else {
#       library <- library %>% select(exo_id, exo_seq, exo_symbol, exo_target, exo_cut) %>% unique()
#     }
#     exorcised <- left_join(input, library, by = "exo_seq")
#   }
#   
#   if("exo_harm" %in% names(exorcised)) {
#     exorcised <- exorcised %>% relocate(exo_id, exo_id_harm, exo_seq, exo_symbol, exo_harm, exo_orig, exo_target, exo_cut)
#   } else {
#     exorcised <- exorcised %>% relocate(exo_id, exo_seq, exo_symbol, exo_target, exo_cut)
#   }
#   
#   fwrite(exorcised, paste0(opt$outdir, "/exorcise.tsv"), sep = "\t")
#   
# }