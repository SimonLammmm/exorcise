#### Base edit ####

# in silico base edit
baseEdit <- function(opt, blats) {
  
  if(file.not.exist.or.zero(blats$file_base_edited)) {
    # Determine base editing window ranges
    mapping <- left_join(fread(blats$file_genomic_ranges_matched), fread(blats$file_genomic_ranges),
                         by = c("seqnames", "start", "end", "strand", "exo_seq", "exo_cut", "assembly", "exo_target"))
    mapping <- mapping %>%
      mutate(be_window_start = case_when(strand == "+" ~ guideBegin + opt$be_window_from,
                                         strand == "-" ~ guideFinal - opt$be_window_to),
             be_window_end = case_when(strand == "+" ~ guideBegin + opt$be_window_to,
                                       strand == "-" ~ guideFinal - opt$be_window_from))
    
    # Determine base editing window content
    baseEdited <- mapping %>% filter(exo_target != "X")
    baseEdited <- GRanges(seqnames = baseEdited$seqnames,
                          ranges = IRanges(start = baseEdited$be_window_start,
                                           end = baseEdited$be_window_end),
                          strand = baseEdited$strand,
                          exo_target = baseEdited$exo_target)
    baseEdited$be_window_seq <- as.character(import.2bit(blats$file_genome, which = baseEdited))
    baseEdited$be_window_seq <- sub("^.", "", baseEdited$be_window_seq) # fix to enforce zero-based half-open ranges
    baseEdited$be_window_seq <- DNAStringSet(baseEdited$be_window_seq)
    
    # In silico base edit all possible edits per guide
    baseEdited <- unique(baseEdited)
    baseEdits <- tibble(
      be_chr = as.character(baseEdited@seqnames),
      exo_target = baseEdited$exo_target,
      strand = as.character(baseEdited@strand),
      start = baseEdited@ranges@start,
      beWindowSeq = as.character(baseEdited$be_window_seq)) %>%
      group_by(exo_target, strand, beWindowSeq) %>%
      reframe(
        be_chr, exo_target, strand, beWindowSeq, start,
        be_position_in_window = case_when(strand == "+" ~ stri_locate_all(beWindowSeq, fixed = opt$be_from[1]),   # find be position relative to window
                                          strand == "-" ~ stri_locate_all(beWindowSeq, fixed = opt$be_from[2])),
        be_position_in_window = be_position_in_window[[1]][,1],
        be_original = substr(beWindowSeq, be_position_in_window, be_position_in_window),
        be_mutation = case_when(be_original == opt$be_from[1] ~ opt$be_to[1],
                                be_original == opt$be_from[2] ~ opt$be_to[2]),
        be_position = be_position_in_window + start - 1
      )
    
    # Calculate bystander mutations
    calculateBystanders <- function(target) {
      bystanders <- tibble()
      # From pairwise to n-wise simultaneous mutations
      for (b in 2:nchar(target$beWindowSeq[1])) {
        if (b <= nrow(target)) {
          combinations = combn(1:nrow(target), b)
          # For each set of n-wise simultaneous mutations, calculate original and mutant sequences
          for (n in 1:ncol(combinations)) {
            comb = combinations[, n]
            comb_position_in_window = c(target$be_position_in_window[min(comb)],
                                        target$be_position_in_window[max(comb)])
            comb_be_original = substr(target$beWindowSeq[1],
                                      comb_position_in_window[1],
                                      comb_position_in_window[2])
            comb_be_mutation = target$beWindowSeq[1]
            for (c in comb) {
              substr(comb_be_mutation, target$be_position_in_window[c], target$be_position_in_window[c]) <- target$be_mutation[c]
            }
            comb_be_mutation = substr(comb_be_mutation,
                                      comb_position_in_window[1],
                                      comb_position_in_window[2])
            combEdits <- tibble(exo_target = target$exo_target[1],
                                strand = target$strand[1],
                                beWindowSeq = target$beWindowSeq[1],
                                be_chr = target$be_chr[1],
                                be_position_in_window = min(comb_position_in_window),
                                be_original = comb_be_original,
                                be_mutation = comb_be_mutation,
                                be_position = be_position_in_window + target$start[1] - 1
            )
            bystanders <- bind_rows(bystanders,
                                    combEdits)
          }
        }
      }
      return(bystanders)
    }
    
    # Calculate bystander mutations  
    bystanders <- foreach(target = unique(baseEdits$exo_target), .combine = "bind_rows") %do% {
      calculateBystanders(target = baseEdits %>% filter(exo_target == target))
    }
    
    # Append to singlet base edits
    baseEdits <- bind_rows(baseEdits, bystanders) %>%
      select(-start)
    
    baseEdits <- baseEdits %>% filter(!is.na(be_position))
    baseEdits2 <- GRanges(seqnames = baseEdits$be_chr,
                          ranges = IRanges(start = baseEdits$be_position,
                                           end = baseEdits$be_position+nchar(baseEdits$be_original)),
                          strand = baseEdits$strand,
                          exo_target = baseEdits$exo_target,
                          be_original = baseEdits$be_original,
                          be_mutation = baseEdits$be_mutation)
    
    
    
    #### In silico splice the entire genome to find CDSs ####
    
    # Obtain all exons in the exome
    exome <- fread(opt$exome, skip = "Starts")
    names(exome) <- sub("^#", "", names(exome))
    colExonStarts <- grep("exonStarts", names(exome), value = T)
    colExonEnds <- grep("exonEnds", names(exome), value = T)
    colExonFrame <- grep("exonFrame", names(exome), value = T)
    colExonChrom <- grep("chrom", names(exome), value = T)
    colExonStrand <- grep("strand", names(exome), value = T)
    colCdsStart <- grep("cdsStart", names(exome), value = T)
    colCdsEnd <- grep("cdsEnd", names(exome), value = T)
    colName <- grep("name$", names(exome), value = T)
    colName2 <- grep("name2", names(exome), value = T)
    inherit <- grep("chr|seqname|exonStarts|exonEnds|str|symbol|name2|name$|cdsStart|cdsEnd|exonFrame|^ranges$|^seqlevels$|^seqlengths$|^isCircular$|^start$|^end$|^width$|^element$", names(as_tibble(exome)), value = T, invert = T)
    exome <- exome %>%
      separate_longer_delim(cols = all_of(c(colExonStarts, colExonEnds, colExonFrame)), delim = ",")
    exome <- exome %>%
      filter(exome[[colExonFrame]] != "")
    exons <- GRanges(seqnames = exome[[colExonChrom]],
                     ranges = IRanges(start = as.numeric(exome[[colExonStarts]]),
                                      end = as.numeric(exome[[colExonEnds]])),
                     strand = exome[[colExonStrand]],
                     cdsStart = exome[[colCdsStart]],
                     cdsEnd = exome[[colCdsEnd]],
                     transcript = exome[[colName]],
                     name2 = exome[[colName2]],
                     exonFrames = exome[[colExonFrame]],
                     exome[inherit])
    
    # Obtain all CDSs in the genome
    cdss <- GRanges(seqnames = exome[[colExonChrom]],
                    ranges = IRanges(start = exome[[colCdsStart]],
                                     end = exome[[colCdsEnd]]),
                    strand = exome[[colExonStrand]],
                    transcript = exome[[colName]],
                    name2 = exome[[colName2]])
    
    # In silico splice
    spliced <- pintersect(exons, cdss)
    
    # Determine transcripts within base editing windows of any guide
    idxExomeTranscriptsAffects <- findOverlaps(baseEdited, exons, ignore.strand = T)
    transcriptsAffected <- unique(exons$transcript[idxExomeTranscriptsAffects@to])
    
    
    # Filter CDS bounds for those affected by guides
    spliced <- spliced %>% filter(transcript %in% transcriptsAffected)
    
    # Obtain CDSs
    spliced$cds <- as.character(import.2bit(opt$genome, which = spliced))
    spliced$cds <- sub("^.", "", spliced$cds)
    spliced$cds <- DNAStringSet(spliced$cds)
    
    # Attach library base edits to spliced exons
    splicedBaseEdits <- findOverlaps(baseEdits2, spliced, ignore.strand = T)
    splicedBaseEdits2 <- tibble(
      exo_target = baseEdits2$exo_target[splicedBaseEdits@from],
      be_original = baseEdits2$be_original[splicedBaseEdits@from],
      be_mutation = baseEdits2$be_mutation[splicedBaseEdits@from],
      be_chr = as.character(baseEdits2@seqnames)[splicedBaseEdits@from],
      be_strand = as.character(baseEdits2@strand)[splicedBaseEdits@from],
      be_position = baseEdits2@ranges@start[splicedBaseEdits@from],
      transcript = spliced$transcript[splicedBaseEdits@to],
      strand = as.character(spliced@strand)[splicedBaseEdits@to],
      name2 = spliced$name2[splicedBaseEdits@to],
      exonFrame = spliced$exonFrames[splicedBaseEdits@to],
      exonStart = spliced@ranges@start[splicedBaseEdits@to],
      exonEnd = spliced@ranges@start[splicedBaseEdits@to] + spliced@ranges@width[splicedBaseEdits@to] - 1,
      cds = as.character(spliced$cds[splicedBaseEdits@to]),
      inherit = as_tibble(spliced)[splicedBaseEdits@to, inherit]
    )
    splicedBaseEdits2 <- splicedBaseEdits2 %>%
      mutate(be_relative_position_in_exon = be_position - exonStart) %>%
      filter(be_relative_position_in_exon >= 0) %>% # remove boundary edits outside of exon (upstream)
      filter(be_relative_position_in_exon < exonEnd - exonStart) # remove boundary edits outside of exon (downstream)
    substr(splicedBaseEdits2$cds, splicedBaseEdits2$be_relative_position_in_exon + 1, splicedBaseEdits2$be_relative_position_in_exon + 1) <- splicedBaseEdits2$be_mutation
    splicedBaseEdits2 <- splicedBaseEdits2 %>%
      filter(cds != "")
    
    
    #### Determine base editing consequences per base edit ####
    
    # Convert exon-relative position to CDS-relative position
    spliced <- spliced %>%
      as_tibble() %>%
      filter(hit == T & exonFrames != -1) %>% # remove ranges not corresponding to CDS
      group_by(transcript) %>%
      arrange(start) %>%
      # Calculate number of positions before each exon's start position
      mutate(widthBefore = c(0, cumsum(width - 1))[1:dplyr::n()]) %>% # enforce zero-based
      mutate(widthBefore = case_when(widthBefore == -1 ~ 0,
                                     T ~ widthBefore)) %>%
      ungroup()
    
    # Convert exon-relative base editing positions to CDS-relative base editing positions
    splicedBaseEdits2 <- splicedBaseEdits2 %>%
      left_join(spliced %>% as_tibble() %>% transmute(exonStart = start, transcript, widthBefore), by = c("exonStart", "transcript")) %>%
      mutate(be_relative_position_in_cds = be_relative_position_in_exon + widthBefore)
    
    # Assemble full-length transcripts without base editing
    spliced2 <- tibble(
      transcript = spliced$transcript,
      name2 = spliced$name2,
      strand = as.character(spliced$strand),
      cds = as.character(spliced$cds),
      exonStart = spliced$start
    ) %>%
      group_by(transcript, name2, strand) %>%
      arrange(exonStart) %>%
      ungroup %>%
      summarise(flseq = paste0(cds, collapse = ""),
                .by = c(transcript, name2, strand))
    
    # In silico base edit
    splicedBaseEdits3 <- splicedBaseEdits2 %>%
      #transmute(exo_target, be_original, be_mutation, be_chr, be_strand, be_position, transcript, strand, name2, be_relative_position_in_cds) %>%
      left_join(spliced2, by = c("transcript", "name2", "strand")) %>%
      filter(be_relative_position_in_cds < nchar(flseq)) # remove base edits 1nt outside of exon bounds
    splicedBaseEdits3$be_flseq <- splicedBaseEdits3$flseq
    for (i in 1:nrow(splicedBaseEdits3)) {
      substr(splicedBaseEdits3$be_flseq[i],
             splicedBaseEdits3$be_relative_position_in_cds[i]+1,
             splicedBaseEdits3$be_relative_position_in_cds[i]+nchar(splicedBaseEdits3$be_original[i])) <- splicedBaseEdits3$be_mutation[i]
    }
    
    # Reverse complement sequences on the minus strand and in silico translate
    # Unedited sequences
    spliced2 <- spliced2 %>%
      mutate(flseq = case_when(strand == "+" ~ flseq,
                               strand == "-" ~ chartr("ATCG", "TAGC", stri_reverse(flseq)))) %>%
      mutate(translated = as.character(translate(DNAStringSet(flseq))))
    
    # Base edited sequences
    splicedBaseEdits3 <- splicedBaseEdits3 %>%
      mutate(flseq = case_when(strand == "+" ~ flseq,
                               strand == "-" ~ chartr("ATCG", "TAGC", stri_reverse(flseq))),
             be_flseq = case_when(strand == "+" ~ be_flseq,
                                  strand == "-" ~ chartr("ATCG", "TAGC", stri_reverse(be_flseq)))) %>%
      # In silico translate
      mutate(be_translated = as.character(translate(DNAStringSet(be_flseq))),
             be_position_in_aa = case_when(strand == "+" ~ floor(be_relative_position_in_cds / 3) + 1,
                                           strand == "-" ~ floor((nchar(be_flseq) - 1 - be_relative_position_in_cds) / 3 + 1)),# AA sequences are ONE-based, NOT ZERO-BASED
             be_modified_in_aa = substr(be_translated, be_position_in_aa, be_position_in_aa + floor((nchar(be_original)-1)/3)))
    
    # Determine if the translated sequence is changed
    # Function to find difference range
    find_difference_range <- function(s1, s2) {
      suppressWarnings({
        # Ensure strings are of equal length
        stopifnot(nchar(s1) == nchar(s2))
        
        # Convert strings to character vectors
        chars1 <- strsplit(s1, "")[[1]]
        chars2 <- strsplit(s2, "")[[1]]
        
        # Find first difference
        first_diff <- which(chars1 != chars2)[1]
        
        # Find last difference
        last_diff <- max(which(chars1 != chars2))
      })
      
      # Return range of differences
      return(list(
        start = first_diff, 
        end = last_diff, 
        length = last_diff - first_diff + 1
      ))
    }
    
    splicedBaseEdits3 <- splicedBaseEdits3 %>%
      rowwise() %>%
      mutate(translated = spliced2$translated[spliced2$transcript == transcript],
             diff = list(find_difference_range(translated, be_translated)),
             be_original_in_aa = case_when(!is.na(diff$start) ~ substr(translated, diff$start, diff$start+diff$length-1),
                                           T ~ substr(translated, be_position_in_aa, be_position_in_aa + floor(nchar(be_original)/3))),
             be_modified_in_aa = case_when(!is.na(diff$start) ~ substr(be_translated, diff$start, diff$start+diff$length-1),
                                           T ~ be_original_in_aa),
             be_position_in_aa = case_when(!is.na(diff$start) ~ diff$start,
                                           T ~ be_position_in_aa)) %>%
      mutate(ntMutationDesc = case_when(strand == "+" ~ paste0(be_original, be_relative_position_in_cds+1, be_mutation),
                                        strand == "-" ~ paste0(chartr("ATCG", "TAGC", be_original), nchar(be_flseq)-be_relative_position_in_cds, chartr("ATCG", "TAGC", be_mutation))),
             aaMutationDesc = paste0(be_original_in_aa, be_position_in_aa, be_modified_in_aa)) %>%
      ungroup()
    
    splicedBaseEdits3 <- splicedBaseEdits3 %>%
      mutate(be_consequence = case_when(be_original_in_aa == be_modified_in_aa ~ "synonymous",
                                        be_modified_in_aa == "*" & be_original != "*" ~ "stopgain",
                                        be_modified_in_aa != "*" & be_original == "*" ~ "stoploss",
                                        be_modified_in_aa != be_original_in_aa ~ "missense"))
    
    #### Summarise final base editing results ####
    
    result <- left_join(
      #x = (mapping %>% transmute(seqnames, start, end, width, strand, exo_symbol, exo_seq, exo_cut, assembly, exome, exo_target) %>% unique()),
      x = mapping %>% unique(),
      y = (splicedBaseEdits3 %>% transmute(exo_target, transcript, transcript_strand = strand,
                                           be_position_in_grch38 = paste0(be_chr, ":", be_position, "-", be_position+nchar(be_original)),
                                           be_nt_mutation = ntMutationDesc, be_aa_mutation = aaMutationDesc, be_consequence)),
      by = "exo_target",
      relationship = "many-to-many"
    ) %>%
      mutate(exo_be = paste0(exo_symbol, ":", be_aa_mutation)) %>% # Add mutation summary column
      unique()
    
    #### To disk
    fwrite(result, blats$file_base_edited, sep = "\t")
  } else {
    log_info("Not recalculating base edits: results already exist, ", blats$file_base_edited, ".")
  }
}
