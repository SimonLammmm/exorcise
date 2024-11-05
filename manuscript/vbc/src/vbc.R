library(data.table)
library(dplyr)
library(tidyr)
library(stats)

fread <- function(x, ...) data.table::fread(file = Sys.glob(paste0(x, "*"))[1], ...)

project_root <- getwd()                        # New files will be written here
project_root <- paste0(project_root, "/")
dir.create(paste0(project_root, "VBC-top20-hs"))
dir.create(paste0(project_root, "VBC-top20-hs/author"))
dir.create(paste0(project_root, "VBC-top20-mm"))
dir.create(paste0(project_root, "VBC-top20-mm/author"))

# Load VBC guides
# Download from <https://www.vbc-score.org/download>, untick "top 6 sgRNAs per gene", tick "genome-wide prediction"
vbc_data_hs <- fread("data/vbc_hg38_all_genes.tsv")
vbc_data_mm <- fread("data/vbc_mm10_all_sgRNAs.tsv")

vbc_data_hs <- vbc_data_hs %>% mutate(seqnames = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\1", Position),   # Columnise Position to contexts
                                      context.start = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\2", Position),
                                      context.end = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\3", Position),
                                      strand = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\4", Position))

vbc_data_hs$context.start <- as.character(vbc_data_hs$context.start)
vbc_data_hs$context.end   <- as.character(vbc_data_hs$context.end)

vbc_data_mm <- vbc_data_mm %>% mutate(seqnames = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\1", Position),   # Columnise Position to contexts
                                      context.start = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\2", Position),
                                      context.end = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\3", Position),
                                      strand = sub("^(chr.+?):(\\d+)-(\\d+)\\((.)\\)$", "\\4", Position))

vbc_data_mm$context.start <- as.character(vbc_data_mm$context.start)
vbc_data_mm$context.end   <- as.character(vbc_data_mm$context.end)

# Obtain top 20 guides per gene
vbc_lib_hs <- vbc_data_hs %>%
  filter(`Likely off-targets OT>90` == 0) %>%
  filter(`Very likely off-targets OT>95` == 0) %>%
  group_by(gene) %>%
  arrange(desc(`VBC score`), desc(`Frameshift ratio inDelphi`)) %>%
  mutate(rank = 1:n()) %>%
  filter(rank <= 20) %>%
  ungroup()
fwrite(vbc_lib_hs, paste0(project_root, "VBC-top20-hs/author/VBC-top20-hs.tsv"), sep = "\t")

vbc_lib_mm <- vbc_data_mm %>%
  filter(`Likely off-targets OT>90` == 0) %>%
  filter(`Very likely off-targets OT>95` == 0) %>%
  group_by(gene) %>%
  arrange(desc(`VBC score`), desc(`Frameshift ratio inDelphi`)) %>%
  mutate(rank = 1:n()) %>%
  filter(rank <= 20) %>%
  ungroup()
fwrite(vbc_lib_mm, paste0(project_root, "VBC-top20-mm/author/VBC-top20-mm.tsv"), sep = "\t")

# Run Exorcise
# Download 2bit genome files from <https://hgdownload.soe.ucsc.edu/downloads.html>
# Download exome annotations for the same assembly from <https://genome.ucsc.edu/cgi-bin/hgTables>, choose NCBI RefSeq All, filter chrom doesn't match *_*
# Download symbol priorities from <https://www.ncbi.nlm.nih.gov/datasets/gene/>, select Symbol and Gene Type
system("exorcise -i VBC-top20-hs/author/VBC-top20-hs.tsv -g 2 -n 1 -z NGG -v data/hg38.2020-09-22.2bit -w data/hg38.refseq.exons.tsv.gz -y data/hsa.priorities.tsv.gz -o VBC-top20-hs/")
system("exorcise -i VBC-top20-mm/author/VBC-top20-mm.tsv -g 2 -n 1 -z NGG -v data/mm10.2021-04-08.2bit -w data/mm10.refseq.exons.tsv.gz -y data/mmu.priorities.tsv.gz -o VBC-top20-mm/")

# Load Exorcise results
vbc_exo_hs <- fread(paste0(project_root, "VBC-top20-hs/exorcise.tsv"))
vbc_exo_hs <- vbc_exo_hs %>% select(all_of(grep("^exo_|assembly", names(vbc_exo_hs), value = T)))
vbc_exo_mm <- fread(paste0(project_root, "VBC-top20-mm/exorcise.tsv"))
vbc_exo_mm <- vbc_exo_mm %>% select(all_of(grep("^exo_|assembly", names(vbc_exo_mm), value = T)))
vbc_lib_hs <- vbc_lib_hs %>% select(-seqnames, -strand) %>% left_join(vbc_exo_hs, by = c("sgRNA" = "exo_seq"))
vbc_lib_mm <- vbc_lib_mm %>% select(-seqnames, -strand) %>% left_join(vbc_exo_mm, by = c("sgRNA" = "exo_seq"))

# Remove missed-target effects
vbc_lib_hs <- vbc_lib_hs %>% filter(exo_symbol != "X")
vbc_lib_mm <- vbc_lib_mm %>% filter(exo_symbol != "X")

# Remove other-locus off-target effects
vbc_lib_hs <- vbc_lib_hs %>%
  filter(exo_target != "X") %>%
  filter(grepl("^chr(\\d+|X|Y):", exo_target)) %>%                              # Ignore non-canonical chromosomes
  group_by(exo_id) %>%                                                          # For each guide
  mutate(Exome_hits = length(unique(exo_target))) %>%                           # Count number of exome hits
  mutate(Kind = case_when(Exome_hits == 1 ~ "Same locus",                       # Record if the guide hits the exome at exactly one locus - same locus
                          Exome_hits > 1 ~ "Other locus")) %>%                  # Record if the guide hits the exome at more than one locus
  filter(Kind != "Other locus")

vbc_lib_mm <- vbc_lib_mm %>% 
  filter(exo_target != "X") %>%
  filter(grepl("^chr(\\d+|X|Y):", exo_target)) %>%                              # Ignore non-canonical chromosomes
  group_by(exo_id) %>%                                                          # For each guide
  mutate(Exome_hits = length(unique(exo_target))) %>%                           # Count number of exome hits
  mutate(Kind = case_when(Exome_hits == 1 ~ "Same locus",                       # Record if the guide hits the exome at exactly one locus - same locus
                          Exome_hits > 1 ~ "Other locus")) %>%                  # Record if the guide hits the exome at more than one locus
  filter(Kind != "Other locus")

# Obtain final top 6 guides per gene
vbc_lib_hs <- vbc_lib_hs %>%
  group_by(gene) %>%
  arrange(desc(`VBC score`), desc(`Frameshift ratio inDelphi`)) %>%
  mutate(rank = 1:n()) %>%
  filter(rank <= 6) %>%
  ungroup()
fwrite(vbc_lib_hs, paste0(project_root, "VBC-top6-hs/author/VBC-top6-hs.tsv"), sep = "\t")

vbc_lib_mm <- vbc_lib_mm %>%
  group_by(gene) %>%
  arrange(desc(`VBC score`), desc(`Frameshift ratio inDelphi`)) %>%
  mutate(rank = 1:n()) %>%
  filter(rank <= 6) %>%
  ungroup()
fwrite(vbc_lib_mm, paste0(project_root, "VBC-top6-mm/author/VBC-top6-mm.tsv"), sep = "\t")