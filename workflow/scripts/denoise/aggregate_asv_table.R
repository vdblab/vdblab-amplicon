#!/usr/bin/env Rscript#
if (!"snakemake" %in% ls()) {
  # for debugging, set the files from the command line
  warning("This is for debug only; usage: Rscript aggregate_asv_table.R  path/to/asvs.fasta,path/to/asvs.fasta path/to/counts.csv,path/to/counts2.csv path/to/trackers.csv,path/to/trackers2.csv")
  scripts_dir <- commandArgs()[grep("--file", commandArgs())]
  source(file.path(dirname(gsub("--file=", "", scripts_dir)), "test_class.R"))
  args <- commandArgs(trailing = TRUE)
  # args <- list("vdb_16S/test_output/dada2/ASV_counts_sample2B..pool1059.csv,vdb_16S/test_output/dada2/ASV_counts_sample2A..pool1059.csv", "vdb_16S/test_output/reports/dada2_stats_sample2A..pool1059.tab,vdb_16S/test_output/reports/dada2_stats_sample2B..pool1059.tab")
  snakemake <- Snakemake(
    input = list("seqtabs" = strsplit(args[1], split = ","), "trackers" = strsplit(args[2], split = ",")),
    output = list(
      "seqtab" = "tmp_seqtab.tab",
      "stats" = "tmp_tracker.tab",
      "fasta" = "tmp_ASVs.fasta",
      "counts" = "tmp_ASVs.counts"
    ),
    params = list("minasvlen" = 200), log = list("tmp.e", "tmp.o"),
    threads = 1
  )
  print(snakemake)
}


loge <- file(snakemake@log[["e"]], open = "wt")
logo <- file(snakemake@log[["o"]], open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

seqtabs <- unlist(snakemake@input[["seqtabs"]])
tracker_files <- unlist(snakemake@input[["trackers"]])




library(tidyverse)
library(dada2)

print("reading sequence tables")
seqtab_obs <- lapply(seqtabs, readRDS)
print(sapply(seqtab_obs , function(x) ncol(x) > 0))
print(sapply(seqtab_obs , function(x) is.null(x)))
# drop empty seqtabs, these could arise from blanks or low-quality samples
seqtab_obs <- seqtab_obs[sapply(seqtab_obs , function(x) ncol(x) > 0)]

if (length(seqtab_obs) > 1){
    print(" merging sequence tables")
    st.all <- mergeSequenceTables(tables=seqtab_obs)
} else if  (length(seqtab_obs) == 1){
    print("only single sequence table passed filtering!")
    st.all <- seqtab_obs
} else{
    stop("No valid sequence tables from this run! Exiting")
}

seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)
write.csv(seqtab, snakemake@output[["seqtab"]])




# write out fasta
ASV_sequences <- colnames(seqtab)
print("number of ASVs")
print(length(ASV_sequences))
sprintf("%s - before writing out ASVs", Sys.time())


get_seq_hash <- function(seq_vector) {
  sapply(tolower(seq_vector), function(x) {
    digest::digest(x, algo = "md5")
  })
}
fasta_strs <- paste(paste0(">", get_seq_hash(ASV_sequences)), ASV_sequences, sep = "\n")
writeLines(fasta_strs, con = snakemake@output[["fasta"]])


# write out counts table
seqtabfinal_counts_raw <- as.data.frame(seqtab)
colnames(seqtabfinal_counts_raw) <- get_seq_hash(ASV_sequences)
seqtabfinal_counts_raw$oligos_id <- rownames(seqtab)

seqtabfinal_counts <- seqtabfinal_counts_raw %>%
  pivot_longer(cols = -oligos_id, names_to = "asv_md5", values_to = "count") %>%
  filter(count > 0)

sprintf("%s - before writing out asv_counts.csv", Sys.time())

write.csv(seqtabfinal_counts, snakemake@output[["counts"]], row.names = F)


# merge and write out tracker for stats
track <- purrr::map(tracker_files, .f = function(x) {
  read.csv(x, sep = "\t")
}) %>% bind_rows()

print(track)
# by=0 means merge on row names
track <- merge(track, data.frame("nonchimera"=rowSums(seqtab)), by=0)
print(track)

write.table(track, file = snakemake@output[["stats"]], row.names = FALSE)
print(str(track))

sprintf("%s - ending script", Sys.time())
