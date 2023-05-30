loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)
library(dplyr)
library(readr)
library(tidyr)

load_sequence_counts <- function(seq_count_paths) {
  sample_seqtabs <- lapply(snakemake@input$seqtab, readRDS)

  # drop empty seq_counts, these could arise from blanks or low-quality samples
  sample_seqtabs <- sample_seqtabs[
    sapply(sample_seqtabs , function(x) ncol(x) > 0)
  ]

  if (length(sample_seqtabs) > 1) {
    seqtab <- mergeSequenceTables(tables = sample_seqtabs)
  } else if (length(sample_seqtabs) == 1) {
    seqtab <- sample_seqtabs
  } else {
    seqtab <- NULL
  }
  seqtab
}

get_seq_hash <- function(seq_vector) {
  sapply(tolower(seq_vector), function(x) {
    digest::digest(x, algo = "md5")
  })
}

write_asv_fasta <- function(asv_sequences, fasta_path) {
  fasta_strs <- paste(
    paste0(">", get_seq_hash(asv_sequences)), asv_sequences, sep = "\n"
  )
  writeLines(fasta_strs, con = fasta_path)
}

sprintf("%s - reading sequence tables", Sys.time())
seqtab <- load_sequence_counts(snakemake@input$seqtab)

if (is.null(seqtab)) stop("No valid sequence tables from this run! Exiting")

sprintf("%s - removeBimeraDenovo", Sys.time())
seqtab <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

sprintf("%s - writing sequence table", Sys.time())
seqtab |>
  as_tibble() |>
  mutate(sample_id = rownames(seqtab)) |>
  relocate(sample_id) |>
  write_tsv(snakemake@output$seqtab)

asv_sequences <- colnames(seqtab)

sprintf("%s - writing ASV fasta", Sys.time())
write_asv_fasta(asv_sequences, snakemake@output$fasta)


seqtabfinal_counts_raw <- as.data.frame(seqtab)
colnames(seqtabfinal_counts_raw) <- get_seq_hash(asv_sequences)
seqtabfinal_counts_raw$oligos_id <- rownames(seqtab)

seqtabfinal_counts <- seqtabfinal_counts_raw %>%
  pivot_longer(cols = -oligos_id, names_to = "asv_md5", values_to = "count") %>%
  filter(count > 0)

sprintf("%s - before writing out asv_counts.csv", Sys.time())

write.csv(seqtabfinal_counts, snakemake@output$counts, row.names = FALSE)
