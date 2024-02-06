loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)
library(dplyr)
library(readr)
library(tidyr)

derep_R1 <- readRDS(snakemake@input[[1]])

dada_R1 <- readRDS(snakemake@input[[2]])

if (snakemake@params$is_paired){
    derep_R2 <- readRDS(snakemake@input[[3]])
    dada_R2 <- readRDS(snakemake@input[[4]])

    merged <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)
    saveRDS(merged, snakemake@output$merged)
}else{
    merged <- dada_R1
}
seqtab <- makeSequenceTable(merged)
saveRDS(seqtab, snakemake@output$seqtab)


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



sprintf("%s - removeBimeraDenovo", Sys.time())
seqtab_nobimera <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)
sprintf("%s - writing sequence table", Sys.time())
seqtab_nobimera_clean <- seqtab_nobimera |>
  as_tibble() |>
  mutate(sample_id = rownames(seqtab_nobimera)) |>
  relocate(sample_id)
write_tsv(seqtab_nobimera_clean, snakemake@output$seqtab)

asv_sequences <- colnames(seqtab_nobimera)

sprintf("%s - writing ASV fasta", Sys.time())
write_asv_fasta(asv_sequences, snakemake@output$fasta)


seqtabfinal_counts_raw <- as.data.frame(seqtab_nobimera)
colnames(seqtabfinal_counts_raw) <- get_seq_hash(asv_sequences)
seqtabfinal_counts_raw$oligos_id <- rownames(seqtab_nobimera)

seqtabfinal_counts <- seqtabfinal_counts_raw %>%
  pivot_longer(cols = -oligos_id, names_to = "asv_md5", values_to = "count") %>%
  filter(count > 0)

sprintf("%s - before writing out asv_counts.csv", Sys.time())

write.csv(seqtabfinal_counts, snakemake@output$counts, row.names = FALSE)





getN <- function(x) sum(getUniques(x))

if (snakemake@params$is_paired){
    derep_R2 <- readRDS(snakemake@input[[3]])
    dada_R2 <- readRDS(snakemake@input[[4]])
    merged <- readRDS(snakemake@input[[5]])
    metrics <- data.frame(
        derepped_R1 = sapply(derep_R1, getN),
        derepped_R2 = sapply(derep_R2, getN),
        denoised_R1 = sapply(dada_R1, getN),
        denoised_R2 = sapply(dada_R2, getN),
        merged = sapply(merged, getN),
        sample_id = names(dada_R1)
    )
} else{
    metrics <- data.frame(
        derepped_R1 = sapply(derep_R1, getN),
        derepped_R2 = rep(NA, length(derep_R1)),
        denoised_R1 = sapply(dada_R1, getN),
        denoised_R2 = rep(NA, length(derep_R1)),
        merged = rep(NA, length(derep_R1)),
        sample_id = names(dada_R1)
    )
}



rownames(metrics) <- names(derep_R1)
write_tsv(metrics, snakemake@output$metrics)
no_chim_counts <- seqtab_nobimera_clean |>
  tidyr::pivot_longer(-sample_id, names_to = "asv", values_to = "count") |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(no_chimeras = sum(count), .groups = "drop")

print(metrics)
print(no_chim_counts)
metrics |>
  dplyr::left_join(no_chim_counts, by = "sample_id") |>
  readr::write_tsv(snakemake@output$metrics)
