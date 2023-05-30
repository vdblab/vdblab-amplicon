loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)
pooling <- snakemake@params$pooling
err <- readRDS(snakemake@input$error)

derep_R1 <- readRDS(snakemake@input$derep_R1)
derep_R2 <- readRDS(snakemake@input$derep_R2)

dada_R1 <- readRDS(snakemake@input$dada_R1)
dada_R2 <- readRDS(snakemake@input$dada_R2)

ff <- 0

merged <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)
saveRDS(merged, snakemake@output$merged)

seqtab <- makeSequenceTable(setNames(list(merged), snakemake@wildcards$sample))
saveRDS(seqtab, snakemake@output$seqtab)


library(dada2)
library(dplyr)
library(readr)
library(tidyr)

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
seqtab_nobimera |>
  as_tibble() |>
  mutate(sample_id = rownames(seqtab_nobimera)) |>
  relocate(sample_id) |>
  write_tsv(snakemake@output$seqtab)

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



sample_metrics <- snakemake@input$sample_metrics |>
  purrr::map(
    ~ readr::read_tsv(
      .,
      col_types = readr::cols(
        .default = "i",
        sample_id = "c"
      )
    )
  ) |>
  dplyr::bind_rows()

no_chim_counts <- readr::read_tsv(snakemake@input$seq_counts) |>
  tidyr::pivot_longer(-sample_id, names_to = "asv", values_to = "count") |>
  dplyr::group_by(sample_id) |>
  dplyr::summarise(no_chimeras = sum(count), .groups = "drop")

sample_metrics |>
  dplyr::left_join(no_chim_counts, by = "sample_id") |>
  readr::write_tsv(snakemake@output$metrics)



getN <- function(x) sum(getUniques(x))

metrics <- data.frame(
  sample_id = snakemake@wildcards$sample,
  derepped_R1 = getN(derep_R1),
  derepped_R2 = getN(derep_R2),
  denoised_R1 = getN(dada_R1),
  denoised_R2 = getN(dada_R2),
  merged = getN(merged)
)
write_tsv(metrics, snakemake@output$metrics)
