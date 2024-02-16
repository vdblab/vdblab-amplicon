loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)
library(readr)
getN <- function(x) sum(getUniques(x))

derep_R1 <- readRDS(snakemake@input$derep_R1)
dada_R1 <- readRDS(snakemake@input$dada_R1)

if (snakemake@params$is_paired){
    derep_R2 <- readRDS(snakemake@input$derep_R2)
    dada_R2 <- readRDS(snakemake@input$dada_R2)
    merged <- readRDS(snakemake@input$merged)
    metrics <- data.frame(
        sample_id = snakemake@wildcards$sample,
        derepped_R1 = getN(derep_R1),
        derepped_R2 = getN(derep_R2),
        denoised_R1 = getN(dada_R1),
        denoised_R2 = getN(dada_R2),
        merged = getN(merged)
    )
} else{
    metrics <- data.frame(
        sample_id = snakemake@wildcards$sample,
        derepped_R1 = getN(derep_R1),
        derepped_R2 = NA,
        denoised_R1 = getN(dada_R1),
        denoised_R2 = NA,
        merged = NA
    )
}

write_tsv(metrics, snakemake@output$metrics)
