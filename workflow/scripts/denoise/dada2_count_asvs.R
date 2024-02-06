loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)

derep_R1 <- readRDS(snakemake@input[[1]])

dada_R1 <- readRDS(snakemake@input$dada_R1[[2]])

if (snakemake@params$is_paired){
    derep_R2 <- readRDS(snakemake@input[[3]])
    dada_R2 <- readRDS(snakemake@input[[4]])

    merged <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)

    seqtab <- makeSequenceTable(setNames(list(merged), snakemake@wildcards$sample))
}else{
    seqtab <- makeSequenceTable(setNames(list(dada), snakemake@wildcards$sample))
}
saveRDS(merged, snakemake@output$merged)
saveRDS(seqtab, snakemake@output$seqtab)
