loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)
derep_R1 <- readRDS(snakemake@input$derep_R1)
dada_R1 <- readRDS(snakemake@input$dada_R1)

if (snakemake@params$is_paired){
    derep_R2 <- readRDS(snakemake@input$derep_R2)
    dada_R2 <- readRDS(snakemake@input$dada_R2)

    merged <- mergePairs(dada_R1, derep_R1, dada_R2, derep_R2, verbose = TRUE)

    seqtab <- makeSequenceTable(setNames(list(merged), snakemake@wildcards$sample))
    saveRDS(merged, snakemake@output$merged)
}else{
    seqtab <- makeSequenceTable(setNames(list(dada_R1), snakemake@wildcards$sample))
}
saveRDS(seqtab, snakemake@output$seqtab)
