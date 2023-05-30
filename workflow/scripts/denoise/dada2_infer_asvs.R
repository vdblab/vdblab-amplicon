loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)

err <- readRDS(snakemake@input$error)

# see https://github.com/benjjneb/dada2/issues/1095
# we could get dereplication is done on the fly, but keeping it allows us to have the QC info
sprintf("%s - derepFastq", Sys.time())

derep <- derepFastq(snakemake@input$fastq, verbose = TRUE)
saveRDS(derep, snakemake@output$derep)

sprintf("%s - dada", Sys.time())
dada_out <- dada(derep, err = err, multithread = snakemake@threads, pool = "none")
saveRDS(dada_out, snakemake@output$dada)
