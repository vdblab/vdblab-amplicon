loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)
pooling <- snakemake@params$pooling
err <- readRDS(snakemake@input$error)

# see https://github.com/benjjneb/dada2/issues/1095
# we could get dereplication is done on the fly, but keeping it allows us to have the QC info
sprintf("%s - derepFastq", Sys.time())

manifest <- read.csv(snakemake@input$manifest, sep = "\t", stringsAsFactors = FALSE)
fq_col <- paste0("R", snakemake@wildcards$dir)

derep <- derepFastq(manifest[[fq_col]], verbose = TRUE)
saveRDS(derep, snakemake@output$derep)

sprintf("%s - dada", Sys.time())
dada_out <- dada(derep, err = err, multithread = snakemake@threads, pool = pooling)
saveRDS(dada_out, snakemake@output$dada)
