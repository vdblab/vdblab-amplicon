loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)

manifest <- read.csv(snakemake@input$manifest, sep = "\t", stringsAsFactors = FALSE)
fq_col <- paste0("R", snakemake@wildcards$dir)

sprintf("%s - learnErrors", Sys.time())
err <- learnErrors(
  manifest[[fq_col]],
  multithread = snakemake@threads,
  nbases = 1e8,
  randomize = TRUE
)

sprintf("%s - writing error profile RDS", Sys.time())
saveRDS(err, snakemake@output$error)

sprintf("%s - plotErrors", Sys.time())
p_err <- plotErrors(err, err_in = TRUE, nominalQ = TRUE)

ggplot2::ggsave(p_err, filename = snakemake@output$error_fig, width = 6, height = 6)
