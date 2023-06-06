#!/usr/bin/env Rscript

# see https://stackoverflow.com/questions/64101921/
loge <- file(snakemake@log[["e"]], open = "wt")
logo <- file(snakemake@log[["o"]], open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)

in_R1 <- snakemake@input[["R1"]]
in_R2 <- snakemake@input[["R2"]]

out_R1 <- snakemake@output[["R1"]]
out_R2 <- snakemake@output[["R2"]]
out_stats <- snakemake@output[["stats"]]
out_fig <- snakemake@output[["figpath"]]

filter_trunclen <- c(
  snakemake@params[["trunclen_R1"]], snakemake@params[["trunclen_R2"]]
)
ncores <- snakemake@threads

packageVersion("dada2")

print(paste0("Running multithreaded steps with ", ncores, " cores"))
sprintf("%s - running FilterAndTrim", Sys.time())

out <- data.frame(filterAndTrim(
  in_R1, out_R1, in_R2, out_R2,
  truncLen = filter_trunclen,
  maxN = 0, maxEE = 2, truncQ = 2, rm.phix = TRUE,
  compress = TRUE, multithread = ncores
))

rownames(out) <- gsub("_noprimers_R1.fastq.gz", "", rownames(out))
print(out)
out$pct_loss = round(100 - (100 * out$reads.out / out$reads.in), 4)
write.table(out , file = out_stats, sep = "\t")

print("saving quality plots")
ggplot2::ggsave(
  plotQualityProfile(c(in_R1, in_R2, out_R1, out_R2)),
  filename = out_fig,
  width = 9,
  height = 5,
  dpi = 200
)

sprintf("%s - done filtering", Sys.time())
