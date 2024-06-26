#!/usr/bin/env Rscript

# see https://stackoverflow.com/questions/64101921/
loge <- file(snakemake@log[["e"]], open = "wt")
logo <- file(snakemake@log[["o"]], open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

library(dada2)

inreads <- snakemake@input[["reads"]]

out_R1 <- snakemake@output[["R1"]]
out_R2 <- snakemake@output[["R2"]]
out_stats <- snakemake@output[["stats"]]
out_fig <- snakemake@output[["figpath"]]
is_paired <- snakemake@params[["is_paired"]]
filter_trunclen <- c(
  snakemake@params[["trunclen_R1"]], snakemake@params[["trunclen_R2"]]
)
maxLen = ifelse(snakemake@params[["maxLen"]] == 0, Inf, snakemake@params[["maxLen"]])
print(maxLen)
ncores <- snakemake@threads

packageVersion("dada2")

print(paste0("Running multithreaded steps with ", ncores, " cores"))
sprintf("%s - running FilterAndTrim", Sys.time())

if (is_paired){
    out <- data.frame(filterAndTrim(
        inreads[1], out_R1, inreads[2], out_R2,
        truncLen = filter_trunclen,
        trimLeft = snakemake@params[["trimLeft"]],
        trimRight = snakemake@params[["trimRight"]],
        maxLen = maxLen,
        maxN = snakemake@params[["maxN"]],
        maxEE = snakemake@params[["maxEE"]],
        truncQ = snakemake@params[["truncQ"]],
        rm.phix = snakemake@params[["rmphix"]],
        rm.lowcomplex = snakemake@params[["rmlowcomplex"]],
        compress = TRUE,
        multithread = ncores
    ))
} else{
    out <- data.frame(filterAndTrim(
        inreads[1], out_R1,
        truncLen = filter_trunclen[1],
        trimLeft = snakemake@params[["trimLeft"]],
        trimRight = snakemake@params[["trimRight"]],
        maxLen = maxLen,
        maxN = snakemake@params[["maxN"]],
        maxEE = snakemake@params[["maxEE"]],
        truncQ = snakemake@params[["truncQ"]],
        rm.phix = snakemake@params[["rmphix"]],
        rm.lowcomplex = snakemake@params[["rmlowcomplex"]],
        compress = TRUE,
        multithread = ncores
    ))
}

rownames(out) <- gsub("_noprimers_R1.fastq.gz", "", rownames(out))
print(out)
out$pct_loss = round(100 - (100 * out$reads.out / out$reads.in), 4)
write.table(out , file = out_stats, sep = "\t")

print("saving quality plots")
toplot <- c(inreads[1], out_R1)
if (is_paired){
    toplot <- c(inreads[1], inreads[2], out_R1, out_R2)
}
ggplot2::ggsave(
  plotQualityProfile(toplot),
  filename = out_fig,
  width = 4.5 * (length(toplot)/2),
  height = 5,
  dpi = 200
)

sprintf("%s - done filtering", Sys.time())
