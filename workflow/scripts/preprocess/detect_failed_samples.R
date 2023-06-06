#!/usr/bin/env Rscript
# TODO: This could be expanded to read in the multiQC aggregation of the fastqc reports

if (!"snakemake" %in% ls()) {
  warning("This is for debug only; usage: Rscript autoexclude.R path/to/trim_report path/to/adapter_report path/to/failures")
  scripts_dir <- commandArgs()[grep("--file", commandArgs())]
  source(file.path(dirname(gsub("--file=", "", scripts_dir)), "test_class.R"))
  args <- commandArgs(trailing = TRUE)
  snakemake <- Snakemake(input = list("trim_report" = args[1], "adapter_contam_report" = args[2]), output = list("out" = args[3]), params = list("min_retained" = .5, "min_read_pairs" = 200, "max_adapter_perc" = 10), log = list("tmp.log"))
  print(snakemake)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

trim_stats <- read.csv(snakemake@input[["trim_report"]], sep = "\t", comment.char = "#")
adapter_stats <- read.csv(snakemake@input[["adapter_contam_report"]], sep = "\t", skip = 3)

samples_high_junk_proportion <- trim_stats[trim_stats$postfilter / trim_stats$prefilter < snakemake@params[["min_retained"]], "sample"]
samples_small_library <- trim_stats[trim_stats$postfilter < snakemake@params[["min_read_pairs"]], "sample"]
samples_adapter_contam <- adapter_stats[adapter_stats$Percentage > snakemake@params[["max_adapter_perc"]], "SampleID"]

# prefix removal reason here for removal here with Autoexclude
print("recommending failing the following samples")
res <- data.frame("sample_id" = character(), "reason_for_failure" = character())
for (i in 1:nrow(trim_stats)) {
  sample <- trim_stats$sample[i]
  tmp_res <- data.frame("sample_id" = sample)
  reason <- ""
  if (sample %in% samples_high_junk_proportion) reason <- paste(sep=";", reason, paste0("Retained less than ", snakemake@params[["min_retained"]], " of reads"))
  if (sample %in% samples_small_library) reason <- paste(sep=";", reason, paste0("Sample contained fewer than  ", snakemake@params[["min_read_pairs"]], " reads"))
  if (sample %in% samples_adapter_contam) reason <- paste(sep=";", reason, paste0("Sample contained over ", snakemake@params[["max_adapter_perc"]], " adapters"))
  if (reason != "") reason <- paste0("Auto-excluder warning:: ", reason)

  tmp_res[, "reason_for_failure"] <- reason
  print(tmp_res)
  if (reason != "") {
    res <- rbind(res, tmp_res)
  }
}
write.table(file = snakemake@output[["out"]], res, quote = FALSE, sep = "\t", row.names = FALSE)
