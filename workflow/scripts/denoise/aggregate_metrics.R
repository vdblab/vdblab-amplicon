loge <- file(snakemake@log$e, open = "wt")
logo <- file(snakemake@log$o, open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")

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
