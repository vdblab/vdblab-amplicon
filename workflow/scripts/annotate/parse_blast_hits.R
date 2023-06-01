loge <- file(snakemake@log[["e"]], open = "wt")
logo <- file(snakemake@log[["o"]], open = "wt")
sink(logo, type = "output")
sink(loge, type = "message")
library(tidyverse)

print("Parsing annotation file")
dt_annotation <- read.csv(snakemake@input[["annotation"]], sep = "\t", header = FALSE, col.names=c("accession_id", "tax_id_species", "taxon"))

# The annotation file will map tax_id_species to a taxonomic annotation.


rownames(dt_annotation) <- dt_annotation$accession_id

blast_names <- c("ASVId", "taxid", "accession", "species", "query_length", "align_length", "nident", "pident", "bitscore", "evalue", "score")
results <- read.csv(snakemake@input[["blast"]], sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = blast_names) %>%
  left_join(dt_annotation, by = c("accession" = "accession_id"))

# results$tax_id_species = dt_annotation[as.character(results$accession),"tax_id_species"]

## following identifies the accessions with the best match, including tie scores
# find the tie score hits
# for the species, replace the first space with underscore, delete second space and beyond, then alphabetical order
blasted <- results %>%
  group_by(ASVId) %>%
  filter(bitscore == max(bitscore)) %>%
  mutate(
    species_short = sub(" ", "_", species),
    species_short = sub(" .*$", "", species_short)
  ) %>%
  arrange(species_short) %>%
  # identify unique names, put them in alphabetical order, and merge them, and also add the % identity
  mutate(
    name = paste(paste(unique(species_short), collapse = ";"), align_length[1], pident[1], sep = ";")
  ) %>%
  slice(1) %>%
  ungroup()



# this next section selects ASVs that should be added back if removed during chimera selection
# see how much of the ASV sequence BLAST was able to align
blasted$length_ratio <- blasted$query_length / blasted$align_length

# examine ASVs with a pident score 97% or better
blasted$pident_97 <- blasted$pident >= 97

# blasted_97<-blasted[blasted$pident_97==TRUE,]
# plot the results
# hist(blasted_97$length_ratio)
# it's a bimodal distribution, so let's use a cut-off of 1.1

print(str(blasted))

blasted$length_ratio_1.1 <- blasted$length_ratio <= 1.1
blasted$passed <- blasted$length_ratio_1.1 == TRUE & blasted$pident_97 == TRUE
print(table(blasted$passed))
# so these are the ASVs that blasted well enough to avoid chimera checking

# The output of blast is send to this file to be load in database.
write.table(blasted %>% filter(passed), snakemake@output[["passed"]], sep = "\t", row.names = FALSE)
write.table(blasted %>% filter(!passed), snakemake@output[["not_passed"]], sep = "\t", row.names = FALSE)
write.table(blasted, snakemake@output[["detailed"]], sep = "\t", row.names = FALSE)



tmp <- blasted %>%
  #   left_join(dt_annotation, by=c("accession"="accession_id")) %>%
  rename(
    asv_temp_id = "ASVId",
    annotation = "taxon",
    blast_pass = "passed"
  ) %>%
  select(
    asv_temp_id,
    blast_pass,
    accession,
    annotation,
    tax_id_species,
    pident,
    evalue,
    score
  )


write.table(tmp, snakemake@output[["annotations"]], sep = "\t", row.names = FALSE)
