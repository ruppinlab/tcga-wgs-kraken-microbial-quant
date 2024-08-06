suppressPackageStartupMessages({
    library(dplyr)
    library(AnnotationHub)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

ah <- AnnotationHub()
meta <- query(ah, "EuPathDB")
meta_df <- data.frame(TaxID = meta$taxonomyid, Species = meta$species)
meta_df <- meta_df %>% distinct()
meta_df <- meta_df %>%
    arrange(TaxID, desc(Species)) %>%
    filter(!duplicated(TaxID))
meta_df <- meta_df %>% arrange(Species)
write.table(
    meta_df,
    file = snakemake@output[[1]],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
