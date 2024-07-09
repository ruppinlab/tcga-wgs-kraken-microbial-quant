suppressPackageStartupMessages({
    library(httr)
})
# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

add_headers("X-Auth-Token" = snakemake@params[["token"]])

GET(
    paste0(
        "https://api.gdc.cancer.gov/slicing/view/",
        snakemake@params[["file_uuid"]],
        "?region=unmapped"
    ),
    write_disk(snakemake@output[[1]], overwrite = TRUE)
)

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
