suppressPackageStartupMessages({
    library(GenomicDataCommons)
})

# Sink the stderr and stdout to the snakemake log file
# https://stackoverflow.com/a/48173272
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

file_query <-
    files() %>%
    GenomicDataCommons::filter(
        cases.project.program.name == "TCGA" &
            cases.samples.sample_type %in% snakemake@params[["sample_types"]] &
            experimental_strategy == "WGS" &
            analysis.workflow_type == "BWA with Mark Duplicates and BQSR"
    ) %>%
    GenomicDataCommons::select(c(
        "file_name",
        "cases.project.project_id",
        "cases.case_id",
        "cases.submitter_id",
        "cases.samples.sample_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
        "cases.samples.portions.analytes.aliquots.aliquot_id",
        "cases.samples.portions.analytes.aliquots.submitter_id"
    ))

file_results <- results_all(file_query)

file_meta <- data.frame(
    file_id = file_results$file_id,
    file_name = file_results$file_name,
    project_id = unlist(
        sapply(file_results$cases, `[[`, "project"),
        recursive = FALSE, use.names = FALSE
    ),
    case_id = sapply(file_results$cases, `[[`, "case_id"),
    case_submitter_id = sapply(file_results$cases, `[[`, "submitter_id"),
    sample_id = sapply(
        sapply(file_results$cases, `[[`, "samples"), `[[`, "sample_id"
    ),
    sample_submitter_id = sapply(
        sapply(file_results$cases, `[[`, "samples"), `[[`, "submitter_id"
    ),
    sample_type = sapply(
        sapply(file_results$cases, `[[`, "samples"), `[[`, "sample_type"
    ),
    aliquot_id = sapply(
        sapply(
            sapply(
                sapply(
                    sapply(
                        file_results$cases, `[[`, "samples"
                    ), `[[`, "portions"
                ), `[[`, "analytes"
            ), `[[`, "aliquots"
        ), `[[`, "aliquot_id"
    ),
    aliquot_submitter_id = sapply(
        sapply(
            sapply(
                sapply(
                    sapply(
                        file_results$cases, `[[`, "samples"
                    ), `[[`, "portions"
                ), `[[`, "analytes"
            ), `[[`, "aliquots"
        ), `[[`, "submitter_id"
    ),
    row.names = file_results$file_id,
    stringsAsFactors = FALSE
)

saveRDS(file_meta, file = snakemake@output[[1]])

# Proper syntax to close the connection for the log file but could be optional
# for Snakemake wrapper
sink(type = "message")
sink()
