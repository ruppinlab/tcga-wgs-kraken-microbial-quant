suppressPackageStartupMessages({
    library(GenomicDataCommons)
    library(yaml)
})

config <- yaml.load_file("./config/config.yaml")

file_query <-
    files() %>%
    GenomicDataCommons::filter(
        cases.project.program.name %in% config$input$gdc$program_names &
            cases.samples.sample_type %in% config$input$gdc$sample_types &
            experimental_strategy == config$input$gdc$exp_strategy &
            analysis.workflow_type %in% config$input$gdc$workflow_types
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
        "cases.samples.portions.analytes.aliquots.submitter_id",
        "analysis.metadata.read_groups.read_length"
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
    # read_length = sapply(
    #     file_results$analysis$metadata$read_groups, `[[`, "read_length"
    # ),
    row.names = file_results$file_id,
    stringsAsFactors = FALSE
)

write.table(
    file_meta,
    file = config$input$gdc$meta_file,
    quote = FALSE, sep = "\t", row.names = FALSE
)
