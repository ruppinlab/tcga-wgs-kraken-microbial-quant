suppressPackageStartupMessages({
    library(dplyr)
    library(GenomicDataCommons)
    library(yaml)
})

config <- yaml.load_file("./config/config.yaml")

stopifnot(GenomicDataCommons::status()$status == "OK")

file_query <-
    files() %>%
    GenomicDataCommons::filter(
        cases.project.program.name %in% config$gdc$program_names &
            cases.samples.sample_type %in% config$gdc$sample_types &
            experimental_strategy == config$gdc$exp_strategy &
            analysis.workflow_type %in% config$gdc$workflow_types
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
    )) %>%
    GenomicDataCommons::expand("analysis.metadata.read_groups")

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

readgrp_meta <- merge(
    file_meta,
    as.data.frame(
        bind_rows(
            setNames(
                file_results$analysis$metadata$read_groups, file_results$file_id
            ),
            .id = "file_id"
        ),
    ),
    by = "file_id"
)

uniq_readgrps <- dplyr::distinct(
    readgrp_meta, file_id, read_length, is_paired_end
)
num_uniq_readgrps <- dplyr::count(
    uniq_readgrps, file_id,
    name = "num_uniq_read_groups"
)
single_readgrps <- semi_join(
    uniq_readgrps, dplyr::filter(num_uniq_readgrps, num_uniq_read_groups == 1),
    by = join_by(file_id)
)

file_meta <- merge(
    file_meta, num_uniq_readgrps,
    by = "file_id", all.x = TRUE
)
file_meta <- merge(
    file_meta, single_readgrps,
    by = "file_id", all.x = TRUE
)
row.names(file_meta) <- file_meta$file_id
file_meta <- arrange(file_meta, file_id)

readgrp_meta <- merge(readgrp_meta, num_uniq_readgrps, by = "file_id")
row.names(readgrp_meta) <- readgrp_meta$read_group_id
readgrp_meta <- arrange(readgrp_meta, read_group_id)

data_dir <- config$gdc$metadata$data_dir
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE, mode = "0755")
file_meta_filename <- config$gdc$metadata$file_meta_filename
cat("Writing", file_meta_filename, "\n")
write.table(
    file_meta,
    file = paste(data_dir, file_meta_filename, sep = "/"),
    quote = FALSE, sep = "\t", row.names = FALSE
)
readgrp_meta_filename <-
    config$gdc$metadata$readgrp_meta_filename
cat("Writing", readgrp_meta_filename, "\n")
write.table(
    readgrp_meta,
    file = paste(data_dir, readgrp_meta_filename, sep = "/"),
    quote = FALSE, sep = "\t", row.names = FALSE
)
uniq_readgrps_filename <-
    config$gdc$metadata$uniq_readgrps_filename
cat("Writing", uniq_readgrps_filename, "\n")
write.table(
    uniq_readgrps,
    file = paste(data_dir, uniq_readgrps_filename, sep = "/"),
    quote = FALSE, sep = "\t", row.names = FALSE
)
