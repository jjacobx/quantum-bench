load_quietly <- function(...) {
    suppressMessages(suppressWarnings(library(...)))
}

load_quietly(glue)
load_quietly(rlang)
load_quietly(tidyverse)

cmd_args = commandArgs(trailingOnly = TRUE)
if (is_empty(cmd_args)) {
    abort("Please provide the name of the job to scrape from. ")
} else {
    job_name <- cmd_args[[1]]
}

DATA_DIR <- file.path("../jobs", job_name, "out")

files_out <- list.files(DATA_DIR, "*.out")
data_table <- tibble()

for (file in files_out) {
    lines <- read_lines(file.path(DATA_DIR, file))

    # Which line seperates metadata and output
    sep <- str_which(lines, "^Running.*:") + 1

    # Match metadata in the form FIELD=value
    meta_lines <- str_extract(lines[1:(sep - 1)], "^[A-Z_]+=.*")
    meta_lines <- na.omit(meta_lines)
    meta_fields <- str_replace(meta_lines, "=.*", "")
    meta_values <- str_replace(meta_lines, ".*=", "")

    data <- as.list(meta_values)
    names(data) <- meta_fields

    out <- str_c(lines[sep:length(lines)], collapse = '\n')
    data <- c(data, OUTPUT = out)

    data_table <- rbind(data_table, as_tibble(data))
}

files_err <- list.files(DATA_DIR, "*.err")
has_errors <- file.size(file.path(DATA_DIR, files_err)) > 0
data_table <- mutate(data_table, HAS_ERRORS = has_errors)

write_csv(data_table, glue(job_name, "-raw.csv"))
