load_quietly <- function(...) {
    suppressMessages(suppressWarnings(library(...)))
}

load_quietly(tidyverse)

tb <- read_csv("itensor-energy-raw.csv", show_col_types = FALSE)
tb <- filter(tb, !HAS_ERRORS) %>% select(-HAS_ERRORS)

extract_runtime <- function(s) {
    t <- str_extract(s, "[0-9.e+]+ ms")
    t <- str_remove(t, " ms")
    as.numeric(t) / 1000
}

extract_energy <- function(s) {
    e <- str_extract(s, "[0-9.]+[A-Z]?$")

    # Multipliers for parsing integers with K/M/G suffixes
    suf <- str_extract(e, "[A-Z]?$")
    mul <- integer(length(suf))
    mul[suf == ""]  <- 1
    mul[suf == "K"] <- 1E3
    mul[suf == "M"] <- 1E6
    mul[suf == "G"] <- 1E9

    e <- str_remove(e, "[A-Z]?$")
    as.numeric(e) * mul
}

tb <- mutate(tb, RUNTIME = extract_runtime(OUTPUT))
tb <- mutate(tb, ENERGY = extract_energy(OUTPUT))

tb <- select(tb, -OUTPUT)

tb <- group_by(tb, PROG, SLURM_NTASKS, QREG_SIZE)

write_csv(tb, "itensor-energy.csv")
