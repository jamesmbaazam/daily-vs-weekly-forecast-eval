library(EpiNow2) # We need to load this even if not needed in the script so the shared_inputs_script doesn't error
library(data.table)

.args <- if (interactive()) {
    .prov <- "GP"
    .tmp <- sprintf(
        c(file.path("local/output",
                    c("forecast_daily_%s.rds",
                      "forecast_weekly_%s.rds",
                      "forecast_rescale_%s.rds"
                    )
        ),
        file.path("local/output", "diagnostics_%s.csv")
        ),
        .prov
    )
    c(.tmp[1:length(.tmp) - 1],
        file.path("./R/pipeline_shared_inputs.R"),
      .tmp[length(.tmp)]
    )
} else {
    commandArgs(trailingOnly = TRUE)
}

# Load helper functions and shared model inputs
source(.args[length(.args) - 1])

# Extract the files
diagnostics_dt_combined <- read_bulk_and_rbind(.args[1:3], "diagnostics")

# Save as csv
write.csv(diagnostics_dt_combined, tail(.args, 1), row.names = FALSE)
