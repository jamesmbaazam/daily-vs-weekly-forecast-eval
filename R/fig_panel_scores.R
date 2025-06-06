library(data.table)
library(ggplot2)
library(patchwork)

.args <- if (interactive()) {
    .prov <- "GP"
    .tmp <- sprintf(
        c(
            file.path("local", "data", c("daily_%s.rds", "weekly_%s.rds")), # cases
            file.path("local", "output", "score_%s.rds"), # scores
            file.path("local", "output",
                      c("forecast_daily_%s.rds",
                        "forecast_weekly_%s.rds",
                        "forecast_rescale_%s.rds"
                      )
            ), # forecasts (also contains timing)
            file.path("local", "output", "diagnostics_%s.csv"), # diagnostics
            file.path("local", "figures", "fig_panel_%s.png") # diagnostics
        ),
        .prov)
    c(.tmp[1:length(.tmp) - 1],
      file.path("./R/pipeline_shared_inputs.R"),
      .tmp[length(.tmp)]
    )
} else commandArgs(trailingOnly = TRUE)

# Load helper functions and shared model inputs
source(.args[length(.args) - 1])

# Load the raw data
# Cases
daily_cases <- readRDS(.args[1])
weekly_cases <- readRDS(.args[2])
# Scores
scores <- readRDS(.args[3])

# Forecasts
forecasts <- read_bulk_and_rbind(.args[4:6], "forecast")

# Runtimes
runtimes <- read_bulk_and_rbind(.args[4:6], "timing")

# Diagnostics
diagnostics_dt <- fread(.args[7])

# Other data
# train_window <- 70 # This has to be the same as the one in the pipeline.R script

#####
#Plots
####
# Cases plot
cases_plt <- ggplot() +
	geom_point(data = daily_cases,
		aes(x = date, y = confirm),
		size = 0.5
	) +
	geom_segment(
		aes(x = date - 6.5, xend = date + 0.5, y = confirm/7, yend = confirm/7),
		data = weekly_cases[!is.na(confirm)],
		color = "firebrick"
	) +
	scale_y_log10(
		"Daily Incidence (log10)", sec.axis = sec_axis(
			~ . * 7, name = "Weekly Incidence (log10)"
		)
	) +
    scale_x_date(NULL, date_breaks = "month", date_labels = "%b '%y")

cases_plt

##########
# Scores
##########
# scores plot
score_plt <-
    # First make a layer with all dates present
    ggplot(data = daily_cases) +
    geom_blank(
        aes(x = date, y = max(confirm))
    ) +
    # Now add the scores data
	geom_line(
	    data = scores[, type := forecast],
	    aes(x = date,
	        y = crps,
	        color = type
	    )
	) +
    geom_point(
        data = scores,
        aes(x = date,
            y = crps,
            color = type
        )
    ) +
    scale_x_date(NULL, date_breaks = "month", date_labels = "%b '%y") +
    scale_y_log10() +
    scale_color_brewer(na.translate = FALSE, palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    facet_wrap(~data, ncol = 1, strip.position = "right") +
    labs(y = "CRPS (log10)",
         linetype = "Data",
         color = "Forecast target"
    )

score_plt

### Diagnostics

# Get the slides and their dates for annotating the diagnostics data
slide_dates_dictionary <- forecasts[type == "weekly", .SD[1], by = "slide", .SDcols = c("date")]

# Add the dates by slide
diagnostics_dt <- diagnostics_dt[
    slide_dates_dictionary,
    on = "slide"
]

# Reshape the ESS per sec columns into a single column for plotting
diagnostics_dt_long <- melt(
    diagnostics_dt,
    measure.vars = c("ess_basic", "ess_bulk", "ess_tail"),
    variable.name = "ess_type",
    value.name = "ess_value"
)

# Calculate ESS per sec
diagnostics_dt_long[, ess_per_sec := ess_value/stan_elapsed_time]

# Shorten ess_type values
# diagnostics_dt_long[, ess_type := gsub("ess_(.*)", "\\1", ess_type)]

# Plot
diagnostics_plt <-
    # First make a layer with all dates present
    ggplot(data = daily_cases) +
    geom_blank(
        aes(x = date, y = max(confirm))
    ) +
    # Now add the diagnostics data
    geom_line(
        data = diagnostics_dt_long[ess_type == "ess_tail"], # Only plotting ESS tail
        aes(x = date,
            y = ess_per_sec,
            color = type
        )
    ) +
    geom_point(
        data = diagnostics_dt_long[ess_type == "ess_tail"], # Only plotting ESS tail
        aes(x = date,
            y = ess_per_sec,
            color = type
        ),
        size = 1
    ) +
    scale_y_log10() +
    scale_x_date(NULL, date_breaks = "month", date_labels = "%b '%y") +
    scale_color_brewer(na.translate = FALSE, palette = "Dark2") +
    labs(
        # title = "Effective sample size per second",
        x = "Date",
        y = "ESS tail per sec (log10)",
        color = "Forecast target",
        linetype = "Data"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

diagnostics_plt

# Patchwork
panel_fig <- (cases_plt + score_plt + diagnostics_plt) &
    plot_layout(ncol = 1, guides = "collect", axes = "collect_x") &
    plot_annotation(title = paste(daily_cases$province[1])) &
    theme_minimal() &
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


ggsave(tail(.args, 1), panel_fig, bg = "white", width = 12, height = 6)
