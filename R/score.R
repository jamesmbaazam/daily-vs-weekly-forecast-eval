
library(data.table)
library(scoringutils)

.args <- if (interactive()) {
	.prov <- "GP"
	.tmp <- c(
	  sprintf(file.path("local", "data", "%s_%s.rds"), c("daily", "weekly"), .prov),
	  sprintf(file.path("local", "output", "forecast_%s_%s.rds"), c("daily", "weekly", "rescale"), .prov),
	  sprintf(file.path("local", "output", "score_%s.rds"), .prov)
	)
	c(.tmp[1:length(.tmp) - 1],
	  file.path("./R/pipeline_shared_inputs.R"),
	  .tmp[length(.tmp)]
	)
} else commandArgs(trailingOnly = TRUE)

# Load helper functions and shared model inputs
source(.args[length(.args) - 1])

# True data
daily_ref_dt <- readRDS(.args[1]) |> setnames("confirm", "true_value")
weekly_ref_dt <- readRDS(.args[2]) |> setnames("confirm", "true_value")
# Forecasts
daily_fore_dt <- readRDS(.args[3])$forecast |> rbindlist() |> setnames("value", "prediction")
weekly_fore_dt <- readRDS(.args[4])$forecast |> rbindlist() |> setnames("value", "prediction")
rescale_fore_dt <- readRDS(.args[5])$forecast |> rbindlist() |> setnames("value", "prediction")

# wherever the true value is NA, we are assuming the prediction should be
# accumulated to wherever the next observation occurs

score_dt <- rbind(
# daily forecast vs daily data
join_and_score(daily_fore_dt, daily_ref_dt) |> dtextract("daily", "daily"),

# daily forecasts vs weekly data
join_and_score(daily_fore_dt, weekly_ref_dt) |> dtextract("daily", "weekly"),

# weekly forecasts vs daily data
join_and_score(weekly_fore_dt, daily_ref_dt) |> dtextract("weekly", "daily"),

# weekly forecasts vs weekly data
join_and_score(weekly_fore_dt, weekly_ref_dt) |> dtextract("weekly", "weekly"),

# weekly-scale forecasts vs weekly data
join_and_score(rescale_fore_dt, weekly_ref_dt) |> dtextract("rescale", "weekly")

)

saveRDS(score_dt, tail(.args, 1))
