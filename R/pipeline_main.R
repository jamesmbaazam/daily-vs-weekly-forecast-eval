library(EpiNow2)
library(epiparameter)
library(data.table)
library(parallel)
library(bayesplot)

.args <- if (interactive()) {
	.scale <- "weekly"
	.prov <- "GP"
	tmp <- sprintf(file.path(
		"local", c("data", "output"),
  		c("%s_%s.rds", "forecast_%s_%s.rds")
	), .scale, .prov)
	c(tmp[1], file.path("R", "pipeline_shared_inputs.R"), tmp[2])
} else commandArgs(trailingOnly = TRUE)

# inflate as.Date, because EpiNow2 seems to prefer Date over IDate
dt <- readRDS(.args[1])[, .(date = as.Date(date), confirm)]

source(.args[2])

train_window <- 7*10
test_window <- 7*2

slides <- seq(0, dt[, .N - (train_window + test_window)], by = test_window)

# Incubation period
# Get from epiparameter package doi:10.3390/jcm9020538
sars_cov_incubation_dist <- epiparameter_db(
    disease = "COVID-19",
    single_epiparameter = TRUE,
    epi_name = "incubation period"
)

incubation_params <- c(
	get_parameters(sars_cov_incubation_dist),
	list(max = round(quantile(sars_cov_incubation_dist, 0.999)))
) # Upper 99.9% range needed for EpiNow2

incubation_period <- do.call(LogNormal, incubation_params)

# Generation period - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9837419/
generation_time <- Gamma(mean = 7.12, sd = 1.72, max = 10)

# Reporting delays - corresponds to mean = 2, sd = 1
reporting_delay <- LogNormal(meanlog = 0.58, sdlog = 0.47, max = 10)

# Total delays
delay <- incubation_period + reporting_delay

# Rt prior - mean = 2, sd = 0.1
rt_prior <- LogNormal(meanlog = 0.69, sdlog = 0.05)

# Check if forecasting is being done for daily data or weekly. Will be used to
# turn on/off week effect below; week effect is off for the weekly accumulated data.
is_daily <- sub("//_*", "", basename(.args[1])) == "daily"

# Observation model
obs <- obs_opts(
  week_effect = ifelse(is_daily, TRUE, FALSE), # turn on week effect if data is on daily scale
  likelihood = TRUE,
  return_likelihood = FALSE
)

###############################
# Pipeline
###############################
view_dt <- fill_missing(
    dt, missing_dates = "accumulate", missing_obs = "accumulate"
)

res_dt <- lapply(slides, \(slide) {
	slice <- view_dt[seq_len(train_window) + slide] |> trim_leading_zero()
	if (slice[, .N > (test_window * 2)]) {
		diagnostics <- data.table(divergent_transitions = 20) # place holder to guarantee entry into while
		ratchets <- -1
		next_stan <- stan
		stan_elapsed_time <- 0
		crude_run_time <- 0

		while(diagnostics$divergent_transitions > 10) {

			ratchets <- ratchets + 1

			# Fit the model
			out <- epinow(
				data = slice,
				generation_time = generation_time_opts(generation_time),
				delays = delay_opts(delay),
				rt = rt_opts(prior = rt_prior),
				forecast = forecast_opts(horizon = test_window, accumulate = 1),
				obs = obs,
				stan = next_stan
			)

			# Extract the diagnostic information
			diagnostics <- get_rstan_diagnostics(out$estimates$fit)
			last_run_time <- sum(
				rstan::get_elapsed_time(out$estimates$fit)
			)
			stan_elapsed_time <- stan_elapsed_time + last_run_time
			crude_run_time <- crude_run_time + out$timing
			next_stan <- ratchet_control(next_stan)

		}
		# Extract the forecast cases
		forecasts <- out$estimates$samples[
			variable == "reported_cases" & type == "forecast",
			.(date, sample, value, slide = slide)
			]

		diagnostics <- diagnostics[, slide := slide]
		diagnostics <- diagnostics[, stan_elapsed_time := stan_elapsed_time]
		#  NB: NEEDS REVIEW: Currently computes total time taken for warmup and sampling for all chains.

		# Combine the forecast, timing and diagnostics
		res_dt <- data.table(
		    forecast = list(forecasts),
		    timing = list(
		        data.table(
		            slide = slide,
		            crude_run_time = crude_run_time,
		            stan_elapsed_time = stan_elapsed_time,
					keep_run_time = last_run_time,
					ratchets = ratchets
		        )
		    ),
		    diagnostics = list(diagnostics),
		    fit = ifelse(stan_elapsed_time < lubridate::duration(3), list(out$estimates$fit), list(NA)) # Only save the fit if the runtime is less than specified secs (for memory reasons; the fits are massive = 27 Gb ish)
		)
		res_dt
	} else {
		empty_forecast <- data.table(
			date = dt[train_window + slide, date + seq_len(test_window)],
			sample = NA_integer_, value = NA_integer_, slide = slide
		)
		data.table(
			forecast = list(empty_forecast),
			timing = list(data.table(
			    slide = slide,
			    crude_run_time = lubridate::as.duration(NA),
			    stan_elapsed_time = lubridate::as.duration(NA),
				keep_run_time = lubridate::as.duration(NA),
				ratchets = NA_integer_
			)),
			diagnostics = list(data.table(
				slide = slide,
				"samples" = NA,
				"max_rhat" = NA,
				"divergent_transitions" = NA,
				"per_divergent_transitions" = NA,
				"max_treedepth" = NA,
				"no_at_max_treedepth" = NA,
				"per_at_max_treedepth" = NA,
				"ess_basic" = NA,
				"ess_bulk" = NA,
				"ess_tail" = NA
			)),
			fit = list(NA)
		)
	}
})

res_dt |> rbindlist() |> saveRDS(tail(.args, 1))
