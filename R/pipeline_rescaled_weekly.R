library(EpiNow2)
library(data.table)
library(parallel)
library(bayesplot)

.args <- if (interactive()) {
    .prov <- "GP"
    .tmp <- sprintf(
        c(
            "local/data/weekly_%s.rds",
            "local/output/forecast_special_%s.rds"
            ),
        .prov
    )
    c(.tmp[1:length(.tmp) - 1],
      file.path("R", "pipeline_shared_inputs.R"),
      .tmp[length(.tmp)]
    )
} else commandArgs(trailingOnly = TRUE)

# Load helper functions and shared model inputs
source(.args[length(.args) - 1])


####################################
# Parameters
####################################
# when changing units, the mean, sd, and max scale the same way
incubation_period <- LogNormal(mean = 5 / 7, sd = 1 / 7, max = 14 / 7)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7201952/
# generation_time <- LogNormal(mean = 5.2, sd = 1.72, max = 10)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9837419/

# Generation period
generation_time <- Gamma(mean = 7.12 / 7, sd = 1.72 / 7, max = 10 / 7) # mean and sd are more intuitive for rescaling

# Reporting delays
reporting_delay <- LogNormal(mean = 2 / 7, sd = 1 / 7, max = 10 / 7) # mean and sd are more intuitive for rescaling

# Total delays
delay <- incubation_period + reporting_delay

# Rt prior
rt_prior <- LogNormal(meanlog = 0.69, sdlog = 0.05) # mean = 2, sd = 0.1

####################################
# Observation model
####################################
obs <- obs_opts(
  week_effect = FALSE,
  likelihood = TRUE,
  return_likelihood = FALSE
)

###############################
# Pipeline
###############################

# Prepare data
# inflate as.Date, because EpiNow2 seems to prefer Date over IDate
dt <- readRDS(.args[1])[, .(date = as.Date(date), confirm)][!is.na(confirm)]

# EpiNow wants to work in terms of days, so we're going to pretend
# as if weeks are days
dt[, orig_date := date]

fake_daily_dates <- seq.Date(
  from = dt$orig_date[1],
  by = "day",
  length.out = length(dt$orig_date)
)

dt$date <- fake_daily_dates

# Slides for fitting
slides <- seq(0, dt[, .N - (train_window_rescaled + test_window_rescaled)], by = test_window_rescaled)

# Fit models
res_dt <- lapply(slides, \(slide) {
    slice <- dt[seq_len(train_window_rescaled) + slide] |> trim_leading_zero()
    # Slides for fitting are in weeks but we need to rescale back to
    # days for aligning with other scales
    slide_rescaled <- slide * 7
    # Fit model
    if (slice[, .N > test_window_rescaled * 2]) {
        # diagnostics place holder to guarantee entry into while
        diagnostics <- data.table(
            divergent_transitions = 20,
            ess_bulk = 200,
            rhat = 2
        )
        ratchets <- -1
        next_stan <- stan
        stan_elapsed_time <- 0
        crude_run_time <- 0
        # Sources for while loop conditions:
        # - rhat <= 1.05: https://search.r-project.org/CRAN/refmans/rstan/html/Rhat.html AND https://arxiv.org/abs/1903.08008
        # - ess_bulk >= 400: https://search.r-project.org/CRAN/refmans/rstan/html/Rhat.html AND https://arxiv.org/abs/1903.08008
        # - divergences <= 10: all we have is that the divergences should be low, so we're assuming 10 here for now. See
        # https://mc-stan.org/learn-stan/diagnostics-warnings.html#divergent-transitions-after-warmup
        # - To prevent the loop from running forever, we also stop refitting after 15 tries and return/process the last fit.
        while(ratchets < 16 &&
              (diagnostics$divergent_transitions > 10 ||
              diagnostics$ess_bulk < 400 ||
              diagnostics$rhat > 1.05)) {
            # The first ratchet counts as 0
            ratchets <- ratchets + 1
            # fit the model
            out <- epinow(
                data = slice,
                generation_time = generation_time_opts(generation_time),
                delays = delay_opts(delay),
                rt = rt_opts(prior = rt_prior),
                forecast = forecast_opts(horizon = test_window_rescaled),
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
            .(date, sample, value, slide = slide_rescaled)
        ]

        diagnostics <- diagnostics[, slide := slide_rescaled]
        diagnostics <- diagnostics[, stan_elapsed_time := stan_elapsed_time]
        #  NB: NEEDS REVIEW: Currently computes total time taken for warmup and sampling for all chains.

        # Combine the forecast, timing and diagnostics
        res_dt <- data.table(
            forecast = list(forecasts),
            timing = list(
                data.table(
                    slide = slide_rescaled,
                    crude_run_time = crude_run_time,
                    stan_elapsed_time = stan_elapsed_time,
                    keep_run_time = last_run_time,
                    ratchets = ratchets
                )
            ),
            diagnostics = list(diagnostics)
        )
        res_dt
    } else {
        empty_forecast <- data.table(
            date = dt[train_window_rescaled + slide, date + seq_len(test_window_rescaled)],
            sample = NA_integer_, value = NA_integer_, slide = slide_rescaled
        )
        data.table(
            forecast = list(empty_forecast),
            timing = list(data.table(
                slide = slide_rescaled,
                crude_run_time = lubridate::as.duration(NA),
                stan_elapsed_time = lubridate::as.duration(NA),
                keep_run_time = lubridate::as.duration(NA),
                ratchets = NA_integer_
            )),
            diagnostics = list(data.table(
                slide = slide_rescaled,
                "samples" = NA_integer_,
                "max_rhat" = NA_integer_,
                "divergent_transitions" = NA_integer_,
                "per_divergent_transitions" = NA_integer_,
                "max_treedepth" = NA_integer_,
                "no_at_max_treedepth" = NA_integer_,
                "per_at_max_treedepth" = NA_integer_,
                "ess_basic" = NA_integer_,
                "ess_bulk" = NA_integer_,
                "ess_tail" = NA_integer_,
                "rhat" = NA_integer_
            ))
        )
    }
}) |> rbindlist()

# Reach into res_dt and update forecast as follows:
# - replace the fake dates with orig_dates by doing a merge on date
# - Remove "confirm" and "date" which was fake
res_dt[,
       forecast := lapply(forecast, function(x) {
           merge(x, dt)
       })]

res_dt[,
       forecast := lapply(forecast, function(x) {
           x[, date := orig_date
           ][, `:=`(
               orig_date = NULL,
               confirm = NULL
           )]
       })]

# Save output
res_dt |> saveRDS(tail(.args, 1))
