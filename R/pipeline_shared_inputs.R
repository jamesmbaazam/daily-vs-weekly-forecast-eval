
#' Trim leading zeros
#'
#' @param init_dt the raw dt
#'
#' @returns
#' @export
#'
#' @examples
trim_leading_zero <- function (init_dt) {

	first_non_zero <- init_dt[, which.max(confirm != 0)]
	if (first_non_zero == 1) {
		return(init_dt)
	} else {
		# while the first non-NA value is a zero, drop that and all leading values
		while(init_dt[!is.na(confirm)][1, confirm == 0]) {
			init_dt <- init_dt[-(1:which.max(confirm == 0))]
		}

		return(init_dt)
	}
}

#' @title Get rstan diagnostics
#' @description
#' Summarise the diagnostic information contained in a `<stanfit>` object. If
#' the object is not a stanfit object, return a data.table with NA values.
#' This function is adapted from the `{epidist}` R package in
#' https://github.com/epinowcast/epidist/pull/175/files
#'
#' @param fit A stanfit object
#'
#' @return A data.table containing the summarised diagnostics
get_rstan_diagnostics <- function(fit) {
	if (inherits(fit, "stanfit")) {
		np <- bayesplot::nuts_params(fit)
		divergent_indices <- np$Parameter == "divergent__"
		treedepth_indices <- np$Parameter == "treedepth__"
		# Calculating ESS (basic, bulk, and tail)
		# ESS can only be calculated on the extracted variable in the form of a matrix with dimensions iterations x chains
		# Extract the infections variable as that is used for forecasting
		reports_posterior <- posterior::extract_variable_array(
		    posterior::as_draws_array(fit),
		    "reports" # NB: NEEDS REVIEW; is it rather infections??
		)
		# Calculate the different types of ess (basic, bulk, and tail)
		fit_ess_basic <- posterior::ess_basic(reports_posterior)
		fit_ess_bulk <- posterior::ess_bulk(reports_posterior)
		fit_ess_tail <- posterior::ess_tail(reports_posterior)

		diagnostics <- data.table(
			"samples" = nrow(np) / length(unique(np$Parameter)),
			"max_rhat" = round(max(bayesplot::rhat(fit), na.rm = TRUE), 3),
			"divergent_transitions" = sum(np[divergent_indices, ]$Value),
			"per_divergent_transitions" = mean(np[divergent_indices, ]$Value),
			"max_treedepth" = max(np[treedepth_indices, ]$Value),
			"ess_basic" = fit_ess_basic,
			"ess_bulk" = fit_ess_bulk,
			"ess_tail" = fit_ess_tail
		)
		diagnostics[, no_at_max_treedepth :=
									sum(np[treedepth_indices, ]$Value == max_treedepth)
		][, per_at_max_treedepth := no_at_max_treedepth / samples]
	} else{
		diagnostics <- data.table(
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
		)
	}
	return(diagnostics[])
}

control_opts <- list(adapt_delta = 0.8, max_treedepth = 10, stepsize = 0.1)

stan <- stan_opts(
	cores = parallel::detectCores() - 2,
	samples = 5000,
	control = control_opts
)

#' Define new parameter values for tuning the model
#'
#' @param stan_cfg Current stan parameter values
#'
#' @returns New stan parameter values
#' @export
#'
#' @examples
ratchet_control <- function(stan_cfg) within(stan_cfg, {
    control <- within(control, {
        adapt_delta <- adapt_delta + (1 - adapt_delta) * 0.5
        max_treedepth <- max_treedepth + 2
        stepsize <- stepsize * 0.5
    })
})