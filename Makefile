
default: endrule

# when testing pipeline, single province target to make
# to override, can provide `ONEPROV=WC` (for example) at end of `make sometarget` invocation
ONEPROV ?= GP

# see example.makefile for notes on how to make this
-include local.makefile

REFDIR ?= local

local:
	mkdir -p $@

# structural definitions
DATDIR := ${REFDIR}/data
FIGDIR := ${REFDIR}/figures
OUTDIR := ${REFDIR}/output

# convenience definitions
# use: $(call R[, optional other arguments])
R = $(strip Rscript $^ $(1) $@)
wget = wget -O $@ $(1)

RENV = .Rprofile

# build renv/library & other renv infrastructure
${RENV}: install.R
	 Rscript --vanilla $^

# for make directory rules
define md
$(1): | ${RENV} ${REFDIR}
	mkdir -p $$@

endef

# define all the necessary directory creation & then `md` them
DIRS := ${DATDIR} ${FIGDIR} ${OUTDIR}

$(foreach dir,${DIRS},$(eval $(call md,${dir})))

# source data
DATAURL := https://raw.githubusercontent.com/dsfsi/covid19za/master/data/covid19za_provincial_cumulative_timeline_confirmed.csv

# get the raw data
${DATDIR}/raw.csv: | ${DATDIR}
	$(call wget,${DATAURL})

# initial organization + saving as binary; no cleaning, only type conversion
# & pivoting to long
${DATDIR}/intermediate.rds: R/import.R ${DATDIR}/raw.csv | ${DATDIR}
	$(call R)

# n.b. raw data also has an UNKNOWN
PROVINCES := GP WC EC KZN FS LP MP NC NW

# define all possible extracts
$(foreach agg,daily weekly,$(foreach tar,${PROVINCES},$(eval EXTRACTS += ${DATDIR}/${agg}_${tar}.rds)))

# extraction rule; also cleans data
${EXTRACTS}: R/extract.R ${DATDIR}/intermediate.rds
	$(call R,$(subst _, ,$(basename $(notdir $@))))

${DATDIR}/daily_RSA.rds: R/aggregate.R $(filter ${DATDIR}/daily_%.rds,${EXTRACTS})
	$(call R)

${DATDIR}/weekly_RSA.rds: R/aggregate.R $(filter ${DATDIR}/weekly_%.rds,${EXTRACTS})
	$(call R)

allextracts: ${EXTRACTS} ${DATDIR}/daily_RSA.rds ${DATDIR}/weekly_RSA.rds

# needs some tweaking, but basically right
${FIGDIR}/incidence.png: R/fig_incidence.R ${DATDIR}/intermediate.rds | ${FIGDIR}
	$(call R)

${FIGDIR}/daily_vs_weekly_%.png: R/fig_daily_vs_weekly.R ${DATDIR}/daily_%.rds ${DATDIR}/weekly_%.rds | ${FIGDIR}
	$(call R)

${FIGDIR}/benchmarks_%.png: R/fig_timing.R ${OUTDIR}/forecast_daily_%.rds ${OUTDIR}/forecast_weekly_%.rds ${OUTDIR}/forecast_rescale_%.rds | ${FIGDIR}
	$(call R)

${FIGDIR}/fig_panel_%.png: \
	R/fig_panel.R \
	${DATDIR}/daily_%.rds \
	${DATDIR}/weekly_%.rds \
	${OUTDIR}/score_%.rds \
	${OUTDIR}/forecast_daily_%.rds \
	${OUTDIR}/forecast_weekly_%.rds \
	${OUTDIR}/forecast_rescale_%.rds \
	${OUTDIR}/diagnostics_%.csv | ${FIGDIR}
	$(call R)

${FIGDIR}/fig_panel_ratchets_%.png: \
	R/fig_panel_ratchets.R \
	${DATDIR}/daily_%.rds \
	${DATDIR}/weekly_%.rds \
	${OUTDIR}/forecast_daily_%.rds \
	${OUTDIR}/forecast_weekly_%.rds \
	${OUTDIR}/forecast_rescale_%.rds | ${FIGDIR}
	$(call R)

${FIGDIR}/score_scatter_%.png: R/fig_crps.R ${OUTDIR}/score_%.rds
	$(call R)

alldvswfigs: $(patsubst %,${FIGDIR}/daily_vs_weekly_%.png,${PROVINCES})

allbenchmarkfigs: $(patsubst %,${FIGDIR}/benchmarks_%.png,${PROVINCES})

allpanelfigs: $(patsubst %,${FIGDIR}/fig_panel_%.png,${PROVINCES})

allratchetfigs: $(patsubst %,${FIGDIR}/fig_panel_ratchets_%.png,${PROVINCES})

allscorescatterfigs: $(patsubst %,${FIGDIR}/score_scatter_%.png,${PROVINCES})

# pattern = some province
DAILYDAT_PAT = ${DATDIR}/daily_%.rds
WEEKLYDAT_PAT = ${DATDIR}/weekly_%.rds
# pattern = province_(daily|weekly|rescale)
FORECAST_PAT = ${OUTDIR}/forecast_%.rds

${FORECAST_PAT}: R/pipeline_main.R ${DATDIR}/%.rds R/pipeline_shared_inputs.R | ${OUTDIR}
	$(call R)

${OUTDIR}/forecast_rescale_%.rds: R/pipeline_rescaled_weekly.R ${DATDIR}/weekly_%.rds R/pipeline_shared_inputs.R | ${OUTDIR}
	$(call R)

${OUTDIR}/score_%.rds: R/score.R ${DATDIR}/daily_%.rds ${DATDIR}/weekly_%.rds ${OUTDIR}/forecast_daily_%.rds ${OUTDIR}/forecast_weekly_%.rds ${OUTDIR}/forecast_rescale_%.rds
	$(call R)

${OUTDIR}/diagnostics_%.csv: R/diagnostics.R ${OUTDIR}/forecast_daily_%.rds ${OUTDIR}/forecast_weekly_%.rds ${OUTDIR}/forecast_rescale_%.rds
	$(call R)

# all targets at once
alldiagnostics: $(patsubst %,${OUTDIR}/diagnostics_%.csv,${PROVINCES} RSA)
allforecasts: $(patsubst %,${OUTDIR}/forecast_daily_%.rds,${PROVINCES} RSA) $(patsubst %,${OUTDIR}/forecast_weekly_%.rds,${PROVINCES} RSA) $(patsubst %,${OUTDIR}/forecast_rescale_%.rds,${PROVINCES} RSA)
allscores: $(patsubst %,${OUTDIR}/score_%.rds,${PROVINCES} RSA)

## Main target
allpanelfigs: $(patsubst %,${FIGDIR}/fig_panel_%.png,${PROVINCES} RSA)

test: $(patsubst %,${OUTDIR}/forecast_daily_%.rds,${ONEPROV}) $(patsubst %,${OUTDIR}/forecast_weekly_%.rds,${ONEPROV}) $(patsubst %,${OUTDIR}/forecast_rescale_%.rds,${ONEPROV})

endrule: allpanelfigs
