# Tuomo.A.Nieminen 2016

# ## Simulation assumptions
#
# - baseline rate: 10 events per thousand person years
# - risk ratio ($RR$) : 2
# - rates are equal for all individuals
# - risk interval length: 14 days (1-14 days after exposure)
# - nonrisk interval length: 46 days (15-60 days after exposure)
# - exposed individuals: 120000 (two birth cohorts)


# SIMULATION
# ----------

# utility functions
# --
rm(list = ls())
setwd("~/TNversio/vaccsafety/Signal/Simulations")
source("simfunctions.R")
# Sequential package functions
source("Sequential/R/Performance.Binomial.R")
source("Sequential/R/CV.Binomial.R")
source("Sequential/R/SampleSize.Binomial.R")


# Risk and nonrisk interval definitions
# --

# Actual risk period
risk_period = c(0, 13)

# case inervals
case_intervals <- list("reference" = c(0,13),
                       "short risk 1" = c(0,11), "short risk 2" = c(0,8),"short risk 3" = c(0,5),
                       "long risk 1" = c(0, 16), "long risk 2" = c(0, 21), "long risk 3" = c(0, 30))
# control intervals
control_intervals <- list("reference" = c(14,41),
                          "short risk 1" = c(12, 12*3 - 1), "short risk 2" = c(9, 9*3 - 1), "short risk 3" = c(6,6*3 - 1),
                          "long risk 1" = c(17, 17*3 - 1), "long risk 2" = c(22, 22*3 - 1), "long risk 3" = c(31, 31*3 - 1))

# names for the interval choices
overlaps <- names(case_intervals)

# the choices as a data frame
cases <- sapply(case_intervals, paste, collapse = "-")
controls <- sapply(control_intervals, paste, collapse = "-")
lengths_case <- sapply(case_intervals, range_length)
lengths_ctrl <- sapply(control_intervals, range_length)
choices <- data.frame("risk interval" = cases, "length risk" = lengths_case, "ctrl interval" =  controls, "length ctrl" = lengths_ctrl)
choices$z <- choices$length.ctrl / choices$length.risk
choices

# Simulation parameters
# ---

# risk ratio
RR <- 1.5

# surveillance stopping parameters for each overlap type
params <- lapply(overlaps, function(overlap) {
  cat("\r computing parameters for ", overlap)
  # choose the case and control intervals and compute period lengths
  case_interval <- case_intervals[[overlap]]
  ctrl_interval <- control_intervals[[overlap]]
  riskdays <- range_length(case_interval)
  nonriskdays <- range_length(ctrl_interval)
  # match ratio
  z <- nonriskdays/riskdays

  # compute the critical value and maximum sample size with type I error 0.05 and type II error 0.1 (power = 0.9)
  N_cv <- SampleSize.Binomial(RR = RR, alpha = 0.05, power = 0.9, z = z)
  max_events <- N_cv$Required_N
  critical_value <- N_cv$cv

  list("cv" = critical_value, "N" = max_events, "z" = z, "case" = case_interval, "ctrl" = ctrl_interval)
})

names(params) <- overlaps

# cohorts and dates of observation
cohorts <- 2010:2015
# baseline rate (per person years)
br <- 20/1000
dates <- seq.Date(from = as.Date(paste0(min(cohorts),"-01-01")),
                  to = Sys.Date(),
                  by = "day")

# number of simulation iterations
N <- 10000
n_indv <- 200000
# Empty matrix for results
results <- data.frame(matrix(0, N, length(overlaps)))
colnames(results) <- overlaps

# simulation main loop
for(i in 1:N) {

  # simulate AE data
  data <- sim_adverse(n_indv = n_indv, baseline_rate = br,  risk_period = risk_period, cohorts = cohorts, RR = RR)

  # for each overlap type run the maxSPRT detection method using the simulated data
  for(overlap in overlaps)  {
    
    cat("\r iter ", i, "overlap ", overlap, "  ")

    # classify the data to cases and controls with the chosen case and control intervals
    data_cl <- classify_data(data, case_interval = params[[overlap]][["case"]], control_interval = params[[overlap]][["ctrl"]])

    # get the detection day of the b-maxSPRT method (NA if not detected)
    detection <- run_BmaxSPRT(data = data_cl, max_events = params[[overlap]][["N"]], critical_value = params[[overlap]][["cv"]], z = params[[overlap]][["z"]])

    # save results
    results[[overlap]][i] <- detection
  }
}

sim <- list(iter = N, case_intervals = case_intervals, control_intervals = control_intervals, dates = dates, results = results, choices = choices,
            risk_period = risk_period)
save(file = "../data/risk_interval_simulation_results.Rda", sim)


# RESULTS
# -------

results <- get(load("risk_interval_simulation_results.Rda"))$results

# detect rates for each interval type
apply(results, 2, function(r) {
  detect_rate <- sum(!is.na(r))/N
  detect_rate
})

# the number of days to detection
start_day <- dates[1]
days2detect <- vapply(results, FUN.VALUE = numeric(N),function(r) {
  -as.numeric(start_day - r)
})

summary(days2detect)

