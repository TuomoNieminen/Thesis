# simulation
# ---------

# simulate event ages using the poisson distribution (# events) and discrete uniform distribution (age at event)
# returns a data.frame with
# indv = individual indicator
# event_age = age at event
# this function is used by sim_adverse()
get_event_ages <- function(n_indv, rate_years, start_age, end_age) {
  if(start_age >= end_age) stop("end_age must come after start day")
  
  period_ages <- start_age:end_age
  period_length <- length(period_ages)
  period_rate <- period_length * rate_years / 365.25
  
  # number of events for each individual (sample from poisson distribution)
  n_events <- rpois(n_indv, period_rate)
  
  nonzero_indx <- which(n_events > 0)
  nonzero_events <- n_events[nonzero_indx]
  # individual ID's
  indvs_with_events <- rep(nonzero_indx, times = nonzero_events)
  
  # age at event for each event. for each individual and then for each event, an event age is sampled without replacement
  # from the possible ages (period ages). This will fail is some individual experienced more events than there are days in the period.
  # however the probability for this is almost zero
  event_ages <- do.call(c, lapply(nonzero_events, function(n) { sample(period_ages, n, replace = F)}))
  
  event_data <- data.frame(indv = indvs_with_events, event_age = event_ages)
  
  indv_no_events <- which(n_events==0)
  
  return(list(event_data = event_data, no_events = indv_no_events))
}

# simulate adverse event data
# returns a data.frame with columns:
#
# indv = incdividual indicator
# event_age = age at event
# vac_age = age at vaccination (constant)
# birthday = patient_birthday
# event_date = date of event
# days2event = number of days from vaccination to event
sim_adverse <- function(n_indv = 60000, baseline_rate = 10/1000, RR = 2, risk_period = c(1,14), obs_period_end_age = 2*365, vac_age = 60, cohorts=2014:2015) {
  
  if(identical(cohorts[1],cohorts[2]) | length(cohorts) == 1) {
    cohorts <- cohorts[1]
  } else {
    cohorts <- cohorts[1]:cohorts[2]
  }
  
  risk_rate <- baseline_rate*RR
  risk_start_age <- vac_age + risk_period[1]
  risk_end_age <- vac_age + risk_period[2]
  
  simdata <- list()
  
  for(cohort in cohorts) {
    
    nonrisk.1 <- get_event_ages(n_indv, rate_years = baseline_rate, start_age = 0, end_age = risk_start_age - 1)
    nonrisk.2 <- get_event_ages(n_indv, rate_years = baseline_rate, start_age = risk_end_age + 1, end_age = obs_period_end_age)
    nonrisk_data <- rbind(nonrisk.1$event_data, nonrisk.2$event_data)
    
    risk <- get_event_ages(n_indv, rate_years = risk_rate, start_age = risk_start_age, end_age = risk_end_age)
    
    event_data <- rbind(nonrisk_data, risk$event_data)
    indvs_events <- unique(event_data$indv)
    
    indvs_no_events <- unique(c(nonrisk.1$no_events, nonrisk.2$no_events, risk$no_events))
    indvs_no_events <- indvs_no_events[!indvs_no_events %in% indvs_events]
    
    # add birthdates
    possible_BDs <- seq(as.Date(paste0(cohort,"-01-01")), as.Date(paste0(cohort,"-12-31")), by="days")
    
    BD_events <- sample(possible_BDs, length(indvs_events), replace = T)
    BD_no_events <- sample(possible_BDs, length(indvs_no_events), replace = T)
    
    if(length(indvs_events)==0) {
      event_vac_age <- numeric(0)
    } else {
      event_vac_age <- vac_age
    }
    
    if(length(indvs_no_events)==0) {
      no_event_vac_age <- numeric(0)
    } else {
      no_event_vac_age <- vac_age
    }
    
    birthdays_events <- data.frame(indv = indvs_events, birthday = BD_events, vac_age = event_vac_age)
    birthdays_no_events <- data.frame(indv = indvs_no_events, birthday = BD_no_events, vac_age = no_event_vac_age)
    
    data <- merge(event_data, birthdays_events, by = "indv", all.x = T, all.y = T)
    
    # add event dates
    data$event_date <- data$birthday + data$event_age
    
    # add days from vaccination to event
    data$days2event <- data$event_age - data$vac_age
    
    data <- dplyr::bind_rows(data, birthdays_no_events)
    
    row.names(data) <- 1:nrow(data)
    
    # add year to the individual indentifier
    data$indv <- as.numeric(paste0(data$indv, cohort))
    
    simdata[[cohort]] <- data
    
  }
  
  dplyr::bind_rows(simdata)
}

# utility functions
# -----------------

# Is an integer inside an interval
is.in.interval <- function(intgr, interval) (intgr >= min(interval)) & (intgr <= max(interval))

# Compute the length of an integer interval
range_length <- function(range) {
  interval <- c(range[1], range[length(range)])
  max(interval) - min(interval) + 1
}


# BmaxSPRT
# --------

# add case / control information to data
classify_data <- function(data, case_interval = c(1, 14), control_interval = c(15, 60)) {
  data$case <- 0
  data$control <- 0
  
  data$case[is.in.interval(data$days2event, case_interval)] <- 1
  data$control[is.in.interval(data$days2event, control_interval)] <- 1
  
  obsend <- max(case_interval, control_interval)
  data$obsend <- data$birthday + data$vac_age + obsend
  
  data
}


# Compute the maxSPRT test statistic (vectorized with ifelse)
LLR_Binomial <- function(c_n, n, z) {
  ifelse(n == 0, 0,
         ifelse(c_n == n, c_n*log(z + 1),
                ifelse((z*c_n /(n - c_n)) <= 1, 0,
                       c_n * log(c_n/n) + (n-c_n) * log((n-c_n)/n)-c_n*log(1/(z+1))-(n-c_n)*log(z/(z + 1)))))
  
}

# do retrospective maxSPRT surveillance
run_BmaxSPRT <- function(data, max_events, critical_value, z) {
  
  # exlude non-events and order by the date where the individual had passed thought the observation window
  data <- data[!is.na(data[["event_date"]]),]
  data <- data[order(data[["obsend"]]),]
  
  # compute cumulative number of cases
  data$n_case <- cumsum(data[["case"]])
  
  # compute cumulative number of total events
  data$n <- data[["n_case"]] + cumsum(data[["control"]])
  
  # exclude rows where total number of events exceeds max_events
  data <- data[data[["n"]] <= max_events, ]
  
  if(data[["n"]][length(data[["n"]])] < max_events) warning("End of surveillance not reached")
  
  # compute the B-maxSPRT test statistic for each (n_events, n) pair
  TS <- LLR_Binomial(data[["n_case"]], data[["n"]], z = z)
  
  # for which indice the test statistic first exceeds the critical value
  day <- which(TS > critical_value)[1]
  
  # return the date of detection (NA if signal not generated)
  data[day, "event_date"]
}

# plots
# -----

interval_plot <- function(risk, case, ctrl, ...) {
  x <- 0:150
  h <- 10
  plot(x,x, type = "n", yaxt = "n", ylab = "", xlab = "time", ...)
  # risk window
  abline(v = h + risk[1], col = "red")
  abline(v = h + risk[2], col = "red")
  rect(h + case[1], 30, h + case[2], 40, col = "red", lwd = 0)
  rect(h + ctrl[1] - 1, 30, h + ctrl[2] - 1, 40, col ="grey80", lwd = 0)
  leg <- paste(c("actual risk period", "case interval", "control interval"), lapply(list(risk, case, ctrl), paste, collapse = "-"))
  legend("topright", legend = leg, col = c("red", "red", "grey80"), lty = 1, lwd = c(1, 10,10))
}
