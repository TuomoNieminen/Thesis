# functions to evaluate the maxSPRT test statistic for Binomial and Poisson data

# auth. Tuomo Nieminen

# data formatting
# ---------------

classify_data <- function(data, case_interval = c(1, 14), control_interval = c(15, 60)) {
  data$case <- 0
  data$control <- 0
  
  data$case[is.in.interval(data$time2hosp, case_interval)] <- 1
  data$control[is.in.interval(data$time2hosp, control_interval)] <- 1
  
  data
}


# maxSPRT test statistic, Poisson likelihood
# mu_t : expected events
# c_t : observed events
LLR.Poisson <- function(mu_t, c_t) {
  if( c_t < mu_t | mu_t == 0 ) {
    return(0)
  } else {
    (mu_t - c_t) + c_t * log(c_t / mu_t)
  }
}

# maxSPRT test statistic, Binomial likelihood
# c_n : observed events in risk group
# n : total observed events
# z : expected proportion of controls to cases
LLR.Binomial <- function(c_n, n, z) {
  
  if(n==0) return(0)
  
  if(c_n == n) {
    return(c_n * log(z + 1))
  }
  if((z * c_n / (n - c_n)) <= 1) {
    return(0)
  }
  TS <-  c_n * log(c_n / n) + (n - c_n) * log((n - c_n) / n) - c_n * log(1 / (z + 1)) - (n - c_n) * log(z / (z + 1))
  TS
}


# compute the maxSPRT test statistic for both Binomial and Poisson models.
compute_LLR <- function(data, cur_date, rate_years, risk_interval, model_type, nonrisk_interval = NULL, match_ratio = NULL) {
  
  data$age_days <- as.numeric(cur_date - data$dob)
  data <- data[data$age_days > data$vac_age + max(risk_interval, nonrisk_interval),]
  
  data$event_date <- data$admidt
  
  n_cases <- sum(data[!is.na(data$event_date) & data$event_date <= cur_date,]$case)
  
  if(model_type=="Poisson") {
    
    days_at_risk <- get_days_in_interval(data, risk_interval)
    
    mu_h0 <- days_at_risk * rate_years / 365.25
    TS <- LLR.Poisson(mu_t = mu_h0, c_t = n_cases)
    
    out <- c(TS, days_at_risk / 365.25 , mu_h0, n_cases)
    
  } else { # model_type = Binomial
    
    
    if(is.null(nonrisk_interval)) stop("nonrisk_interval must be defined")
    
    n_controls <- sum(data[!is.na(data$event_date) & data$event_date <= cur_date,]$control)
    n <- n_cases + n_controls
    TS <- LLR.Binomial(c_n = n_cases, n = n, z = match_ratio)
    rate_ratio <- ifelse(n==n_cases, Inf, (match_ratio * n_cases) / (n - n_cases))
    
    out <- c("test_stat" = TS , "z" = match_ratio, "controls" = n_controls, "cases" = n_cases, "n" = n, "RR" = round(rate_ratio, 2))
  }
  
  return(out)
}


# calculate the total nuber of days individuals spent within an age interval relative to their age at vaccination
# this is used to calculate the times spent in risk and nonrisk periods.
get_days_in_interval<- function(data, interval) {
  interval_length <- range_length(interval)
  interval_start_age <- data$vac_age + interval[1]
  days_spent_in_interval <- data$age_days - interval_start_age + 1
  days_spent_in_interval[days_spent_in_interval < 0] <- 0
  days_spent_in_interval[days_spent_in_interval > interval_length ] <- interval_length
  sum(days_spent_in_interval)
}



# utility functions
# -----------------


events2persondays <- function(events, rate_year) {
  365.25* events / rate_year
}

persondays2calendardays <- function(persondays, interval_size, cohort_size) {
  round(365.25 / interval_size * persondays / cohort_size )
}

is.in.interval <- function(intgr, interval) (intgr >= min(interval)) & (intgr <= max(interval))

seqR <- function(range) {
  seq(from = range[1], to = range[2])
}

range_length <- function(range) {
  length(seqR(range))
}


# plotting
# --------
  
signalPlot <- function(x, Y, ...) {
  x <- x[1:nrow(Y)]
  y <- Y[, "ts"]
  cv <- Y[,"c"][1]
  par(mar = c(4, 4, 4, 14), las = 2)
  ylim <- c(0, max(y, cv) + 1)
  detectionday <- which(y >= cv)[1]
  plot(x, y, type = "l" , xlab = "LLR", ylab = "day", ylim = ylim,...)
  abline(h = cv, col = "red")
  abline(v = x[detectionday], col = "green")
  par(xpd = T)
  legend(max(x),ylim[2],legend = c("maxSPRT test statistic", "critical value", 
                               paste0("signal: ", x[detectionday]," (day ",detectionday,")")),
         col=c("black","red","green"), lty =1)
  par(xpd = F)
}

signalPlot2 <- function(x, Y, stop = T , ...) {
  
  # Info for plot
  TS <- Y[, "ts"]
  cv <- Y[,"c"][1]
  detectionday <- which(TS >= cv)[1]
  end_day <- ifelse(is.na(detectionday) | !stop, length(TS) - 1, detectionday + 1)
  day <- x[1:end_day]
  RR <- Y[1:end_day, "RR"]
  TS <- TS[1:end_day]

  # graphica ladjustments
  par(mar = c(6, 4, 4, 4), las = 2)
  ylim <- c(0, max(TS, cv) + 3)
  
  # plot day vs test statistic
  plot(day, TS, type = "l" , xlab = "day", ylab = "LLR", ylim = ylim,...)
  
  # add rate ratio estimate
  lines(day, RR, lty = 2, col = "cornflowerblue")
  
  # horizontal line at critical value and 1 ( for rate ratio)
  abline(h = cv, col = "red")
  abline(h = 1, col = "grey95")
  # add critical value and detection day
  abline(v = day[detectionday], col = "green")
  
  # add legend
  par(xpd = T)
  legend("topright",legend = c("maxSPRT test statistic", "critical value", "rate ratio",
                                   paste0("signal: ", day[detectionday]," (day ",detectionday,")")),
         col=c("black","red","cornflowerblue", "green"), lty =c(1, 1, 2, 1))
  par(xpd = F)
}

ccPlot <- function(x, Y,  ...) {
  x <- x[1:nrow(Y)]
  y1 <- Y[, 4]
  y2 <- Y[, 5]
  par(mar = c(4, 4, 4, 10), las = 2)
  ylim <- c(0, max(y1, y2) + 10)
  plot(x, y1, type = "l", col = "red", ylim = ylim, ...)
  lines(x, y2, type = "l", col ="green")
  par(xpd=T)
  legend(max(x), ylim[2], legend = c("controls","cases"), col = c("green","red"), lty=1)
  par(xpd=F)
}

ppPlot <- function(x, Y, ...) {
  x <- x[1:nrow(Y)]
  p_0 <- 1/(Y[, "z"] +1)
  p <- Y[, "cases"] / Y[,"n"]
  plot(x, p_0, type = "l", col="blue", las = 2, ylim = c(0,1), ylab = "proportion of cases", xlab="", ...)
  lines(x, p, type="l", col="red")
  legend("topright", legend=c("expected","observed"), col=c("blue","red"), lty=1)
}
