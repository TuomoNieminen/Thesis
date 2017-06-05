# apply the maxSPRT method retrospectively to mpr, rota and pcv vaccination data and seizure hospitalization data

# auth. Tuomo Nieminen

# published at http://opus.thl.fi/roko/stat/vaccsafety/Seizure.signal.detection.html

source("~/TNversio/vaccsafety/EUSAFE/EUSAFE_functions.R")
source("~/TNversio/vaccsafety/Signal/maxSPRT.R")

library(eRos)
library(eRosPublish)

op <- par(no.readonly = T)


picpath <- "/group/rokostat/www/vaccsafety/"
out <- list()


# INFO

out$About <- paste0(
  "This is an experiment with a signal detection method related to vaccine safety surveillance. 
  The method used is the binomial maxSPRT method introduced by Kulldorff (Kulldorff 2011).")

out$Data <- paste0(
  "The vaccinations of interest are the first doses of MPR, PCV and Rota vaccinations for children under 2 years of age born during 2010 - 2014. 
  For the MPR group, the 2014 cohort was exluded because the expected age at vaccination is higher and no sufficient data is yet available.",
  "The possible adverse events of interest are hospitalizations with icd10 diagnosis codes related to seizures. ",
  "For both vaccinations and hospitalizations, age intervals that cover first dose vaccinations of the chosen vaccines were concidered: 
mpr (250, 650), rota (30, 140), pcv (60, 200). ",
  "For more details on the icd10 codes see a previous <a href='https://terho.thl.fi/wiki01/display/rokobiom/Case-series+analyysi'>case series analysis.</a>
  </p>")

out$Method <- paste0(
  
  "<p>",
  "The method used is the binomial variant of maxSPRT (Kulldorff 2011). 
From an epidemiological study standpoint, it is a self-controlled case-control method, 
where the surveillance period for each child is divided into case and control periods. 
An event occuring during the case period is classified as a case and an event occuring during the control period as a control.
The method uses only data on cases. </p>",
  
  "<p>",
  "In the binomial maxSPRT the assumption under the null hypothesis (H0) is that 
the probability of an event classified as a case is proportional to the time spent in the case interval. 
For example, if the case period is 7 days and the control period 14 days and the person with an event has fully passed through both intervals, 
then the (H0) probability of the event classified as a case is 7 / (7 + 14) = 1/3. 
The alternative hypothesis (H1) is that the probability of a case is higher, which corresponds to an added risk of an event during the risk period. </p>",
  
  "<p>",
  "The test statistic (ts) is the maximized likelihood ratio L(p|H1) / L(p|H0). 
Critical values were computed numerically using the Sequential R-package. 
For each day of surveillance the test statistic is calculated and compared to the critical value. 
If the test statistic equals or exceeds the critical value at any point, H1 is accepted and a signal is detected. </p>",
  
  "<p>",
  "You can further explore the method with a <a href= 'http://shiny.app.thl.fi/maxSPRT/' >maxSPRT web-app</a>.</p>")

out$`Further assumptions` <- paste0(
  
  "<p>",
  "An estimate for the baseline incidence for each age group described above was used to compute the 
expected number of events during surveillance under H0, and this was then used as the upper limit of surveillance (N).
In this experiment the incidences were calculated by summing the events during the ages of interest (excluding the risk period) 
and estimating the population at risk (vaccinated) using exact population data and the estimated vaccination coverages of the chosen vaccines. </p>",
  
  "<p>",
  "The case interval was chosen to be 1 - 14 days after vaccination and the control interval 15 - 42 days after vaccination. 
  Choosing the control interval to be after vaccination gets around the healthy vaccinee problem, 
  but also assumes that the possible effect of the vaccination is transient.</p>",
  "<p>",
  "Type I error (false signal) probability is 0.05. </p>")


#### ANALYSIS

# HARDCODED VARIABLES

# target cohorts
cohorts <- 2010:2014
cohorts_older <- 2010:2013

# birth data from http://www.stat.fi/til/synt/2014/synt_2014_2015-04-14_tie_001_fi.html
pops<- c("2010" = 60980,	"2011" = 59961, "2012" = 59493, "2013" = 58134, "2014" = 57232)
popsize <- sum(pops)
popsize_older <- sum(pops[paste(cohorts_older)]) 

start_day = "2010-01-01"

# age groups loosely based on age distributions here: https://terho.thl.fi/wiki01/display/rokobiom/Case-series+analyysi
mpr_age  <- c(250, 650)
rota_age <- c(30, 140)
pcv_age  <- c(60, 200)

# vaccine coverages. source: http://shiny.app.thl.fi/childhoodvaccination/
mpr_coverage  <- 0.95
rota_coverage <- 0.93
pcv_coverage  <- 0.94

# case and control intervals
case_interval <- c(1, 14)
ctrl_interval <- c(15, 42)

match_ratio <- range_length(ctrl_interval) / range_length(case_interval)
obs_days <- range_length(case_interval) + range_length(ctrl_interval)


#### HOSPITALIZATION DATA

# hospitalization data (seizure diagnoses)
diagnoses <- c("A" = "Viral infections \n of the central nervous system \n icd10 A858,A86,A87,A88,A89",
               "G" = "Inlfammatory diseases \n of the central nervous system \n icd10 G038,G039,G04,G05",
               "R" = "General symptoms and signs \n icd10 R291,R55,R560,R568")
hosp <- get_hosp(diagnoses, cohorts = cohorts, start_day = start_day)



#### VACCINATION DATA

# vaccination data (mpr, rota, pcv)
mpr  <- get_vac_adv(52, "MPR",  cohorts = cohorts_older, start_day = start_day)
rota <- get_vac_adv(55, "Rota", cohorts = cohorts, start_day = start_day)
pcv  <- get_vac_adv(40, "PCV",  cohorts = cohorts, start_day = start_day)

mpr  <- mpr[is.in.interval(mpr$vac_age, mpr_age),]
rota <- rota[is.in.interval(rota$vac_age, rota_age),]
pcv  <- pcv[is.in.interval(pcv$vac_age, pcv_age),]

# Add closest vaccination date
mpr_hosp  <- add_vaccination(hosp, mpr, ageinterval = mpr_age)
rota_hosp <- add_vaccination(hosp, rota, ageinterval = rota_age)
pcv_hosp  <- add_vaccination(hosp, pcv, ageinterval = pcv_age)

# add case / control information
mpr_hosp  <- classify_data(mpr_hosp, case_interval = case_interval, control_interval = ctrl_interval)
rota_hosp <- classify_data(rota_hosp, case_interval = case_interval, control_interval = ctrl_interval)
pcv_hosp  <- classify_data(pcv_hosp, case_interval = case_interval, control_interval = ctrl_interval)

# non-case hospitalizations during agegroups to estimate baseline incidence
events_mprgroup  <- sum(mpr_hosp$case==0)
events_rotagroup <- sum(rota_hosp$case==0)
events_pcvgroup  <- sum(pcv_hosp$case==0)

# exlude the nonvaccinated for later analysis
mpr_hosp  <- mpr_hosp[mpr_hosp$vaccinated, ]
rota_hosp <- rota_hosp[rota_hosp$vaccinated,]
pcv_hosp  <- pcv_hosp[pcv_hosp$vaccinated, ]


#### INCIDENCE and expected events

# # Estimate the incidence of the event for power calculations

# person years spent at non-case age groups
mpr_PY  <- popsize_older * (range_length(mpr_age) - range_length(case_interval)) / 365.25
rota_PY <- popsize* (range_length(rota_age) - range_length(case_interval)) / 365.25
pcv_PY  <- popsize * (range_length(pcv_age) - range_length(case_interval)) / 365.25

# incidence per personyears by agegroup.
incidence_mprgroup  <- events_mprgroup / mpr_PY
incidence_rotagroup <- events_rotagroup / rota_PY
incidence_pcvgroup  <- events_pcvgroup / pcv_PY

# expected total events during observation periods
mpr_expected  <- round( popsize_older * mpr_coverage *  obs_days/365.25 * incidence_mprgroup )
rota_expected <- round( popsize * rota_coverage * obs_days/365.25 * incidence_rotagroup )
pcv_expected  <- round( popsize * pcv_coverage *  obs_days/365.25 * incidence_pcvgroup )


#### maxSPRT METHOD

library(Sequential)

# critical values
cv_mpr  <- CV.Binomial(mpr_expected,  alpha = 0.05, M = 1, z = match_ratio)$cv
cv_rota <- CV.Binomial(rota_expected, alpha = 0.05, M = 1, z = match_ratio)$cv
cv_pcv  <- CV.Binomial(pcv_expected,  alpha = 0.05, M = 1, z = match_ratio)$cv

# Performance (power) assuming a risk ratio of 1.5 
RR <- 1.5
power_mpr <- Performance.Binomial(N = mpr_expected, M = 1, cv = cv_mpr, z = match_ratio, RR = RR)$Power
power_rota <- Performance.Binomial(N = rota_expected, M = 1, cv = cv_rota, z = match_ratio, RR = RR)$Power
power_pcv <- Performance.Binomial(N = pcv_expected, M = 1, cv = cv_pcv, z = match_ratio, RR = RR)$Power


conditions <- round(rbind(MPR = c(N = mpr_expected, c = cv_mpr, power = power_mpr, alpha = 0.05),
                           ROTA= c(N = rota_expected, c = cv_rota, power = power_rota, alpha = 0.05),
                           PCV = c(N = pcv_expected, c = cv_pcv, power = power_pcv, alpha = 0.05)),2)

out$`Surveillance conditions` <- paste(eosMatrix(conditions),
                                       "<p> The power is calculated for assumed risk ratio of 1.5 </p>")

# MONITORING

dates <- seq.Date(from = as.Date(paste0(min(cohorts),"-01-01")), to = max(hosp$admidt) + obs_days, by = "days")

# mpr_end_day  <- round(length(cohorts)*365.25) + max(mpr_age) + max(ctrl_interval)
# rota_end_day <- round(length(cohorts)*365.25) + max(rota_age) + max(ctrl_interval)
# pcv_end_day  <- round(length(cohorts)*365.25) + max(pcv_age) + max(ctrl_interval)

mpr_test  <- matrix(data = NA, nrow = length(dates), ncol = 6, dimnames = list("day"=NULL, c("ts","z","controls","cases", "n", "RR")))
rota_test <- matrix(data = NA, nrow = length(dates), ncol = 6, dimnames = list("day"=NULL, c("ts","z","controls","cases", "n", "RR")))
pcv_test  <- matrix(data = NA, nrow = length(dates), ncol = 6, dimnames = list("day"=NULL, c("ts","z","controls","cases", "n", "RR")))

# compute the log-likelihood-ratio (LLR) test statistic for each day (then compare to critical value)
# use the expected number of events as the upper length of surveillance.

for(i in seq_along(dates)) {
  day <- dates[i]
  
  if(max(mpr_test[,"n"], na.rm = T) < mpr_expected ) {
    mpr_test[i, ]  <- compute_LLR(mpr_hosp,  cur_date = day, rate_years = incidence_mprgroup,  risk_interval = case_interval, nonrisk_interval = ctrl_interval, match_ratio = match_ratio, model_type = "Binomial")
  }
  if(max(rota_test[,"n"], na.rm = T) < rota_expected ) {
    rota_test[i, ] <- compute_LLR(rota_hosp, cur_date = day, rate_years = incidence_rotagroup, risk_interval = case_interval, nonrisk_interval = ctrl_interval, match_ratio = match_ratio, model_type = "Binomial")
  }
  if(max(pcv_test[,"n"], na.rm = T) < pcv_expected ) {
    pcv_test[i, ]  <- compute_LLR(pcv_hosp,  cur_date = day, rate_years = incidence_pcvgroup,  risk_interval = case_interval, nonrisk_interval = ctrl_interval, match_ratio = match_ratio, model_type = "Binomial")
  }
}

mpr_test <- data.frame(mpr_test[complete.cases(mpr_test),])
rota_test <- data.frame(rota_test[complete.cases(rota_test),])
pcv_test <- data.frame(pcv_test[complete.cases(pcv_test),])

mpr_test$day <- rownames(mpr_test)
rota_test$day <- rownames(rota_test)
pcv_test$day <- rownames(pcv_test)

mpr_test$signal <- mpr_test[,"ts"] >= cv_mpr
rota_test$signal <- rota_test[, "ts"] >= cv_rota
pcv_test$signal <- pcv_test[,"ts"] >= cv_pcv

mpr_test$c <- cv_mpr
rota_test$c <- cv_rota
pcv_test$c <- cv_pcv

#### RESULTS

# Summary statistics
test_results <- list(MPR = mpr_test, Rota = rota_test, PCV = pcv_test)
detection <- lapply(test_results, function(df){
  day <- which(df$signal)[1]
  if(is.na(day)) day <- nrow(df)
  detection_info <- df[day, c("signal","day", "ts","cases","controls", "n", "RR")]
  detection_info
})

detection <- do.call(rbind, detection)
detection$ts <- round(detection$ts, 2)
colnames(detection)[3] <- "LLR"

application_results <- list(dates = dates, results = test_results, conditions = conditions, summary = detection)
# save(file = "/home/tnin/TNversio/vaccsafety/Gradu_TN/Presentation/application_results.Rda", application_results)

out$`End of surveillance summaries` <- eosMatrix(detection)


# Plots

op <- par(no.readonly = T)

# critical values and signal detection days
png(file=paste0(picpath, "seizure_signal.png"), width = 900, height = 900)
par(mfrow = c(3, 1), cex = 1.3)
signalPlot(dates, mpr_test,  main = "MPR")
signalPlot(dates, rota_test, main = "Rota")
signalPlot(dates, pcv_test,  main ="PCV")
par(op)
dev.off()

out$`Rota, MPR, PCV seizure signal detection` <- eosThumb("seizure_signal.png")


# cumulative cases and controls
png(file=paste0(picpath, "seizure_case_proportions.png"), width = 900, height = 900)
par(mfrow=c(3,1), cex = 1.3)
ppPlot(dates, mpr_test,  main = "MPR")
ppPlot(dates, rota_test, main = "Rota")
ppPlot(dates, pcv_test,  main ="PCV")
par(op)
dev.off()

out$`Rota, MPR, PCV expected and observed proportions of cases` <- eosThumb("seizure_case_proportions.png")



##### PUBLISH #####

source("/group/rokostat/repository/common/headers_common.R")
script.name <- paste0(getwd(), "/seizure_signal_detection.R")
outpath <- "/group/rokostat/www/vaccsafety/"
Publish(out, path = outpath, title = "Seizure signal detection")
