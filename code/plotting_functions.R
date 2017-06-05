# plotting
# --------

# plots
# -----

interval_plot <- function(risk, case, ctrl, ...) {
  x <- 0:150
  h <- 10
  plot(x,x, type = "n", yaxt = "n", ylab = "", xlab = "time", ...)
  # risk window
  abline(v = h + risk[1], col = "red")
  abline(v = h + risk[2], col = "red")
  rect(h + case[1], 30, h + case[2], 40, col = "grey50", lwd = 0)
  rect(h + ctrl[1] - 1, 30, h + ctrl[2] - 1, 40, col ="grey80", lwd = 0)
  leg <- paste(c("actual risk period", "risk interval", "control interval"), lapply(list(risk, case, ctrl), paste, collapse = "-"))
  legend("topright", legend = leg, col = c("red", "grey50", "grey80"), lty = 1, lwd = c(1, 8,8))
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


