# plotting
# --------

# plots
# -----

interval_plot <- function(risk, case, ctrl, ...) {
  x <- 0:150
  h <- 10
  ylim <- c(0, 150)
  plot(x, x, type = "n", yaxt = "n", ylab = "", ylim = ylim,  xlab = "time", ...)
  # risk window
  abline(v = h + risk[1], col = "red")
  abline(v = h + risk[2], col = "red")
  rect(h + risk[1], ylim[1] - 5,  h + risk[2], ylim[2] + 5, col = "salmon", border = NA, density = 15)
  rect(h + case[1], 30, h + case[2], 40, col = "grey50", lwd = 0)
  rect(h + ctrl[1] - 1, 30, h + ctrl[2] - 1, 40, col ="grey80", lwd = 0)
  leg <- paste(c("actual risk period", "risk interval", "control interval"), lapply(list(risk, case, ctrl), paste, collapse = "-"))
  legend("topright", legend = leg, col = c("red", "grey50", "grey80"), lty = 1, lwd = c(1, 8,8), bty = "n")
}

signalPlot2 <- function(x, Y, stop = T , right = 16, cex = 0.7,  ...) {
  
  # Info for plot
  TS <- Y[, "ts"]
  cv <- Y[,"c"][1]
  detectionday <- which(TS >= cv)[1]
  end_day <- ifelse(is.na(detectionday) | !stop, length(TS) - 1, detectionday + 1)
  day <- x[1:end_day]
  RR <- Y[1:end_day, "RR"]
  TS <- TS[1:end_day]
  
  # graphical adjustments
  par(mar = c(6,4,4, right), las = 2, cex = cex)
  ylim <- c(0, max(TS, cv) + 5)
  xlim <- c(min(day), max(day) + 30)
  
  # plot day vs test statistic
  plot(day, TS, type = "n" , xlab = "day", ylab = "LLR", xlim = xlim, ylim = ylim, bty = "l", ...)
  
  # horizontal line at 1 (reference for RR)
  abline(h = 1, col = "grey90")
  
  # horizontal line at critical value 
  abline(h = cv, col = "red")
  
  # rate ratio estimate
  lines(day, RR, lty = 2, col = "cornflowerblue")
  
  # test statistic
  lines(day, TS)

  # add critical value and detection day
  abline(v = day[detectionday], col = "green")
  margin <- as.numeric(xlim[2] - xlim[1])*0.05
  
  # add legend
  legend(xlim[2] + margin , ylim[2] + 2, xpd = T, legend = c("maxSPRT test statistic", "critical value", "rate ratio", "RR = 1",
                               paste0("signal: ", day[detectionday]," (day ",detectionday,")")), bty = "n",
         col=c("black","red","cornflowerblue", "grey95", "green", "grey90"), lty =c(1, 1, 2, 1, 1))
}


