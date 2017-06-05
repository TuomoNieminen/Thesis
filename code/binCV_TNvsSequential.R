# compute critical values using sequential package method and method implemented by Tuomo Nieminen

setwd("~/TNversio/vaccsafety/Signal/Simulations")
source("critical_value_binom.R")

library(Sequential)
N <- c(10, 15, 20)
A <- c(0.05, 0.01, 0.001)
Z <- c(2, 3)

methods <- c("Sequential","TuomoN")

critical <- array(dim = c(length(N), length(A), length(Z), 2), dimnames = list("N"=N, "alpha"=A, "z"=Z, "method"=methods))

# Sequential package method

seq_time <- system.time({
  for(i in seq_along(N)) {
    for(j in seq_along(A)) {
      for(k in seq_along(Z)) {
        critical[i, j, k, 1] <- CV.Binomial(N = N[i], alpha = A[j], z = Z[k])$cv
      }
    }
  }
})

# TuomoN method

tn_time <- system.time({
  for(i in seq_along(N)) {
    for(j in seq_along(A)) {
      for(k in seq_along(Z)) {
        critical[i, j, k, 2] <- cv_binom(N = N[i], alpha = A[j], z = Z[k])$cv
      }
    }
  }
})
  
critical

cat("Sequential time: \n")
seq_time
cat("TuomoN time: \n")
tn_time

