# calculating the critical values for the binomial maxSPRT

# auth: Tuomo Nieminen

#' Compute the maxSPRT test statistic
#'
#' Computes the maxSPRT test statistic for Binomial and Poisson likelihoods. Vectorized for efficiency.
#'
#' @param c_n : The observed events in the risk group. Can be a vector giving the cumulative observed events.
#' @param n : The ttotal number of observed events. can be a vector of same length as c_n giving the cumulative total number of events.
#' @param z : The expected proportion of controls to cases.
#' This is defined by the lengths of the risk and nonrisk periods : length(nonrisk) / length(risk)
#'
#' @return The value of the test statistic
#'
#' @examples
#' df <- data.frame(n_case = cumsum(rpois(n = 1000, lambda = 1)), n_ctrl = cumsum(rpois(n = 1000, lambda = 2)))
#' df$n <- rowSums(df)
#' TS <- LLR_Binomial(df$n_case, df$n, z = 3)
#' str(TS)
#' # num [1:1000] 1.3863 0.0521 0.0231 0 0.0362 ...
#'
#' @export
LLR_Binomial <- function(c_n, n, z) {
  ifelse(n == 0, 0,
         ifelse(c_n == n, c_n*log(z + 1),
                ifelse((z*c_n /(n - c_n)) <= 1, 0,
                       c_n * log(c_n/n) + (n-c_n) * log((n-c_n)/n)-c_n*log(1/(z+1))-(n-c_n)*log(z/(z + 1)))))
  
}

# Compute the state space S, transition matrix P and initial probabilities for the binomial markov chain. 
# P is a square matrix with all possible transitions (n, k) -> (n, k) ;  n = 0, 1, ..., N ; k = 0, 1, .., n.
#' @return List with 3 elements.
#' - v = initial probabilities for the states
#' - S = the state space as a 2 column matrix
#' - P = transition matrix
get_transition <- function(N, p) {
  
  # help function to build the sample space S
  sequence2 <- function(nvec) unlist(lapply(nvec, seq_len)) - 1
  
  # all possible states (n, k) in a 2-column matrix (S)
  states <- rep.int(0:N, times = 1:(N+1))
  S <- cbind(n = states, k = sequence2(1:(N+1)))
  n_states <- length(states)
  
  # base for the transition matrix P
  P <- matrix(data = 0, nrow = n_states, ncol = n_states, dimnames = list("n"= states, "k" = states))
  
  # the transition matrix is sparse, only having values c(1-p, p) in a certain place in each row. 
  # below we determine the places of those values.
  
  # states from which is it possible to leave
  m <- n_states - (N+1)

  # column index of the values c(1-p, p) for each row
  increase <- rep(1, m)
  increase[!duplicated(states[1:m])] <- 2
  increase[1] <- 1
  j <- cumsum(increase) + 1
  
  # indexes for 1-p and p as a 2-column matrix
  p1_indx <- matrix(c(1:m, j), ncol = 2)
  p2_indx <- p1_indx
  p2_indx[, 2] <- p2_indx[, 2] + 1
  
  # insert the values 1-p and p to P
  P[p1_indx] <- 1-p
  P[p2_indx] <- p
  
  # initial probabilitites of the states are c(1, 0, ... , 0)
  initial <- c(1, rep(0, n_states-1))
  
  return(list(v = initial, P = P, S = S))
}


# compute probabilities for each possible state (n, k) after N transitions, given by p = v %*% P^(N)
#' @return a vector of probabilities
get_state_probs <- function(P, N, v) {
  PS <- P
  for(i in 1:(N-1)) { 
      PS <- PS %*% P
  }
  probs <- v %*% PS
  probs
}


# Iteratively determine critical values for binomial maxSPRT.
# N = maximum observations, z = ratio of controls to cases, alpha =  type I error

# example usage:

# cv_binom(N = 10, z = 2)
# $cv
# [1] 2.903923
# 
# $alpha
# [1] 0.04770614
cv_binom <- function(N , z , alpha = 0.05) {
  
  if(N > 25) warning("For N this large, computation will be slow", immediate. = T)
  
  # probability of a case under the null hypothesis
  p <- 1/(1+z)
  
  # transition matrix for all possible states (n, k)
  TR <- get_transition(N, p = p)
  
  # a two column matrix of all possible (n, k) states
  S <- TR$S
  
  # test statistic value for all possible states
  L <- mapply(LLR_Binomial, S[, "k"], S[,"n"], z = z)
  
  # unique values of the test statistic
  test_stats <- sort(unique(L), decreasing = T)

  # for each possible test_stat value (starting from the highest), check which test statistic are >= than that test_statistic. 
  # then compute the typeI error by summing the probabilities of those critical states, using the transition matrix P
  new_a <- 0; i <- 1
  
  message("Iterating...")
  options(digits = 5)
  while(new_a <= alpha & i <= length(test_stats)) {
    ts <- test_stats[i]
    
    cat("\rTrying critical value ", ts, "... ")
    
    # logical vector indicating the states which reject the null
    C <- L >= ts
    
    # transition matrix containing the probabilities of transfering from state (n, k) to another state.
    P <- TR$P
    
    # absorbing states
    P[C, ] <- 0
    diag(P)[C] <- 1
    
    # get probabilities of being in each state after N steps
    p <- get_state_probs(P, N, v = TR$v)
    
    # sum the probabilities corresponding to the rejection states
    new_a <- sum(p[C])
    
    a <- ifelse(new_a > alpha, a, new_a)
    
    cat("alpha is ", new_a)
    i <- i  + 1
  }
  message("Returning the previous critical value and alpha")
  
  return(list(cv=test_stats[i-2], alpha = a))
}


