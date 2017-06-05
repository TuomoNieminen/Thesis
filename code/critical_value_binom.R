# calculating the critical values for the binomial maxSPRT

# auth: Tuomo Nieminen

# maxSPRT test statistic, Binomial likelihood (maximized likelihood ratio)
# k : observed events in risk group (cases). can be a vector
# n : total observed events
# z : expected proportion of controls to cases
LLR.Binomial.vec <- function(k, n, z) {
  
  if(n==0) return(rep(0, length(k)))
  
  TS <- numeric(length(k))
  
  equal <- k == n
  TS[equal] <- k[equal] * log(z + 1)
  k <- k[!equal]
  
  expected <- ( z * k / (n - k) ) <= 1
  TS[!equal][expected] <- 0
  k <- k[!expected]
  
  TS[!equal][!expected] <-  k * log(k / n) + (n - k) * log((n - k) / n) - k * log(1 / (z + 1)) - (n - k) * log(z / (z + 1))
  TS
}

# help function to build the sample space S
sequence2 <- function(nvec) unlist(lapply(nvec, seq_len)) - 1

# compute the transition matrix P for the binomial markov chain. 
# P is a square matrix with all possible transitions (n, k) -> (n, k) ;  n = 0, 1, ..., N ; k = 0, 1, .., n.
get_transition <- function(N, p) {
  
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


# compute probabilities for each possible state (n, k) after N steps, given by p = v %*% P^(N)
get_state_probs <- function(P, N, v) {
  PS <- P
  for(i in 1:(N-1)) { 
      PS <- PS %*% P
  }
  probs <- v %*% PS
  probs
}


# iteratively determine critical values for binomial maxSPRT.
# N = maximum observations, z = ratio of controls to cases, a =  alpha.

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
  L <- mapply(LLR.Binomial.vec, S[, "k"], S[,"n"], z = z)
  
  # unique values of the test statistic
  test_stats <- sort(unique(L), decreasing = T)

  # for each possible test_stat value (starting from the highest), check which test statistic are >= than that test_statistic. 
  # then compute the typeI error by summing the probabilities of those critical states, using the transition matrix P
  new_a <- 0; i <- 1
  while(new_a <= alpha & i <= length(test_stats)) {
    ts <- test_stats[i]
    
    # logical vector indicating the states which reject the null
    C <- L >= ts
    
    # transition matrix containing the probabilities of transfering from state (n, k) to another state.
    P <- TR$P
    
    # absorbing states
    P[C, ] <- 0
    diag(P)[C] <- 1
    
    #critical values of the absorbing states
    C_S <- S[C, ]
    
    # get probabilities of being in each state after N steps and sum the probabilities corresponding to the rejection states
    p <- get_state_probs(P, N, v = TR$v)
    new_a <- sum(p[C])
    
    a <- ifelse(new_a > alpha, a, new_a)
    i <- i  + 1
  }
  
  return(list(cv=ts, alpha = a))
}


