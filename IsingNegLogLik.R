

# Negative log-likelihood function for the Ising model with fixed parameters
neg_log_likelihood_ising <- function(X, 
                                     paras_opt, 
                                     J_fixed, 
                                     h_fixed) {
  # Basic Info
  p <- ncol(X)
  m <- nrow(X)
  
  # Get which parameters are fixed
  v_paras_fixed <- c(J_fixed[upper.tri(J_fixed)], 
                     h_fixed)
  paras_fixed <- !is.na(v_paras_fixed)
  
  # Get full parameter vector
  n_int <- p*(p-1)/2
  paras_full <- rep(NA, p+n_int)
  paras_full[paras_fixed] <- v_paras_fixed[paras_fixed]
  paras_full[!paras_fixed] <- paras_opt
  
  # Get parameters into matrix format
  J <- matrix(NA, p, p)
  J[upper.tri(J)] <- paras_full[1:n_int]
  J[lower.tri(J)] <- J[upper.tri(J)]
  diag(J) <- 0
  h <- paras_full[(n_int+1):length(paras_full)]
  
  # Compute energy term for observed data
  interaction_term <- rowSums((X %*% J) * X) / 2  # Avoid double counting
  threshold_term <- X %*% h
  log_likelihood <- sum(interaction_term + threshold_term)
  
  # Compute partition function Z (exact for small p)
  all_configs <- expand.grid(rep(list(c(-1, 1)), p))
  all_configs <- as.matrix(all_configs)
  energy_all_configs <- rowSums((all_configs %*% J) * all_configs) / 2 +
    all_configs %*% h
  Z <- sum(exp(energy_all_configs))  # Partition function
  
  # Negative log-likelihood
  - (log_likelihood - m * log(Z))
}


