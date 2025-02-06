

# Function to estimate J and h while allowing fixed values
estimate_ising_mle <- function(X, 
                               J_fixed = NULL, 
                               h_fixed = NULL, 
                               method = "BFGS") {
  
  p <- ncol(X)
  
  # In case of no constraints
  if(is.null(J_fixed)) J_fixed <- matrix(NA, p, p)
  if(is.null(h_fixed)) h_fixed <- rep(NA, p)
  
  # Initial Values
  n_int <- p*(p-1)/2
  J <- matrix(NA, p, p)
  J[upper.tri(J)] <- runif(n_int)
  h <- runif(p)
  
  # Collapse into vector
  paras <- c(J[upper.tri(J)], 
             h)
  
  # Collapse fixed parameters in to vector
  v_paras_fixed <- c(J_fixed[upper.tri(J)], 
                     h_fixed)
  paras_fixed <- !is.na(v_paras_fixed)
  
  paras_opt <- paras[!paras_fixed]

  # Optimize using `optim()`
  result <- optim(paras_opt, 
                  neg_log_likelihood_ising, 
                  X = X, 
                  J_fixed = J_fixed, 
                  h_fixed = h_fixed, 
                  method = method, 
                  control = list(maxit = 1000))
  
  # Get full parameter vector
  n_int <- p*(p-1)/2
  paras_full <- rep(NA, p+n_int)
  paras_full[paras_fixed] <- v_paras_fixed[paras_fixed]
  paras_full[!paras_fixed] <- result$par
  
  # Get parameters into matrix format
  J <- matrix(NA, p, p)
  J[upper.tri(J)] <- paras_full[1:n_int]
  J[lower.tri(J)] <- J[upper.tri(J)]
  diag(J) <- 0
  h <- paras_full[(n_int+1):length(paras_full)]
  
  outlist <- list(J = J, 
                  h = h, 
                  logLik = -result$value, 
                  convergence = result$convergence)
  
  return(outlist)
} # eoF












