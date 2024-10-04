library(rstan)
library(CompressiveNMF)
library(LaplacesDemon)


source("R/generate_synthetic_mutation_counts.R")
source("R/generate_synthetic_mutation_counts.R")
source("R/plot_signature.R")
source("R/initialize_prior.R")


DirichletReg_stan <- function(y, X, a){
  stan_data <- list(S = nrow(y), I = ncol(y), J = ncol(X), X = X, y = y, a = a)
  model_file <- "stan/simplex_regression.stan"  
  fit_stan_sig1 <- stan(file = model_file, 
                        data = stan_data, 
                        iter = 2000, chains = 1)
  res <- extract(fit)
  return(res)
}



log_posterior <- function(w, phi, logSumY, X, a, a0, b0) {
  S <- length(logSumY); I <- nrow(X)
  # Find dot product
  dot_w_X <- X %*% w
  # log density
  log_target <- S * lgamma(phi) - S * sum(lgamma(phi * dot_w_X)) + sum((phi * dot_w_X - 1) * logSumY) +
    (a - 1) * sum(log(w)) + dgamma(phi, a0, b0, log = TRUE)
  return(log_target)
} 

DirichletReg <- function(y, X, a = 1, a0 = 1, b0 = 1,
                         nsamples = 2000, burnin = 2000, 
                         method = "max_density", fixed_phi = FALSE, phi = 1000,
                         alpha_proposal = 100, V = 6e-7){
  # Data dimension
  S <- ncol(y); I <- nrow(y); J <- ncol(X)
  # Useful quantities
  logSumY <- colSums(log(y))
  
  # Storing values
  W <- matrix(NA, nrow = nsamples, ncol = J)
  PHI <- rep(NA, nsamples)
  ACCEPT <- matrix(NA, nrow = nsamples, ncol = 2)
  colnames(ACCEPT) <- c("w", "phi")
  
  # Initial value for the sampler
  w <- rep(1/J, J)
  #log_phi <- mean_log_phi <- log(phi)
  mean_phi <- phi
  var_phi <- 1
  
  # Run the sampler
  pb <- txtProgressBar(style=3)
  for(r in 1:(nsamples + burnin)) {
    setTxtProgressBar(pb, r/(nsamples + burnin))
    acc_w <- acc_phi <- 0
    
    #------------------------------------------------ Step 1. Update w
    # Propose the new value
    if(method == "max_density"){
      # Propose the maximum density approach
      alpha_new <- dens.max.dir(center = w, V = V)  
    } else {
      alpha_new <- alpha_proposal * w
    }
    # Propose a new value
    w_new <- c(rdirichlet(1, alpha = alpha_new)) + 1e-12
    w_new <- w_new/sum(w_new)
    # Calculate log target under the two values
    lp_old <- log_posterior(w, phi, logSumY, X, a, a0, b0)
    lp_new <- log_posterior(w_new, phi, logSumY, X, a, a0, b0)
    
    # Calculate MH ratio correction for symmetry
    if(method == "max_density"){
      alpha_old <- dens.max.dir(center = w_new, V = V)  
      sym_correction <- LaplacesDemon::ddirichlet(w, alpha_new, log = TRUE) - 
        LaplacesDemon::ddirichlet(w_new, alpha_old, log = TRUE)
    } else {
      sym_correction <- LaplacesDemon::ddirichlet(w, alpha_proposal * w_new, log = TRUE) - 
        LaplacesDemon::ddirichlet(w_new, alpha_proposal * w, log = TRUE)
    }
    
    # Calculate metropolis ratio
    log_MH_ratio_w <- lp_new - lp_old + sym_correction
    
    if(log(runif(1)) <= log_MH_ratio_w){
      w <- w_new
      acc_w <- 1
    }
    
    #------------------------------------------------ Step 2. Update Phi
    if(!fixed_phi){
      if(r > 1){
        # mean_log_phi_new <- (r - 1) / r * mean_log_phi + log_phi / r
        # var_phi <- (r - 2) / (r - 1) * var_phi + 
        #   1 / (r - 1) * log_phi^2 - 
        #   r / (r-1) * mean_log_phi_new^2 + 
        #   mean_log_phi^2
        # mean_log_phi <-  mean_log_phi_new
        mean_phi_new <- (r - 1) / r * mean_phi + phi / r
        var_phi <- (r - 2) / (r - 1) * var_phi + 
          mean_phi^2 + 
          1 / (r - 1) * phi^2 - 
          r / (r-1) * mean_phi_new^2
        mean_phi <-  mean_phi_new
      }
      
      # Propose the new value for phi
      var_proposal <- ifelse(r >= 11, (2.38)^2 * (var_phi + 0.05), 25)
      #log_phi_new <- rnorm(1, mean = phi, sd = sqrt(var_proposal))
      phi_new <- rnorm(1, mean = phi, sd = sqrt(var_proposal))
      
      # Calculate log_posterior in both cases
      # lp_old_phi <- log_posterior(w, exp(log_phi), logSumY, X, a, a0, b0)
      # lp_new_phi <- log_posterior(w, exp(log_phi_new), logSumY, X, a, a0, b0)
      lp_old_phi <- log_posterior(w, phi, logSumY, X, a, a0, b0)
      lp_new_phi <- log_posterior(w, phi_new, logSumY, X, a, a0, b0)
      
      # Calculate MH ratio
      log_MH_ratio_phi <- lp_new_phi - lp_old_phi
      
      if(log(runif(1)) <= log_MH_ratio_phi){
        #log_phi <- log_phi_new
        phi <- phi_new
        acc_phi <- 1
      }
      #phi <- exp(log_phi)
    }
    
    #------------------------------------------------ Step 3. Save the output
    if(r > burnin){
      W[r - burnin, ] <- w
      PHI[r - burnin] <- phi
      ACCEPT[r - burnin, 1] <- acc_w
      ACCEPT[r - burnin, 2] <- acc_phi
    }
  }
  close(pb)
  colnames(W) <- colnames(X)
  return(list(w = W, phi = PHI, accept = ACCEPT))
}