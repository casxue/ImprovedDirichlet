Compute_sensitivity_precision <- function(R_hat, R_true, cos_cutoff = 0.9){
  # Sensitivity: proportion of ground truth signatures that were estimated with 
  #              sufficiently high similarity
  sig_sens <- sapply(1:ncol(R_true), function(i) 
    max(sapply(1:ncol(R_hat), function(j) cosine(R_true[, i], R_hat[, j]))))
  # Precision: proportion of estimated signatures that were sufficiently similar to 
  #            ground truth signatures
  sig_prec <- sapply(1:ncol(R_hat), function(i) 
    max(sapply(1:ncol(R_true), function(j) cosine(R_true[, j], R_hat[, i]))))
  return(c("Sensitivity" = mean(sig_sens > cos_cutoff), "Precision" = mean(sig_prec > cos_cutoff)))
}

# Function to find the true value of Lambda
get_Lambda_Comp <- function(resComp,  samples = NULL) {
  Lambda <- 0
  id <- resComp$selected_chain
  if(is.null(samples)){
    nsamples <- dim(resComp$mcmc_out[[id]]$Signatures)[1]
    samples <- 1:nsamples
  }
  for (i in samples) {
    Lambda <- Lambda + resComp$mcmc_out[[id]]$Signatures[i, , ] %*% resComp$mcmc_out[[id]]$Weights[i, , ]
  }
  return(Lambda / length(samples))
}

get_cosine_similarity <- function(matchedSign){
  sims <- sapply(1:ncol(matchedSign$R_hat), function(i) cosine(matchedSign$R_hat[, i], matchedSign$R_true[, i]))
  return(sims) 
}

# Match the mutational signatures with the true ones using the Hungarian algorithm 
# to calculate the average cosine similarity 
match_MutSign <- function(R_true, R_hat) {
  # Need to make sure that R_hat and R_true are the same
  k_true <- ncol(R_true)
  k_hat <- ncol(R_hat)
  k_tot <- max(c(k_hat, k_true))
  I <- nrow(R_true)
  mat0 <- matrix(100, nrow = I, ncol = abs(k_hat - k_true)) # 
  if (k_hat > k_true) {
    colnames(mat0) <- paste0("new_extra", 1:ncol(mat0))
    R_true <- cbind(R_true, mat0)
  } else if (k_hat < k_true) {
    R_hat <- cbind(R_hat, mat0)
  }
  
  # Match mutational signatures using the Hungarian algorithm
  CosMat <- matrix(1, k_tot, k_tot)
  for (i in 1:k_tot) {
    for (j in 1:k_tot) {
      CosMat[i, j] <- 1 - cosine(R_true[, i], R_hat[, j])
    }
  }
  match <- RcppHungarian::HungarianSolver(CosMat)$pairs[, 2]
  R_hat_matched <- R_hat[, match]
  
  # Change 100 with 0s
  R_hat_matched[R_hat_matched == 100] <- 0
  R_true[R_true == 100] <- 0
  
  return(list("R_hat" = R_hat_matched, "R_true" = R_true, "match" = match))
}


match_Theta <- function(Theta_true, Theta_hat, match){
  # Need to make sure that Theta_hat and Theta_true are the same
  k_true <- nrow(Theta_true)
  k_hat <- nrow(Theta_hat)
  k_tot <- max(c(k_hat, k_true))
  J <- ncol(Theta_true)
  mat0 <- matrix(0, nrow = abs(k_hat - k_true), ncol = J)
  if (k_hat > k_true) {
    rownames(mat0) <- paste0("new_extra", 1:nrow(mat0))
    Theta_true <- rbind(Theta_true, mat0)
  } else if (k_hat < k_true) {
    Theta_hat <- rbind(Theta_hat, mat0)
  }
  Theta_hat <- Theta_hat[match, ]
  return(Theta_hat)
}
# Calculate the MSE between the true Theta and the estimated one
compute_RMSE_Theta <- function(Theta_true, Theta_hat, match){
  # Need to make sure that Theta_hat and Theta_true are the same
  k_true <- nrow(Theta_true)
  k_hat <- nrow(Theta_hat)
  k_tot <- max(c(k_hat, k_true))
  J <- ncol(Theta_true)
  mat0 <- matrix(0, nrow = abs(k_hat - k_true), ncol = J)
  if (k_hat > k_true) {
    rownames(mat0) <- paste0("new_extra", 1:nrow(mat0))
    Theta_true <- rbind(Theta_true, mat0)
  } else if (k_hat < k_true) {
    Theta_hat <- rbind(Theta_hat, mat0)
  }
  Theta_hat <- Theta_hat[match, ]
  # Now, match the matrices
  is_novel <- grepl("new", rownames(Theta_true))
  # RMSE between true cosmic 
  if(any(!is_novel)){
    if(sum(!is_novel) == 1){
      Theta_temp <- t(Theta_hat[!is_novel, ])
      Theta_temp_true <- t(Theta_true[!is_novel, ])
    } else {
      Theta_temp <- as.matrix(Theta_hat[!is_novel, ])
      Theta_temp_true <- as.matrix(Theta_true[!is_novel, ])
    }
    rmse_cosmic <- mean(sapply(1:sum(!is_novel), function(i) 
      sqrt(mean((Theta_temp[i, ] - Theta_temp_true[i, ])^2))))
  } else {
    rmse_cosmic  <- -1 
  }
  # RMSE between novel
  if(any(is_novel)){
    if(sum(is_novel) == 1){
      Theta_temp <- t(Theta_hat[is_novel, ])
      Theta_temp_true <- t(Theta_true[is_novel, ])
    } else {
      Theta_temp <- as.matrix(Theta_hat[is_novel, ])
      Theta_temp_true <- as.matrix(Theta_true[is_novel, ])
    }
    rmse_novel <- mean(sapply(1:sum(is_novel), function(i) 
      sqrt(mean(Theta_temp[i, ] - Theta_temp_true[i, ])^2)))
  } else {
    rmse_novel  <- -1 
  }
  rmse_all <- c(rmse_novel, rmse_cosmic)
  
  return(c("rmse_novel_Theta" = rmse_novel, "rmse_cosmic_Theta" = rmse_cosmic, 
           "rmse_Weights" = mean(rmse_all[rmse_all!=-1])))
  return(rmse)
}


compute_RMSE_Signature <- function(R_hat, R_true){
  is_novel <- grepl("new", colnames(R_true))
  # RMSE between true cosmic 
  if(any(!is_novel)){
    rmse_cosmic <- mean(sapply(1:sum(!is_novel), function(i) 
      sqrt(mean((as.matrix(R_hat[, !is_novel])[, i] - as.matrix(R_true[, !is_novel])[, i])^2))))
  } else {
    rmse_cosmic  <- -1 
  }
  # RMSE between novel
  if(any(is_novel)){
    rmse_novel <- mean(sapply(1:sum(is_novel), function(i) 
      sqrt(mean((as.matrix(R_hat[, is_novel])[, i] - as.matrix(R_true[, is_novel])[, i])^2))))
  } else {
    rmse_novel  <- -1 
  }
  rmse_all <- c(rmse_novel, rmse_cosmic)
  
  return(c("rmse_novel_Sig" = rmse_novel, "rmse_cosmic_Sig" = rmse_cosmic, 
           "rmse_Signatures" = mean(rmse_all[rmse_all!=-1])))
}

Postprocess_Compressive <- function(resComp, data) {
  # Step 1 - calculate the number of inferred signatures
  K <- length(resComp$RelWeights)
  # Step 2 - Calculate RMSE wrt to the count matrix and the rate matrix
  Lambda <- get_Lambda_Comp(resComp)
  Lambda_true <- data$Rmat %*% data$Theta
  rmse_Lambda <- sqrt(mean((Lambda_true - Lambda)^2))
  rmse_Counts <- sqrt(mean((data$X - Lambda)^2))
  # Step 3 - calculate the cosine similarity between the true and the inferred signatures
  R_true <- apply(data$Rmat, 2, function(x) x / sum(x))
  R_hat <- resComp$Signatures
  matchedSign <- match_MutSign(R_true = R_true, R_hat = R_hat)
  cos_sim <- mean(get_cosine_similarity(matchedSign))
  # Step 4 - calculate the RMSE between Theta and the rest
  Theta_hat <- resComp$Weights
  rmse_R <- compute_RMSE_Signature(R_hat = matchedSign$R_hat, R_true = matchedSign$R_true)  
  rmse_Theta <- compute_RMSE_Theta(Theta_true = data$Theta, Theta_hat = Theta_hat, matchedSign$match)  
  # Step 5 - calculate the sensitivity and precision
  sens_prec  <- Compute_sensitivity_precision(R_hat = R_hat, R_true)
  return(list(
    Lambda = Lambda,
    R_hat = R_hat,
    Theta_hat = Theta_hat,
    Mu_hat = resComp$RelWeights,
    signatures = matchedSign,
    results = c("K" = K, 
                "rmse_Lambda" = rmse_Lambda, 
                "rmse_Counts" = rmse_Counts, 
                rmse_R, 
                rmse_Theta, 
                sens_prec,
                "cos_sim" = cos_sim, 
                "time" = resComp$time)
  ))
}

