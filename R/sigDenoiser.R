library(tidyverse)
library(rstan)
library(CompressiveNMF)
library(LaplacesDemon)
library(coda)
library(nnls)
library(JuliaCall)

# Source Julia here
julia <- julia_setup(installJulia = TRUE)
#julia_command('import Pkg; Pkg.add("Distributions")')
#julia_command('import Pkg; Pkg.add("SpecialFunctions")')
#julia_command('import Pkg; Pkg.add("LinearAlgebra")')
julia_source("Julia/max_density.jl")

# Source useful functions here
source("R/generate_synthetic_mutation_counts.R")
source("R/generate_synthetic_mutation_counts.R")
source("R/plot_signature.R")
source("R/initialize_prior.R")

DirichletReg_stan <- function(y, X, a){
  stan_data <- list(S = nrow(y), I = ncol(y), J = ncol(X), X = X, y = y, a = a)
  model_file <- "stan/simplex_regression.stan"  
  fit <- stan(file = model_file, 
                        data = stan_data, 
                        iter = 2000, chains = 1)
  res <- rstan::extract(fit)
  return(res)
}

plot.trace <- function(W){
  plot(W[, 1], type = "l", ylim = c(0, max(W) + 1e-3))
  for(j in 2:ncol(W)){
    lines(W[, j])
  }
}

rdir <- function(nsamples, alpha){
  # Generate random samples from a Dirichlet distribution
  if(any(alpha == 0)){
    stop("Entries of alpha need to be positive")
  }
  n <- length(alpha)
  out <- matrix(rgamma(nsamples * n, alpha, 1), ncol = n, nrow = nsamples, byrow = TRUE) + 1e-9
  #res <- pmax(t(apply(out, 1, function(x) x/sum(x))), 1e-9)
  res <- t(apply(out, 1, function(x) x/sum(x)))
  return(res)
}

ddir <- function(x, alpha){
  # log of the density of a Dirichlet distribution
  lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) * log(x))
}

log_posterior <- function(w, phi, N, logSumY, X, a, a0, b0) {
  I <- nrow(X)
  # Find dot product
  dot_w_X <- X %*% w
  # log density
  log_target <- N * lgamma(phi) - N * sum(lgamma(phi * dot_w_X)) + sum((phi * dot_w_X - 1) * logSumY) +
    sum((a - 1) * log(w)) + (a0 - 1) * log(phi) - b0 * phi
  return(log_target)
} 

# sum(apply(y, 1, function(x) ddir(x, alpha = 800 * X %*% w)))
# 
# log_posterior(w, phi, N, logSumY, X, a, a0, b0)
# log_posterior(w_new, phi, N, logSumY, X, a, a0, b0)


sample_w <- function(w, phi, N, logSumY, X, a, a0, b0, V, alpha_proposal, method = "max_density", useJulia = TRUE){
  if(method == "max_density"){
    # Propose the maximum density approach
    if(useJulia) {
      alpha_new <- julia_call("max_density_for_dirichlet", w, V)
    } else {
      alpha_new <- dens.max.dir(center = w, V = V)  
    }
  } else if (method == "fixed") {
    alpha_new <- alpha_proposal * w
  }
  
  # Propose a new value
  w_new <- c(rdir(1, alpha = alpha_new))
  # Calculate log target under the two values
  lp_old <- log_posterior(w, phi, N, logSumY, X, a, a0, b0)
  lp_new <- log_posterior(w_new, phi, N, logSumY, X, a, a0, b0)
  
  # Calculate MH ratio correction for symmetry
  if(method == "max_density") {
    if(useJulia) {
      alpha_old <- julia_call("max_density_for_dirichlet", w_new, V)
    } else {
        alpha_old <- dens.max.dir(center = w_new, V = V)  
      }
    sym_correction <- ddir(w, alpha_new) - ddir(w_new, alpha_old)
    
  } else if (method == "fixed") {
    sym_correction <- ddir(w, alpha_proposal * w_new) - ddir(w_new, alpha_proposal * w)
  }
  
  # Calculate metropolis ratio
  log_MH_ratio_w <- lp_new - lp_old + sym_correction
  
  # Return the true value
  if(log(runif(1)) <= log_MH_ratio_w){
    w <- w_new
  }
  
  return(w)
}

sample_phi <- function(phi, w, N, logSumY, X, a, a0, b0, maxit = 10, tol = 1e-8) {
  dot_X_w <- X %*% w
  # Use the newton method to find the optimal proposal
  x <- phi
  doubleSum_y <- sum(dot_X_w * logSumY)
  A <- a0 + N/2 #a0 - x * N * trigamma(x) + x * N * sum(dot_X_w^2 * trigamma(x * dot_X_w))
  B <- b0 #b0 + (A - a0) / x - N * digamma(x) + N * sum(dot_X_w * digamma(x * dot_X_w)) - doubleSum_y
  # Get A and B for the proposal
  for(m in 1:maxit){
    x <- A/B
    A <- a0 - x^2 * N * trigamma(x) + x^2 * N * sum(dot_X_w^2 * trigamma(x * dot_X_w))
    B <- b0 + (A - a0) / x - N * digamma(x) + N * sum(dot_X_w * digamma(x * dot_X_w)) - doubleSum_y
    if (abs(x/(A/B) - 1) < tol){
      break
    }
  }
  # Propose from the gamma distribution
  phi_new <- rgamma(1, A, B)
  # Calculate the value of the log-posterior at both points
  lp_old <- log_posterior(w, phi, N, logSumY, X, a, a0, b0)
  lp_new <- log_posterior(w, phi_new, N, logSumY, X, a, a0, b0)
  sym_correction <- dgamma(x = phi, A, B, log = TRUE) - dgamma(x = phi_new, A, B, log = TRUE)
  log_MH_ratio <- lp_new - lp_old + sym_correction
  
  if(log(runif(1)) <= log_MH_ratio){
    phi <- phi_new
  }
  
  return(phi)
}
  
  
DirichletReg <- function(y, X, a = 1, phi = 1000, fixed_phi = FALSE, 
                         a0 = 1, b0 = 1, nsamples = 2000, burnin = 2000, 
                         method = "max_density", useJulia = TRUE,
                         alpha_proposal = 1e5, V = 6e-7, init_nnls = TRUE){
  if(!method %in% c("max_density", "fixed")){
    stop("method must either be either 'max_density' or 'fixed'")
  }
  # Data dimension
  N <- nrow(y); p <- ncol(X)
  # Useful quantities
  logSumY <- colSums(log(y))
  
  # Storing values
  W <- matrix(NA, nrow = nsamples, ncol = p)
  colnames(W) <- colnames(X)
  PHI <- rep(NA, nsamples)
  ACCEPT <- matrix(NA, nrow = nsamples, ncol = 2)
  colnames(ACCEPT) <- c("w", "phi")
  
  # Initial value for the sampler at the normalized nnls solution. 
  if(init_nnls){
    w <- nnls(X, colMeans(y))$x + 1e-8
    w <- pmax(w / sum(w), 1e-9)
  } else {
    w <- rep(1/p, p)#c(rdir(1, phi * rep(1/p, p)))
  }
  
  # Run the sampler
  pb <- txtProgressBar(style=3)
  for(r in 1:(nsamples + burnin)) {
    setTxtProgressBar(pb, r/(nsamples + burnin))
    #------------------------------------------------ Step 1. Update w
    w <- sample_w(w, phi, N = N, logSumY = logSumY, X = X, a = a, a0 = a0, useJulia = useJulia,
                  b0 = b0, V = V, alpha_proposal = alpha_proposal, method = method)
    #------------------------------------------------ Step 2. Update phi
    if(!fixed_phi){
      phi <- sample_phi(phi = phi, w = w, N = N, logSumY = logSumY, X = X, a = a, a0 = a0, b0 = b0)
    }
    #------------------------------------------------ Step 3. Save the output
    if(r > burnin){
      W[r - burnin, ] <- w
      PHI[r - burnin] <- phi
    }
  }
  close(pb)
  return(list(w = W, phi = PHI, accept = 1 - coda::rejectionRate(coda::as.mcmc(cbind("w" = W[, 1], "phi" = PHI)))))
}  



# ---------------- test the function
# Load the data
cosmic_3.0 <- read.COSMIC.sigs(version = "v3")

X <- pmax(as.matrix(cosmic_3.0), 1e-13)

# Data parameters
I <- nrow(X)
p <- ncol(X)

# Generate the data from above (mix-up between SBS2 and SBS13)
prob <- .3
sig_true <- prob * X[, 1] + (1-prob) * X[, 5]
plot.SBS.signature(sig_true)


set.seed(10)
N <- 200
y <- rdirichlet(n = N, alpha = 1000 * sig_true) + 1e-8
y <- t(apply(y, 1, function(x) x/sum(x)))
colnames(y) <- rownames(X)

a <- rep(1, p)

#---- Part 1 - stan
out_stan <- DirichletReg_stan(y, X, a)
plot(out_stan$phi, type = "l")
plot(colMeans(out_stan$probs))

# Part 2 - fixed case
out_fixed <- DirichletReg(y, X, a, nsamples = 1000, burnin = 0, method = "fixed", 
                          alpha_proposal = 1000, phi = 100, 
                          fixed_phi = FALSE, init_nnls = FALSE)
plot(density(out_fixed$phi))
lines(density(out_stan$phi), col = "red")
out_fixed$accept
head(sort(colMeans(tail(out_fixed$w), 20), decreasing = TRUE))
plot.trace(out_fixed$w)


par(mfrow = c(1,2))
plot(out_fixed$w[, 1], type = "l"); acf(out_fixed$w[, 1]) 
plot(out_fixed$w[, 5], type = "l"); acf(out_fixed$w[, 5]) 
par(mfrow = c(1,1))

# Part 3 - max density case
out_max <- DirichletReg(y, X, a, nsamples = 1000, burnin = 0, phi = 100,
                        method = "max_density", V = 6e-7, init_nnls = FALSE, fixed_phi = FALSE)
plot(density(out_max$phi))
lines(density(out_stan$phi), col = "red")
plot(colMeans(out_max$w))
plot.trace(out_fixed$w)
out_max$accept

par(mfrow = c(1,2))
plot(out_max$w[, 1], type = "l"); acf(out_max$w[, 1]) 
plot(out_max$w[, 5], type = "l"); acf(out_max$w[, 5]) 
par(mfrow = c(1,1))


# Compare against stan
fit.CompNMF <- readRDS(file = "data/fit21Breast.rds.gzip")
sig.MCMC.chain <- fit.CompNMF$mcmc_out[[1]]$Signatures
plot.SBS.signature(fit.CompNMF$Signatures)

# Select the signature
y <- sig.MCMC.chain[, , 10]
colnames(y) <- rownames(X)
a = rep(0.5, ncol(X))

# Part 1 - Stan
out_stan <- DirichletReg_stan(y, X[, 1:20], a[1:20])
plot(out_stan$phi, type = "l")
plot(colMeans(out_stan$probs))
res_stan <- colMeans(out_stan$probs)
names(res_stan) <- colnames(X)[1:20]
round(res_stan, 6)
plot.trace(out_stan$probs)


sum(dens.max.dir(out_max$w[1000, ], V = 1e-8))

# Part 2 - mean method
out_fixed <- DirichletReg(y = y, X = X[, 1:20], a = a[1:20], nsamples = 10000, burnin = 0, method = "fixed", 
                          alpha_proposal = 1e5, phi = 1, 
                          fixed_phi = FALSE, init_nnls = TRUE)
plot(density(out_fixed$phi))
lines(density(out_stan$phi), col = "red")
plot.trace(out_fixed$w)
out_fixed$accept
plot(colMeans(out_fixed$w))
plot.trace(out_fixed$w)

cbind(colMeans(out_fixed$w), colMeans(out_stan$probs), abs(colMeans(out_fixed$w) - colMeans(out_stan$probs)))

plot(colMeans(tcrossprod(out_fixed$w, X)), colMeans(y))
plot(colMeans(tcrossprod(out_stan$probs, X)), colMeans(y))

sqrt(mean((colMeans(tcrossprod(out_fixed$w, X)) - colMeans(y))^2))
sqrt(mean((colMeans(tcrossprod(out_stan$probs, X)) - colMeans(y))^2))


# Part 3 - max density case
out_max <- DirichletReg(y = y, X = X[, 1:20], a = a[1:20], nsamples = 5000, burnin = 0, method = "max_density", 
                        alpha_proposal = 1e5, phi = 1, V = 1e-9, useJulia = FALSE,
                        fixed_phi = FALSE, init_nnls = TRUE)
plot(density(out_max$phi))
lines(density(out_fixed$phi), col = "red")
plot(colMeans(out_max$w))
out_max$accept
plot.trace(out_max$w)

plot(colMeans(tcrossprod(out_fixed$w, X[, 1:20])), colMeans(y))
plot(colMeans(tcrossprod(out_max$w, X[, 1:20])), colMeans(y))
plot(colMeans(tcrossprod(out_stan$probs, X[, 1:20])), colMeans(y))

sqrt(mean((colMeans(tcrossprod(out_fixed$w, X[, 1:20])) - colMeans(y))^2))
sqrt(mean((colMeans(tcrossprod(out_stan$probs, X[, 1:20])) - colMeans(y))^2))
sqrt(mean((colMeans(tcrossprod(out_max$w, X[, 1:20])) - colMeans(y))^2))

plot(colMeans(out_max$w), colMeans(out_stan$probs))
plot(colMeans(out_fixed$w), colMeans(out_stan$probs))


plot(colMeans(out_max$w), colMeans(out_stan$probs))
abline(a = 0, b = 1)



out_fixed$w 

plot(out_max$w[1000, ]); points(c(rdir(1, alpha = c(dens.max.dir(out_max$w[1000, ], V = 6e-8)))), col = "red")


# Check results
res <- extract(fit)
round(colMeans(res$probs), 5)[c(1, 5)]

plot(res$phi, type = "l")






res <- DirichletReg(y, X, nsamples = 10000, burnin = 10, method = "fixed", alpha_proposal = 1e6)

res$accept


set.seed(10)
out1 <- DirichletReg(y, X, a = 1, a0 = 1, b0 = 1,
                    nsamples = 2000, 
                    burnin = 2000, 
                    method = "fixed", fixed_phi = TRUE, phi = 800, 
                    V = 1e-9, alpha_proposal = 1e5)

colMeans(out1$w)[c(1,5)]
colMeans(out1$accept)
plot(out1$phi, type = "l")
plot(out1$w[, 1], type = "l")
plot(out1$w[, 5], type = "l")
head(sort(colMeans(out1$w), decreasing = TRUE))
coda::effectiveSize(out1$w)


out2 <- DirichletReg(y, X, a = 1, a0 = 1, b0 = 1,
                    nsamples = 1000, 
                    burnin = 10000, 
                    method = "max_density", fixed_phi = FALSE, phi = 800, 
                    V = 1e-8, alpha_proposal = 1e5)

colMeans(out2$w)[c(1,5)]
colMeans(out2$accept)
plot(out2$phi, type = "l")
plot(out2$w[, 1], type = "l")
plot(out2$w[, 5], type = "l")
head(sort(colMeans(out2$w), decreasing = TRUE))
coda::effectiveSize(out2$w)


##############
# 21 breast cancer application

# I estimate my model on the 21 breast cancer data, with no informative prior
# it will mix-up SBS2 and SBS13. We will use this method to de-couple it.

mut <- t(read.table(system.file("extdata","21_breast_cancers.mutations.txt", package="signeR"), header=TRUE, check.names=FALSE))
transformed_names <- apply(stringr::str_split(rownames(mut), ":", simplify = TRUE), 1, function(x) {
  y <- str_split(x[2], "", simplify = TRUE)
  y[2] <- paste0("[", x[1],"]")
  paste0(y, collapse = "")
})

MutMat <- mut[order(transformed_names), ]
rownames(MutMat) <- sort(transformed_names)
# Load the data
cosmic_3.0 <- read.COSMIC.sigs(version = "v3")
X <- as.matrix(cosmic_3.0)
MutMat <- MutMat[rownames(X), ]

# Estimate the signatures via CompressiveNMF
set.seed(42)
fit.nmf <- CompressiveNMF(X = MutMat, K = 20, nsamples = 1000, burnin = 2500, nchains = 1, seed = 42)
saveRDS(fit.nmf, file = "data/fit21Breast.rds.gzip", compress = "gzip")

y <- fit.nmf$mcmc_out[[1]]$Signatures[, , 10]

rownames(fit.nmf$Signatures) <- rownames(MutMat)

plot.SBS.signature(fit.nmf$Signatures)

# We see that the last signature has been mixed-up
sig1 <- fit.nmf$mcmc_out[[1]]$Signatures[, , 10]
sig2 <- fit.nmf$mcmc_out[[1]]$Signatures[, , 6]
colnames(sig1) <- rownames(X)
# Let's start Signature 1
plot.SBS.signature(fit.nmf$Signatures[, 4])
plot.SBS.signature(colMeans(sig1))

#-------------- STAN
fit.stan_sig1 <- DirichletReg_stan(sig1, X = X, a = rep(1, ncol(X)))
colnames(fit.stan_sig1$probs) <- colnames(X)
head(sort(colMeans(fit.stan_sig1$probs), decreasing = TRUE))
plot(fit.stan_sig1$phi, type = "l")

plot(fit.stan_sig1$probs[, 18], type = "l")
dim(fit.stan_sig1$probs)
out <- X %*% t(fit.stan_sig1$probs)
plot.SBS.signature(rowMeans(out))

sol_nnls <- nnls(A = X, b = colMeans(sig1))$x
names(sol_nnls) <- colnames(X)
sol_nnls <- sol_nnls/sum(sol_nnls)
plot(sol_nnls, colMeans(fit.stan_sig1$probs))
abline(a = 0, b = 1)

plot(rowMeans(out), X %*% sol_nnls)
abline(a = 0, b = 1)



#-------------- Simple MH proposal
p <- .3
sig_true <- p * X[, 1] + (1-p) * X[, 5]
plot.SBS.signature(sig_true)


set.seed(10)
S <- 100
y <- rdirichlet(n = S, alpha = 1000 * sig_true) + 1e-8
y <- t(apply(y, 1, function(x) x/sum(x)))

fit.MH_sig1 <- DirichletReg(y = y, X = X, a = 1, a0 = 1, b0 = 1,
                      nsamples = 5000, 
                      burnin = 10000, 
                      method = "fixed", fixed_phi = TRUE, phi = 400, 
                      V = 1e-9, alpha_proposal = 1e5)
round(head(sort(colMeans(fit.MH_sig1$w), decreasing = TRUE)), 5)

colMeans(fit.MH_sig1$accept)

out <- X %*% t(fit.MH_sig1$w)
plot.SBS.signature(X[, "SBS12"])

par(mfrow = c(1,1))
plot(fit.MH_sig1$w[, 1], type = "l")
plot(fit.MH_sig1$w[, 5], type = "l")
plot(stats::acf(as.matrix(fit.MH_sig1$w[, 1]), plot = FALSE))
plot(stats::acf(as.matrix(fit.MH_sig1$w[, 2]), plot = FALSE))



#-------------- Maximum density proposal
set.seed(10)
fit.MaxDens_sig1 <- DirichletReg(y, X, a = 1, a0 = 1, b0 = 1,
                            nsamples = 10000, 
                            burnin = 20000, 
                            method = "max_density", fixed_phi = TRUE, phi = 819, 
                            V = 2e-9, alpha_proposal = 1e5)
head(sort(colMeans(fit.MaxDens_sig1$w), decreasing = TRUE))
plot(fit.MaxDens_sig1$phi, type = "l")

par(mfrow = c(2,2))
plot(fit.MaxDens_sig1$w[, 2], type = "l", main = "SBS2, MaxDens, V = 1e-9")
plot(fit.MaxDens_sig1$w[, 7], type = "l", main = "SBS7a,  MaxDens, V = 1e-9")
plot(stats::acf(as.matrix(fit.MaxDens_sig1$w[, 2]), plot = FALSE), main = "SBS2, MaxDens, V = 1e-9")
plot(stats::acf(as.matrix(fit.MaxDens_sig1$w[, 7]), plot = FALSE), main = "SBS7a, MaxDens, V = 1e-9")

colMeans(fit.MaxDens_sig1$accept)


#--------------- Comparison with NNLS
nnls_out <- nnls(A = X, b = colMeans(sig1))$x
nnls_out <- nnls_out/sum(nnls_out)
par(mfrow = c(1,1), pty="s")
plot(nnls_out[-c(1,5)], colMeans(fit.MaxDens_sig1$w)[-c(1,5)], ylim =c(0, 1), xlim = c(0,1), asp=1, 
     ylab = "DirichletReg", xlab = "nnls")
points(nnls_out[1], colMeans(fit.MaxDens_sig1$w)[1], col = "red")
points(nnls_out[5], colMeans(fit.MaxDens_sig1$w)[5], col = "blue")
abline(a = 0, b = 1)



coda::autocorr.plot(coda::as.mcmc(fit.MH_sig1$w[, c(1,5)]))
stats::acf(as.matrix(fit.MH_sig1$w[, 1]))
stats::acf(as.matrix(fit.MaxDens_sig1$w[, 1]))

coda::autocorr.plot(coda::as.mcmc(fit.MaxDens_sig1$w[, 1]))
coda::autocorr.plot(coda::as.mcmc(fit.MaxDens_sig1$w[, 5]))



#### Let's use the model


# Load the data
cosmic_3.0 <- read.COSMIC.sigs(version = "v3")

X <- as.matrix(cosmic_3.0)

# Data parameters
I <- nrow(X)
J <- ncol(X)

# Generate the data from above (mix-up between SBS2 and SBS13)
p <- 0
sig_true <- p * X[, 1] + (1-p) * X[, 5]
plot.SBS.signature(sig_true)

set.seed(10)
S <- 50
y <- rdirichlet(n = S, alpha = 1000 * sig_true) + 1e-8
y <- t(apply(y, 1, function(x) x/sum(x)))

y <- sig1[1:S, ]
colnames(y) <- rownames(X)
plot.SBS.signature(colMeans(y))
a <- rep(0.01, J)
#a[5] <- 0.01

# Compile the model
stan_data <- list(S = S, I = I, J = J, X = X, y = y, a = a)
model_file <- "stan/simplex_regression.stan"  
init_fun <- list(list(probs = rep(1/J, J), phi = 100))
fit <- stan(file = model_file, data = stan_data, iter = 2000, chains = 1, 
            init = init_fun)
res_stan <- rstan::extract(fit)
w <- round(colMeans(res_stan$probs), 5)
names(w) <- colnames(X)
sort(w)


## Illustation
apply(X, 2, function(x) - sum(x * log(x)))


cut(apply(X, 2, function(x) - sum(x * log(x))), breaks = 5)

out <- list()
for(i in 1:ncol(X)){
  x <- X[, i]
  alpha <- c(1:50) * 10
  res <- sapply(alpha, function(a) apply(rdir(nsamples = 200, alpha = a*x), 1, function(y) lsa::cosine(x, y)))
  out[[i]] <- res
}
names(out) <- colnames(X)
out_means <- as.data.frame(lapply(out, colMeans))

entropies <- apply(X, 2, function(x) - sum(x * log(x)))
values <- cut(entropies, breaks = 4, labels = FALSE)
cut(entropies, breaks = 4)

flatness <- apply(X, 2, function(x) 1/(96 * sum(x^2)))
values <- cut(flatness, breaks = 4, labels = FALSE)
cut(flatness, breaks = 4)

plot(alpha, out_means[, 1], ylim = c(0.35, 1), type = "l",col = values[1], ylab="cosine similarity with mean")
for(i in 2:20){
  lines(alpha, out_means[, i], col = values[i])
}

out <- list()
V <- 1/c(10, 100, 1000, 10000, 100000)
for(i in c(1:20)){
  x <- X[, i]
  res <- sapply(V, function(a) apply(rdir(nsamples = 1000, alpha = dens.max.dir(x, a)), 1, function(y) lsa::cosine(x, y)))
  out[[i]] <- res
}
#names(out) <- colnames(X)
out_means <- as.data.frame(lapply(out, colMeans))

plot(log(V), out_means[, 1], type = "l",col = values[1], ylab="cosine similarity with mean")
for(i in 2:ncol(out_means)){
  lines(log(V), out_means[, i], col = values[i])
}





sort(colMeans(out_means))











