# Load packages
library(tidyverse)
#library(JuliaCall)
library(scales)

#julia <- julia_setup(installJulia = TRUE)
  
# Source useful files
source("R/generate_synthetic_mutation_counts.R")
source("R/plot_signature.R")

# Useful function
rdir <- function(nsamples, alpha){
  # Generate random samples from a Dirichlet distribution
  if(any(alpha == 0)){
    stop("Entries of alpha need to be positive")
  }
  n <- length(alpha)
  out <- matrix(rgamma(nsamples * n, alpha, 1), ncol = n, nrow = nsamples, byrow = TRUE) + 1e-13
  #res <- pmax(t(apply(out, 1, function(x) x/sum(x))), 1e-9)
  res <- t(apply(out, 1, function(x) x/sum(x)))
  return(res)
}

# Load data
cosmic_3.4 <- read.COSMIC.sigs(version = "v3.4")
#X <- pmax(as.matrix(cosmic_3.4), 1e-10)
X <- apply(cosmic_3.4, 2, function(x) x/sum(x))
write.csv(X, file = "data/cosmic3.4.csv")

# Fixe the number of samples to run the analysis
nsamples <- 2000

#------------------------------------------------------ Mean method
set.seed(42)
res_mean <- vector(mode = "list", length = ncol(X))
seq_alphas <- pracma::logseq(x1 = 0.001, x2 = 5000, n = 500)
for(j in 1:ncol(X)){
  print(paste0("Mean method - ", j))
  x <- X[, j]
  alpha_new <- sapply(seq_alphas, function(y) x * y)#sapply(1:300, function(y) x * y )
  cosine_samples <- apply(alpha_new, 2, function(y) 
    apply(rdir(nsamples, y), 1, function(z) 1 - lsa::cosine(x, z)))
  res_mean[[j]] <- cosine_samples
}
names(res_mean) <- colnames(X)

#------------------------------------------------------ Maximum density
#alpha_max_density <- julia_eval('using MAT; matread("data/alpha_cosmic_cosine.mat")["results_cosine"]')
#saveRDS(alpha_max_density, file = "data/alpha_max_density_COSMICv3.4.rds.gzip", compress = "gzip")
alpha_max_density <- readRDS("data/alpha_max_density_COSMICv3.4.rds.gzip")
dimnames(alpha_max_density) <- list(rownames(X), 1:dim(alpha_max_density)[2], colnames(X))

set.seed(42)
res_max_density <- vector(mode = "list", length = ncol(X))
for(j in 1:ncol(X)){
  print(paste0("Density max - ", j))
  x <- X[, j]
  cosine_samples <- apply(alpha_max_density[, , j], 2, function(y) 
    apply(rdir(nsamples, y), 1, function(z) 1 - lsa::cosine(x, z)))
  res_max_density[[j]] <- cosine_samples
}
names(res_max_density) <- colnames(X)

# Compute the average cosine similarity for the two methods
avg_mean <- as.matrix(data.frame(lapply(res_mean, function(x) 
  apply(x, 2, function(y) c("mean" = mean(y))))))
avg_max_density <- as.matrix(data.frame(lapply(res_max_density, function(x) 
  apply(x, 2, function(y) c("mean" = mean(y))))))

saveRDS(avg_mean, file = "output/Cosines_COSMIC_mean.rds")
saveRDS(avg_max_density, file = "output/Cosines_COSMIC_max_density.rds")

# Find the value for alpha and theta to obtain a target_avg_cos on average 
# across all signatures
target_avg_cos <- sort(seq(0.0, 0.3, 0.025), decreasing = TRUE)
target_avg_cos[length(target_avg_cos)] <- 0.01

#------------------------------ Mean Method
alpha_vals_target_mean <- sapply(target_avg_cos, 
                                 function(x) 
                                   which(!rowMeans(avg_mean) > x)[1])
out_target_mean <- NULL
for(i in 1:length(target_avg_cos)){
  out_target_mean <- rbind(out_target_mean, 
                           avg_mean[alpha_vals_target_mean[i], ])
}

#------------------------------ Density maximization
alpha_vals_target_max_dens <- unname(sapply(target_avg_cos, 
                                            function(x) 
                                              which(rowMeans(avg_max_density) > x)[1]))
out_target_max_dens <- NULL
for(i in 1:length(target_avg_cos)){
  out_target_max_dens <- rbind(out_target_max_dens,
                               avg_max_density[alpha_vals_target_max_dens[i], ])
}


# Make the plot 
pdf(file = "plots/MutationalSignatures.pdf",  width = 5.84, height = 5.11)

par(mfrow = c(1,1))
plot(target_avg_cos, 
     apply(out_target_mean, 1, function(x) quantile(x, 0.75, na.rm = TRUE)), 
     col = 2, type = "l", lwd = 1.7, xlim = c(0.01, 0.3), 
     xlab = "Mean cosine error", 
     ylab = "Average mean cosine error")
lines(target_avg_cos, 
      apply(out_target_mean, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), 
      col = 2, lwd = 1.7)
lines(target_avg_cos, 
     apply(out_target_max_dens, 1, function(x) quantile(x, 0.75, na.rm = TRUE)), 
     col = 4, lwd = 1.7)
lines(target_avg_cos, 
      apply(out_target_max_dens, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), 
      col = 4, lwd = 1.7)
abline(a = 0, b = 1, col = "darkgrey", lty = "dashed", lwd = 1.2)
legend(0.01, 0.4, bty = "n", x.intersp = 0.75, y.intersp = 1.5, lwd = 1.5, seg.len = 1.25, cex = 1, 
       legend = c("mean method", "maximum density"), col = c(2, 4))

dev.off()

# 
# 
# par(mfrow = c(1,1))
# plot(target_avg_cos, rowMeans(res1), type = "l", ylim = c(0, 0.43), xlab = "Target avg. error", xlim = c(0.01, 0.3), ylab = "Across signatures avg. error")
# lines(target_avg_cos, apply(res1, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), lty = "dashed", )
# lines(target_avg_cos, apply(res1, 1, function(x) quantile(x, 0.75, na.rm = TRUE)), lty = "dashed")
# 
# 
# 
# #------------------------------------------------------ Fixed concentration
# set.seed(42)
# res_fixed <- vector(mode = "list", length = ncol(X))
# seq_alphas <- pracma::logseq(x1 = 0.001, x2 = 5000, n = 500)
# for(j in 1:ncol(X)){
#   print(paste0("Fixed - ", j))
#   x <- X[, j]
#   alpha_new <- sapply(seq_alphas, function(y) x * y)#sapply(1:300, function(y) x * y )
#   cosine_samples <- apply(alpha_new, 2, function(y) apply(rdir(nsamples, y), 1, function(z) 1 - lsa::cosine(x, z)))
#   res_fixed[[j]] <- cosine_samples
# }
# names(res_fixed) <- colnames(X)
# 
# 
# #------------------------------------------------------ Constrained KL
# alpha_kl <- julia_eval('using MAT; matread("data/alpha_cosmic_KL.mat")["results_kl"]')
# dimnames(alpha_kl) <- list(rownames(X), 1:dim(alpha_kl)[2], colnames(X))
# # Run the method
# set.seed(42)
# res_kl <- vector(mode = "list", length = ncol(X))
# for(j in 1:ncol(X)){
#   print(paste0("KL - ", j))
#   x <- X[, j]
#   res <- matrix(NA, nrow = nsamples, ncol = dim(alpha_kl)[2])
#   alpha_new <- alpha_kl[, ,j]
#   cosine_samples <- apply(alpha_new, 2, function(y) apply(rdir(nsamples, y), 1, function(z) 1 - lsa::cosine(x, z)))
#   res_kl[[j]] <- cosine_samples
# }
# names(res_kl) <- colnames(X)
# 
# #------------------------------------------------------ Constrained concentration
# alpha_conc <- julia_eval('using MAT; matread("data/alpha_cosmic_concentration.mat")["results_conc"]')
# dimnames(alpha_conc) <- list(rownames(X),  1:dim(alpha_conc)[2], colnames(X))
# 
# set.seed(42)
# res_conc <- vector(mode = "list", length = ncol(X))
# for(j in 1:ncol(X)){
#   print(paste0("Concentration - ", j))
#   x <- X[, j]
#   cosine_samples <- apply(alpha_conc[, , j], 2, function(y) apply(rdir(nsamples, y), 1, function(z) 1 - lsa::cosine(x, z)))
#   res_conc[[j]] <- cosine_samples
# }
# names(res_conc) <- colnames(X)
# 
# 
# #------------------------------------------------------ Constrained cosine
# alpha_cosine <- julia_eval('using MAT; matread("data/alpha_cosmic_cosine.mat")["results_cosine"]')
# dimnames(alpha_cosine) <- list(rownames(X),  1:dim(alpha_cosine)[2], colnames(X))
# 
# set.seed(42)
# res_cosine <- vector(mode = "list", length = ncol(X))
# for(j in 1:ncol(X)){
#   print(paste0("Cosine - ", j))
#   x <- X[, j]
#   cosine_samples <- apply(alpha_cosine[, , j], 2, function(y) apply(rdir(nsamples, y), 1, function(z) 1 - lsa::cosine(x, z)))
#   res_cosine[[j]] <- cosine_samples
# }
# names(res_cosine) <- colnames(X)
# 
# 
# 
# #------ Analyze the output
# avg_cos_kl <- lapply(res_kl, function(x) apply(x, 2, function(y) c("mean" = mean(y),
#                                                                    "stdev" = sd(y),
#                                                                    "lowCI" = quantile(y, 0.25), 
#                                                                    "highCI" = quantile(y, 0.75))))
# 
# avg_cos_conc <- lapply(res_conc, function(x) apply(x, 2, function(y) c("mean" = mean(y),
#                                                                        "stdev" = sd(y),
#                                                                        "lowCI" = quantile(y, 0.25), 
#                                                                        "highCI" = quantile(y, 0.75))))
# 
# avg_cos_cosine <- lapply(res_cosine, function(x) apply(x, 2, function(y) c("mean" = mean(y),
#                                                                        "stdev" = sd(y),
#                                                                        "lowCI" = quantile(y, 0.25), 
#                                                                        "highCI" = quantile(y, 0.75))))
# 
# avg_cos_fixed <- lapply(res_fixed, function(x) apply(x, 2, function(y) c("mean" = mean(y),
#                                                                          "stdev" = sd(y),
#                                                                          "lowCI" = quantile(y, 0.25), 
#                                                                          "highCI" = quantile(y, 0.75))))
# 
# 
# mean_kl <- as.matrix(data.frame(lapply(avg_cos_kl, function(x) x[1,])))
# mean_conc <- as.matrix(data.frame(lapply(avg_cos_conc, function(x) x[1,])))
# mean_cosine <- as.matrix(data.frame(lapply(avg_cos_cosine, function(x) x[1,])))
# mean_fixed <- as.matrix(data.frame(lapply(avg_cos_fixed, function(x) x[1,])))
# 
# par(mfrow = c(1, 4))
# # Mean
# plot(rowMeans(mean_fixed), type="l", main = "fixed", ylab = "avg. cosine error", xlab = "alpha", ylim = c(0, 1))
# lines(apply(mean_fixed,1, function(x) quantile(x, 0.25) ), col = "red", lty = "dashed")
# lines(apply(mean_fixed,1, function(x) quantile(x, 0.75) ), col = "red", lty = "dashed")
# 
# plot(rowMeans(mean_conc), type="l", main = "concentration", ylab = "avg. cosine error", xlab = "theta", ylim = c(0, 1))
# lines(apply(mean_conc,1, function(x) quantile(x, 0.25) ), col = "red", lty = "dashed")
# lines(apply(mean_conc,1, function(x) quantile(x, 0.75) ), col = "red", lty = "dashed")
# 
# plot(rowMeans(mean_kl), type="l", main = "KL", ylab = "avg. cosine error", xlab = "theta", ylim = c(0, 1))
# lines(apply(mean_kl,1, function(x) quantile(x, 0.25) ), col = "red", lty = "dashed")
# lines(apply(mean_kl,1, function(x) quantile(x, 0.75) ), col = "red", lty = "dashed")
# 
# theta_range <- seq(0.001, 1.5, length.out = 300)
# plot(theta_range, rowMeans(mean_cosine), type="l", main = "cosine", ylab = "avg. cosine error", xlab = "theta", ylim = c(0, 1))
# lines(theta_range, apply(mean_cosine,1, function(x) quantile(x, 0.25) ), col = "red", lty = "dashed")
# lines(theta_range, apply(mean_cosine,1, function(x) quantile(x, 0.75) ), col = "red", lty = "dashed")
# 
# 
# par(mfrow = c(1, 1), pty = "s")
# plot(theta_range, rowMeans(mean_cosine), type="l", main = "cosine", ylab = "avg. cosine error", xlab = "theta", ylim = c(0, 1), xlim = c(0, 1), asp = 1)
# lines(theta_range, apply(mean_cosine,1, function(x) quantile(x, 0.25) ), col = "red", lty = "dashed")
# lines(theta_range, apply(mean_cosine,1, function(x) quantile(x, 0.75) ), col = "red", lty = "dashed")
# 
# par(mfrow = c(1, 1), pty = "s")
# plot(theta_range, rowMeans(mean_cosine), type="l", main = "cosine", ylab = "avg. cosine error", xlab = "theta", ylim = c(0, .1), xlim = c(0, .1), asp = 1)
# lines(theta_range, apply(mean_cosine,1, function(x) quantile(x, 0.25) ), col = "red", lty = "dashed")
# lines(theta_range, apply(mean_cosine,1, function(x) quantile(x, 0.75) ), col = "red", lty = "dashed")
# abline(a = 0, b = 1, col = "grey")
# 
# 
# # Invert here the functions and produce a plot for all four methods
# par(mfrow = c(1,1))
# target_avg_cos <- sort(seq(0.0, 0.3, 0.025), decreasing = TRUE)
# target_avg_cos[length(target_avg_cos)] <- 0.01
# 
# cx <- rowMeans(mean_fixed)
# alpha_vals_target_fixed <- sapply(target_avg_cos, function(x) which(!cx > x)[1])
# 
# res1 <- NULL
# for(i in 1:length(target_avg_cos)){
#   res1 <- rbind(res1, mean_fixed[alpha_vals_target_fixed[i], ])
# }
# 
# plot(target_avg_cos, rowMeans(res1), type = "l", xlim = c(0, 0.3), ylim = c(0, 0.44))
# lines(target_avg_cos, apply(res1, 1, function(x) quantile(x, 0.5, na.rm = TRUE)), col = "red", lty = "dashed")
# lines(target_avg_cos, apply(res1, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), col = "red", lty = "dashed")
# lines(target_avg_cos, apply(res1, 1, function(x) quantile(x, 0.85, na.rm = TRUE)), col = "red", lty = "dashed")
# 
# 
# cx <- rowMeans(mean_kl)
# alpha_vals_target_kl <- unname(sapply(target_avg_cos, function(x) which(!cx > x)[1]))
# 
# res2 <- NULL
# for(i in 1:length(target_avg_cos)){
#   res2 <- rbind(res2, mean_kl[alpha_vals_target_kl[i], ])
# }
# 
# 
# plot(target_avg_cos, rowMeans(res2), type = "l", ylim = c(0, 0.7))
# lines(target_avg_cos, apply(res2, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), col = "red", lty = "dashed")
# lines(target_avg_cos, apply(res2, 1, function(x) quantile(x, 0.85, na.rm = TRUE)), col = "red", lty = "dashed")
# 
# 
# cx <- rowMeans(mean_conc)
# alpha_vals_target_conc <- unname(sapply(target_avg_cos, function(x) which(!cx > x)[1]))
# 
# res3 <- NULL
# for(i in 1:length(target_avg_cos)){
#   res3 <- rbind(res3, mean_conc[alpha_vals_target_conc[i], ])
# }
# 
# 
# plot(target_avg_cos, rowMeans(res2), type = "l", ylim = c(0, 0.7))
# lines(target_avg_cos, apply(res2, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), col = "red", lty = "dashed")
# lines(target_avg_cos, apply(res2, 1, function(x) quantile(x, 0.85, na.rm = TRUE)), col = "red", lty = "dashed")
# 
# 
# 
# cx <- rowMeans(mean_cosine)
# alpha_vals_target_cosine <- unname(sapply(target_avg_cos, function(x) which(cx > x)[1]))
# 
# res4 <- NULL
# for(i in 1:length(target_avg_cos)){
#   res4 <- rbind(res4, mean_cosine[alpha_vals_target_cosine[i], ])
# }
# 
# 
# par(mfrow = c(1,1))
# plot(target_avg_cos, rowMeans(res1), type = "l", ylim = c(0, 0.43), xlab = "Target avg. error", xlim = c(0.01, 0.3), ylab = "Across signatures avg. error")
# lines(target_avg_cos, apply(res1, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), lty = "dashed", )
# lines(target_avg_cos, apply(res1, 1, function(x) quantile(x, 0.75, na.rm = TRUE)), lty = "dashed")
# 
# # lines(target_avg_cos, rowMeans(res2), col = "red")
# # lines(target_avg_cos, apply(res2, 1, function(x) quantile(x, 0.2, na.rm = TRUE)), col = "red", lty = "dashed")
# # lines(target_avg_cos, apply(res2, 1, function(x) quantile(x, 0.8, na.rm = TRUE)), col = "red", lty = "dashed")
# # 
# # lines(target_avg_cos, rowMeans(res3), col = "blue")
# # lines(target_avg_cos, apply(res3, 1, function(x) quantile(x, 0.2, na.rm = TRUE)), col = "blue", lty = "dashed")
# # lines(target_avg_cos, apply(res3, 1, function(x) quantile(x, 0.8, na.rm = TRUE)), col = "blue", lty = "dashed")
# 
# lines(target_avg_cos, rowMeans(res4), col = "forestgreen")
# lines(target_avg_cos, apply(res4, 1, function(x) quantile(x, 0.15, na.rm = TRUE)), col = "forestgreen", lty = "dashed")
# lines(target_avg_cos, apply(res4, 1, function(x) quantile(x, 0.75, na.rm = TRUE)), col = "forestgreen", lty = "dashed")
# 
# legend(x = "topleft", lty = 1, col= c("black","red", "blue", "forestgreen"),  legend=c("Fixed", "KL", "Concentration", "Cosine")) 
# 
