# This file runs a mutational signature analysis on synthetic data to test the 
# effectiveness of the optimally-tuned Dirichlet prior
library(CompressiveNMF)
library(tidyverse)
library(LaplacesDemon)
library(lsa)
library(foreach)
library(doParallel)

# Source useful functions
source("R/plot_signature.R")
source("R/generate_synthetic_mutation_counts.R")
source("R/initialize_prior.R")


#-------------------------------------------------------------- Useful functions
# dens.max.dir <- function (center, V = 0.01, renormalize = TRUE) {
#   if (round(sum(center), 5) != 1 || sum(center < 0) > 0)  {
#     stop("center is not a valid point in the simplex")
#   }
#   if (renormalize && sum(center == 1e-16) > 0)  {
#     center[center == 1e-16] <- 1e-16
#     center <- center / sum(center)
#   }
#   
#   i.star <- which.min(sapply(center, function (c) {min(c, 1 - c)}))[1]
#   lst <- dens.max.var(center[i.star], V)
#   conc <- lst$a + lst$b
#   
#   u <- numeric()
#   for (i in 1:length(center)) {
#     if (i == i.star)  {
#       u <- c(u, lst$a / conc)
#     } else  {
#       u <- c(u, dens.max.conc(center[i], conc)$a / conc)
#     }
#   }
#   return(u * conc)
# } 

# dens.max.var <- function (center, V = 0.01) {
#   roots <- uniroot.all(h.prime, lower = u.min(V), upper = u.max(V), center = center, V = V)
#   ind <- sapply(roots, prop.dens, loc = center, V = V) %>% which.max()
#   if (length(ind) > 1)  ind <- ind[1]
# 
#   u <- roots[ind]
#   alpha <- u * (1 - u) / V - 1
#   return(list(a = alpha * u, b = alpha * (1 - u)))
# } 

simulate.count.matrix <- function(loadings, sigs, overdispersion = 0){
  I <- nrow(sigs); J <- ncol(loadings)
  X <- matrix(rnbinom(I * J, size = 1/overdispersion, mu = c(sigs %*% loadings)), nrow = I, ncol = J)
  colnames(X) <- colnames(loadings)
  rownames(X) <- rownames(sigs)
  return(X)
}

tune_prior <- function(sigs, mode = "maxdens", V = 6e-7, alpha_fixed = 1000){
  if(mode == "maxdens"){
    S <- apply(sigs, 2, function(x) suppressWarnings(dens.max.dir(x, V = V)))
  } else if (mode == "fixed") {
    S <- alpha_fixed * sigs
  } else {
    stop("Incorrect mode specified")
  }
  rownames(S) <- rownames(sigs)
  return(S)
}
#--------------------------------------------------------------

# Load the files
path <- "data/SynSigGen_data/S.2.Rsq.0.6/"

path_counts <- file.path(path, "ground.truth.syn.catalog.csv")
path_loadings <- file.path(path, "ground.truth.syn.exposures.csv")
path_sigs <- file.path(path,"/ground.truth.syn.sigs.csv")

# open files
df_counts <- read.csv(path_counts)
df_sigs <- read.csv(path_sigs)
df_loadings <- read.csv(path_loadings)

# Make matrices
X <- as.matrix(df_counts[, -c(1,2)])
trinucl <- str_split(df_counts$Trinucleotide, "", simplify = TRUE)
trinucl[, 2] <- paste0("[", df_counts$Mutation.type, "]")
rownames(X) <- apply(trinucl, 1, function(x) paste0(x, collapse = ""))
# signatures
sigs <- apply(as.matrix(df_sigs[, -c(1,2)]), 2, function(x) x/sum(x))
rownames(sigs) <- apply(trinucl, 1, function(x) paste0(x, collapse = ""))
# loadings 
loadings <- as.matrix(df_loadings[, -1])
rownames(loadings) <- df_loadings[, 1]
  
# Reorder the columns and rows to have the same order across all files
loadings <- loadings[c("SBS1", "SBS5"), ]
sigs <- sigs[, c("SBS1", "SBS5")]

plot(log(loadings[1,]), log(loadings[2, ]))
plot(loadings[1,], loadings[2, ])

sigs_all <- read.COSMIC.sigs(version = "v3")
sigs_all <- as.matrix(apply(sigs_all, 2, function(x) x/sum(x)))
sigs_temp <- sigs_all#[, c("SBS1", "SBS5")]#[, c("SBS1", "SBS2", "SBS3", "SBS13", "SBS5")]
# Tune the infomative prior
S_maxdens <- tune_prior(sigs_temp, V = 6e-7, mode = "maxdens")
S_fixed <- tune_prior(sigs_temp, alpha_fixed = 242, mode = "fixed")

i <- 3
samples <- rdirichlet(2000, S_maxdens[, i])
CI <- t(apply(samples, 2, function(x) quantile(x, c(0.05, 0.95))))
p1 <- plot.SBS.signature(sigs_temp[, i], CI = CI) + ggtitle("Max_dens")

samples <- rdirichlet(2000, S_fixed[, i])
CI <- t(apply(samples, 2, function(x) quantile(x, c(0.05, 0.95))))
p2 <- plot.SBS.signature(sigs_temp[, i], CI = CI) + ggtitle("Fixed")
ggpubr::ggarrange(p1, p2, ncol = 1)
colSums(S_maxdens)

# Run the model
nsamples <- 1000
burnin <- 2000

I <- nrow(sigs)
J <- 50


set.seed(42)
id_sample <- sort(sample(1:ncol(X), J))
X_sim <- simulate.count.matrix(loadings[, id_sample], sigs, overdispersion = 0)
plot(X[, id_sample], X_sim)


out_noprior <- CompressiveNMF(X_sim, K = 20, 
                              nsamples = nsamples, 
                             burnin = burnin, 
                             alpha = 0.5)
CompressiveNMF:::print.CompressiveNMF(out_noprior)

rownames(out_noprior$Signatures) <- rownames(X)
plot.SBS.signature(out_noprior$Signatures)

out_maxdens <- CompressiveNMF(X_sim, K = 0, nsamples = nsamples, 
                              burnin = burnin, 
                              S = S_maxdens, 
                              betah = rep(1, ncol(S_maxdens)),
                              alpha = 0.5, swap_prior = FALSE)
CompressiveNMF:::print.CompressiveNMF(out_maxdens)

rownames(out_maxdens$Signatures) <- rownames(X)
plot.SBS.signature(out_maxdens$Signatures)


out_fixed <- CompressiveNMF(X_sim, K = 0, nsamples = nsamples, 
                              burnin = burnin, 
                              S = S_fixed, 
                              betah = rep(1, ncol(S_maxdens)),
                              alpha = 0.5, swap_prior = FALSE)
CompressiveNMF:::print.CompressiveNMF(out_fixed)

rownames(out_fixed$Signatures) <- rownames(X)
plot.SBS.signature(out_fixed$Signatures)


sqrt(mean(out_maxdens$Signatures - sigs)^2)
sqrt(mean(out_fixed$Signatures - sigs)^2)

