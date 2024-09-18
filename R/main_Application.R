# This file runs a mutational signature analysis on synthetic data to test the 
# effectiveness of the optimally-tuned Dirichlet prior

library(CompressiveNMF)
library(tidyverse)
library(LaplacesDemon)
library(lsa)
library(foreach)
library(doParallel)

# Source the files to run the simulation
source("R/generate_synthetic_mutation_counts.R")
source("R/plot_signature.R")
source("R/initialize_prior.R")

# Simulate the data
## I = n. of channels in each signature
## J = n. of patients
## K = n. of COSMIC signatures
# Load signatures
sigs <- read.COSMIC.sigs(version = "v3")  ## I x K
# Load pcawg data
counts_pcawg <- read.csv("data/WGS_PCAWG.96.ready.tsv", sep = "\t", row.names = 1)  ## I x J
names_pcawg <- str_split(colnames(counts_pcawg), pattern = "\\.\\.", simplify = TRUE)
# Select a cancer 
table(names_pcawg[, 1])
cancer <- "Biliary.AdenoCA"
id_cancer <- which(names_pcawg[, 1] == cancer)
counts <- as.matrix(counts_pcawg[, id_cancer])

# Obtain the true loadings
sig_to_include <- paste0("SBS", c("1", "2", "3", "5", "9", "12", "13", "15", 
                                  "17a", "17b", "18", "21", "22", "24", "32", "40", "44"))
loadings <- loadings.nnls(counts, sigs[, sig_to_include])
# Simulate the counts
set.seed(10)
X <- simulate.counts(loadings, sigs[, sig_to_include], "poisson")

# Run 4 models

# low variance
samples <- 100
burnin <- 500


S1 <- do.call("cbind", apply(sigs, 2, function(x) suppressWarnings(dens.max.dir(x, V = 0.01))))

S2<- do.call("cbind", apply(sigs, 2, function(x) suppressWarnings(dens.max.dir(x, V = 0.0001))))


list <- apply(sigs, 2, function(x) suppressWarnings(dens.max.dir(x, V = 0.0001)))

lapply(list, length)

sigs <- read.COSMIC.sigs(version = "v3") 
dens.max.dir(center = sigs[, "SBS25"])


out <- CompressiveNMF(X, K = K, S = S, 
                      betah = rep(1, ncol(S)), 
                      nsamples = nsamples, 
                      burnin = burnin, 
                      swap_prior = FALSE, 
                      ncores = 1,
                      nchains = 4)

# high variance

# large precision, fixed 

# small precision, fixed 












# K <- nrow(loadings)
# full.normalized.loadings <- loadings / sum(loadings)
# 
# CompressiveNMF:::plo
# 
# rowMeans(prelim.loadings)
# colSums(prelim.loadings)
# apply(prelim.loadings, 1, function(x) sum(x>0))/35
# 
# sig_to_include <- paste0("SBS", c("1", "2", "3", "5", "9", "12", "13", "15", "17a", "17b", "18", "21", "22", "24", "32", "40", "44"))
# sigs[, sig_to_include]
# 
# 
# nmf_update.KL.w(v, w, h, nbterms = 0L, ncterms = 0L, copy = TRUE)
# 
# loadings_old <- loadings.nnls(counts, sigs[, sig_to_include]) 
# loadings <- matrix(rgamma(17 * 35, 2, 2), nrow = 17)#loadings.nnls(counts, sigs[, sig_to_include])
# diff <- 1
# while(diff > 1e-4){
#   loadings_new <- NMF:::nmf_update.KL.h_R(v = counts, w = sigs[, sig_to_include], h = loadings)
#   diff <- max(abs(loadings_new - loadings))
#   loadings <- loadings_new
# }
# 
# roun
# plot(apply(round(loadings, 4), 1, function(x) sum(x>0))/35,apply(loadings_old, 1, function(x) sum(x>0))/35)
# 
# 
# apply(round(loadings_old, 4), 1, function(x) sum(x>0))/35
# apply(round(loadings, 4), 1, function(x) sum(x>0))/35
# loadings
# 
# apply(round(loadings, 4), 1, function(x) sum(x>0))/35
# rowMeans(loadings)
# 
# #------------------------------------------ Functions to load the data
# load_simulation_data <- function(path){
#   X <- as.matrix(read_tsv(path))
#   X <- X[, -1]
#   colnames(X) <- paste0("sample", 1:ncol(X))
#   rownames(X) <- rownames(cosmic_v2)
#   return(X)
# }
# 
# load_priors <- function(path){
#   S <- as.data.frame(read_csv(path, show_col_types = FALSE))
#   rownames(S) <- S[, 1]
#   return(as.matrix(S[, -1]))
# }
# 
# #------------------------------------------
# 
# # Load the cosmic signatures
# cosmic_v2 <- as.data.frame(read_csv("data/cosmic_v2.csv"))
# rownames(cosmic_v2) <- cosmic_v2[, 1]
# cosmic_v2 <- as.matrix(cosmic_v2[, -1])
# 
# # Load the priors
# priors <- list("tuned" = load_priors("data/cosmic_simulation/prior_tuned_alpha.csv") * cosmic_v2, 
#                "fixed" = load_priors("data/cosmic_simulation/prior_fixed_alpha.csv"),
#                "density_max" = load_priors("data/cosmic_simulation/prior_density_max.csv"))
# 
# # Load the data for the simulation
# data <- list("basic" = load_simulation_data("data/cosmic_simulation/synthetic-38-lung-adenoca-all-seed-1.tsv"), 
#              "perturbed" = load_simulation_data("data/cosmic_simulation/synthetic-38-lung-adenoca-all-seed-1-perturbed-0.0025.tsv"),
#              "overdispersed" = load_simulation_data("data/cosmic_simulation/synthetic-38-lung-adenoca-all-seed-1-overdispersed-2.0.tsv"))
# 
# 
# # Run all methods for all data
# cases_all <- data.frame(expand_grid(prior = names(priors), data = names(data), K = c(0, 10)))
# registerDoParallel(nrow(cases_all))
# 
# nsamples <- 2000
# burnin <- 10000
# 
# run <- FALSE
# if(run) {
#   set.seed(10, kind = "L'Ecuyer-CMRG")
#   output <- foreach(i = 1:nrow(cases_all)) %dopar% {
#     prior <- cases_all$prior[i]
#     data_type <- cases_all$data[i]
#     K <- cases_all$K[i]
#     X <- data[[data_type]]
#     S <- priors[[prior]]
#     out <- CompressiveNMF(X, K = K, S = S, betah = rep(1, ncol(S)), 
#                           nsamples = nsamples, burnin = burnin, swap_prior = FALSE, 
#                           ncores = 1, 
#                           nchains = 4)
#     out
#   }
# }
# 
# 
# saveRDS(output, file = "output/cosmic_simulation.rds.gzip", compress = "gzip")
# 
# # Read the output now!
# output <- readRDS(file = "output/cosmic_simulation.rds.gzip")
# loadings <- as.matrix(read_tsv("data/cosmic_simulation/synthetic-38-lung-adenoca-all-loadings.tsv", col_names = FALSE))
# Rmat <- as.matrix(cosmic2[, rowSums(loadings) != 0])
# Theta_true <- as.matrix(loadings[rowSums(loadings) != 0, ])
# rownames(Theta_true) <- colnames(Rmat)
# 
# res <- data.frame()
# for(i in 1:length(output)){
#   out <- output[[i]]
#   prior <- cases_all$prior[i]
#   data_type <- cases_all$data[i]
#   K <- cases_all$K[i]
#   X <- data[[data_type]]
#   data_temp <- list(X = X, Theta = Theta_true, Rmat = Rmat)
#   res_temp <- data.frame(t(Postprocess_Compressive(out, data_temp)$results))
#   res <- rbind(res, res_temp)
# }
# 
# cbind(cases_all, res)
# 
# Xout$Signatures
# Rhat <- apply(out$mcmc_out[[out$selected_chain]]$Signatures, c(2, 3), mean)
# colnames(Rhat) <- colnames(out$mcmc_out[[out$selected_chain]]$Mu)
# rownames(Rhat) <- rownames(cosmic_v2)
# 
# plot(Rhat[, 1], cosmic_v2[, 1])
# 
# Compressive
# 
# rowMeans(loadings)
# 
# 
# 
# 
# 
# 
# 
# 
# cosine(cosmic_v2)
# 
# cases_all
# CompressiveNMF:::plot_matrix_signature_v2(output[[13]]$Signatures[, 1:5])
# 
# CompressiveNMF:::plot_matrix_signature_v2(cosmic_v2[, 1:10])
# CompressiveNMF:::plot_matrix_signature_v2(cosmic_v2[, 11:20])
# CompressiveNMF:::plot_matrix_signature_v2(cosmic_v2[, 21:30])
# 
# 
# cosmic2 <- as.data.frame(read_table("COSMIC_v2_SBS_GRCh37.txt"))
# rownames(cosmic2) <- cosmic2[, 1]
# cosmic2 <- cosmic2[rownames(cosmic_v2),-1]
# 
# rownames(cosmic_v2) == rownames(cosmic2)
# 
# plot(cosmic_v2[, 30], cosmic2[, 30])
# hist(c(as.matrix(cosmic_v2 - cosmic2)))
# CompressiveNMF:::print.CompressiveNMF(output[[14]])
# CompressiveNMF:::print.CompressiveNMF(output[[13]])
# 
# CompressiveNMF:::print.CompressiveNMF(output[[1]])
# CompressiveNMF:::print.CompressiveNMF(output[[2]])
# CompressiveNMF:::print.CompressiveNMF(output[[3]])
# CompressiveNMF:::print.CompressiveNMF(output[[4]])
# 
# 
# load("../CompressiveNMF/data/Cosmic_data_no_artifacts.rdata")
# hist(cosine(as.matrix(cosmic_data[, -c(1:3)]))[lower.tri(cosine(as.matrix(cosmic_data[, -c(1:3)])))])
# 
# hist(cosine(as.matrix(cosmic_v2))[lower.tri(cosine(as.matrix(cosmic_v2)))])
# hist(cosine(as.matrix(cosmic_data[, -c(1:3)]))[lower.tri(cosine(as.matrix(cosmic_data[, -c(1:3)])))])
# 
# 
# library(pheatmap)
# 
# pheatmap(cosine(as.matrix(cosmic_data[, -c(1:3)])), cluster_rows = FALSE, cluster_cols = FALSE)
# 
# pheatmap(cosine(as.matrix(cosmic_v2)), cluster_rows = FALSE, cluster_cols = FALSE)
# 
# 
# 
# cosmic2[, ]
# 
# 
# 
# 
# 
# 
# 
# 
# # # Test the methods
# # X <- data$basic
# # out_tuned <- CompressiveNMF(X, K = 0, S = cosmic_v2, betah = priors$tuned[1, ], 
# #                             nsamples = 1000, burnin = 3000, swap_prior = FALSE)
# # 
# # out_density_max <- CompressiveNMF(X, K = 10, S = priors$density_max, betah = rep(1, ncol(priors$density_max)),
# #                                   nsamples = 500, burnin = 2000, swap_prior = FALSE)
# # 
# # out_density_max$RelWeights
# # plot(priors$fixed[, 1], cosmic_v2[, 1])
# # plot()
# # 
# # priors$fixed[, 1]/cosmic_v2[, 1]
# # 
# # out_tuned$Signatures
# # CompressiveNMF:::print.CompressiveNMF(out_tuned)
# # CompressiveNMF:::print.CompressiveNMF(out_density_max)
# # out_tuned$RelWeights
# # 
# # plot(280.0912 * cosmic_v2[,1])
# # mean(apply(rdirichlet(1000, alpha = 280.0912 *cosmic_v2[,1]), 1, function(x) cosine(cosmic_v2[, 1], x)))
# # cosine()
# # 
# # loadings <- as.matrix(read_tsv("data/cosmic_simulation/synthetic-38-lung-adenoca-all-loadings.tsv", col_names = FALSE))
# # loadings
# # 
# # 
# # sum(rowSums(loadings) > 0)
# # ncol(out_tuned$Signatures)
# # out_tuned$RelWeights
# # 
# # plot(out_tuned$mcmc_out[[1]]$Mu[, 5])
# # 
# # plot(X, out_tuned$Signatures %*% out_tuned$Weights)
# # plot(X, out_density_max$Signatures %*% out_density_max$Weights)
# # 
# # out_density_max$Weights[1, ] - 
# # CompressiveNMF:::plot_matrix_signature_v2(out_density_max$Signatures)
# # CompressiveNMF:::plot_matrix_signature_v2(out_tuned$Signatures, add_cosine_cosmic = FALSE)
# # 
# # as.matrix(rowSums(X))
# # CompressiveNMF:::plot_matrix_signature_v2(as.matrix(rowSums(X)), add_cosine_cosmic = FALSE)
# # 
# # w <- apply(out_density_max$mcmc_out[[1]]$Weights, c(2,3), mean)
# # dim(w)  
# # plot(w, loadings)  
# # 
# # 
# # 
