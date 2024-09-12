library(dplyr)
library(latex2exp)
library(gtools)
library(LaplacesDemon)
library(scales)
library(Ternary)
library(parallel)


##### Initialize empirical stuff

## Simulate some signatures 
I <- 3
K <- 50

set.seed(1)
sigs <- rdirichlet(K - 2, rep(1, I)) %>% t()
sigs <- cbind(sigs, c(1 / 3, 1 / 3, 1 / 3), c(1, 0, 0))

sim.overdispersed <- function (concentration, sig, n = 1e4, sigs.matrix = cosmic.v3) {
    return(rdirichlet(n, concentration * sigs.matrix[,sig]))
}

## Expectation of second-order Taylor approximation to cosine error
taylor.exp <- function (sig, param) {
    if (sum(sig) != 1)  sig <- sig / sum(sig)
    return(1 / (2 * (param + 1) * norm(sig, type = "2")^2) - sum(sig^3) / (2 * (param + 1) * norm(sig, type = "2")^4))
}

## Concentration parameter based on expectation of second-order Taylor approximation to cosine error
param.taylor <- function (sig, delta) {
    if (sum(sig) != 1)  sig <- sig / sum(sig)
    return(1 / (2 * delta * sum(sig^2)) * (1 - sum(sig^3) / sum(sig^2)) - 1)
}

## Maximum attainable cosine error for a given reference signature
max.err <- function (sig) {
    if (sum(sig) != 1)  sig <- sig / sum(sig)
    return(1 - min(sig) / sqrt(sum(sig^2)))
}

## Simulated and compute cosine errors for a single signature and fixed concentration
emp.cos.error <- function (sig, concentration, M = 1e5) {
    samps <- rdirichlet(M, sig * concentration)
    apply(samps, 1, function (u) {1 - sum(u * sig) / sqrt(sum(u^2) * sum(sig^2))}) %>% mean(na.rm = TRUE) %>% return()
}

## Simulate and compute cosine errors for a matrix of signatures and fixed concentration
plot.error <- function (concentration, sig, n = 1e4, sigs.matrix = sigs, show.plots = TRUE) {
    samples <- sim.overdispersed(concentration, sig, n, sigs.matrix)
    cos.error <- sweep(samples, MARGIN = 2, sigs.matrix[,sig], `*`) %>% rowSums()
    cos.error <- 1 - cos.error / norm(sigs.matrix[,sig], type = "2") / apply(samples, 1, norm, type = "2")
    if (show.plots) {
        hist(cos.error, xlim = c(0, 0.1), main = paste0("Cosine Error, ", 
                                                        colnames(sigs.matrix)[sig], 
                                                        ", Concentration paramemter = ", 
                                                        concentration, ", n = ", n))
    }
    return(list(samples, cos.error))
}

sig.flatness <- 1 / apply(sigs, 2, norm, type = "2")^2

## Simulate
set.seed(1)
params <- c(10, 25, seq(50, 500, by = 50))
slopes.sim <- slopes.taylor <- c()
for (param in params) {
    cos.error <- c()
    for (s in 1:K) {
        cos.error <- plot.error(param, s, show.plots = FALSE)[[2]] %>% mean() %>% c(cos.error, .)
    } 
    slopes.sim <- lm(formula = cos.error ~ sig.flatness)$coefficients[2] %>% c(slopes.sim, .)

    taylor <- apply(sigs, 2, taylor.exp, param)
    slopes.taylor <- lm(formula = taylor ~ sig.flatness)$coefficients[2] %>% c(slopes.taylor, .)
}

## Summarize
lm.sim <- lm(formula = log(slopes.sim) ~ log(params))
lm.taylor <- lm(formula = log(slopes.taylor[6:12]) ~ log(params[6:12]))

conc.param <- function (sig, delta, approx = FALSE) {
    param = 0
    if (approx) {
        param = exp(lm.taylor$coefficients[1] - log (delta * norm(sig, type = "2")^2)) - exp(-1 * lm.taylor$coefficients[2])
    } else {
        param = exp(lm.sim$coefficients[1] - log (delta * norm(sig, type = "2")^2)) - exp(-1 * lm.sim$coefficients[2])
    }
    return(unname(param))
}

## Variance-based concentration parameter 
param.variance <- function (u, delta, total.var = FALSE) {
    if (length(u) == 1 || total.var)  {
        return(sum(u * (1 - u)) / delta - 1)
    } else  {
        return(min(u) * (1 - min(u)) / delta - 1)
    }
}


## Effective sample size
ESS <- function (vec, M = length(vec), trunc = M - 1) {
    ac <- acf(vec, lag.max = trunc, plot = FALSE)[1:trunc]
    delta <- 1:trunc
    Teff <- M / (1 + 2 * sum((1 - delta / M) * ac[[1]]))
    return(Teff)
}


## Load COSMIC v2 signatures
cosmic <- read.csv("COSMIC_v2_SBS_GRCh37.txt", sep = "\t", header = T, row.names = 1)

channel.order <- c()
pyrs <- c("C", "T")
nucs <- c("A", "C", "G", "T")
for (i in pyrs)  {
    for (j in setdiff(nucs, i))  {
        for (k in nucs)  {
            for (l in nucs)  {
                channel.order <- c(channel.order, sprintf("%s[%s>%s]%s", k, i, j, l))
            }
        }
    }
}

cosmic <- cosmic[channel.order,]


