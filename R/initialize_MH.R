library(dplyr)
library(gtools)
library(LaplacesDemon)
library(scales)
library(Ternary)
library(parallel)
library(foreach)
library(rootSolve)
library(rstan)
library(latex2exp)


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
    ac <- acf(vec, lag.max = trunc, plot = FALSE)[1:trunc][[1]]
    if (sum(is.nan(ac)) > 0)  {
        ac <- rep(1, length(ac))
    }
    delta <- 1:trunc
    Teff <- M / (1 + 2 * sum((1 - delta / M) * ac))
    return(Teff)
}


## Find quantile of v for Beta distribution with mean u/N and concentration alpha
beta.quantile <- function (u, v, alpha, N) {
    u <- u / N
    return(pbeta(v, alpha * u, alpha * (1 - u)))
}


## Find mean u for Beta distribution with median c and concentration alpha
arg.median.beta <- function (c, alpha, precision = 5) {
    M <- 10^precision
    u <- mean(binsearch(beta.quantile, c(0, M), v = c, alpha = alpha, N = M, target = 0.5)$where) / M
    return(u)
}


## Find variance for Beta distribution with given median and log concentration parameter floor.log.conc + n / N
variance.from.median.concentration <- function (n, median, floor.log.conc, N) {
    alpha <- 2^(floor.log.conc + n / N)
    u <- arg.median.beta(median, alpha)
    var <- u * (1 - u) / (alpha + 1)
    return(var)
}

## Identify mean & concentration for Beta distribution based on median & variance via coordinate descent
proposal.med.var <- function (c, var, precision = 3) {
    var.target <- var
    conc <- 2
    
    u <- arg.median.beta(c, conc)
    var <- u * (1 - u) / (conc + 1)
    
    if (var > var.target) {
        while (var > var.target) {
            conc <- conc * 2
            u <- arg.median.beta(c, conc)
            var <- u * (1 - u) / (conc + 1)
        }
        floor.log.conc <- log2(conc) - 1
    } else if (var < var.target) {
        while (var < var.target) {
            conc <- conc / 2
            u <- arg.median.beta(c, conc)
            var <- u * (1 - u) / (conc + 1)
        }
        floor.log.conc <- log2(conc)
    }
    
    M <- 10^precision
    log.conc <- mean(binsearch(variance.from.median.concentration, c(0, M), median = c, floor.log.conc = floor.log.conc, N = M, target = var.target)$where) / M + floor.log.conc
    conc <- 2^log.conc
    u <- arg.median.beta(c, conc)
    var <- u * (1 - u) / (conc + 1)
    
    return(list(mean = u, alpha = conc, var.actual = var))
}


## Minimum and maximum valid mean for a given variance, with small buffer
u.min <- function (S, eps = 1e-8)  {
    return((1 - sqrt(1 - 4 * S)) / 2 + eps)
}

u.max <- function (S, eps = 1e-8)  {
    return((1 + sqrt(1 - 4 * S)) / 2 - eps)
}


## Evaluate density at target point loc under a Beta with mean u and specified variance or concentration 
prop.dens <- function (u, loc = 0.001, var = 0, conc = 0, log = TRUE) {
    if (var > 0 && conc == 0) {
        alpha <- u * (1 - u) / var - 1
    } else if (var == 0 && conc > 0) {
        alpha <- conc
    } else {
        print("Need valid variance or concentration.")
        return(NA)
    }
    return(dbeta(loc, alpha * u, alpha * (1 - u), log = log))
}


## Objective function for density maximization
h.prime <- function (u, center = 0.001, var = 0.01)  {
    return((log(center) - psigamma(u^2 * (1 - u) / var - u)) * (u * (2 - 3 * u) / var - 1) +
           (log(1 - center) - psigamma(u * (1 - u)^2 / var - (1 - u))) * ((3 * u^2 - 4 * u + 1) / var + 1) +
           psigamma(u * (1 - u) / var - 1) * (1 - 2 * u) / var)
}


######## Proposal distributions

proposal.params <- list()

## 1A: mean, fixed concentration
proposal.params[[1]] <- function (u, alpha = 2) {
    return(list(a = alpha * u, b = alpha * (1 - u)))
}

## 2A: mean, fixed variance
proposal.params[[2]] <- function (u, var = 0.1) {
    alpha <- max(param.variance(u, var), 1e-8)
    return(list(a = alpha * u, b = alpha * (1 - u)))
}

## 3A: mean, adaptive variance
proposal.params[[3]] <- function (u) {
    sd <- min(u, 1 - u, sqrt(0.1))
    alpha <- max(param.variance(u, sd^2), 1e-8)
    return(list(a = alpha * u, b = alpha * (1 - u)))
}

## 1B: median, fixed concentration
proposal.params[[4]] <- function (med, alpha = 2) {
    u <- arg.median.beta(med, alpha)
    return(list(a = alpha * u, b = alpha * (1 - u)))
}

## 2B: median, fixed variance
proposal.params[[5]] <- function (med, var = 0.1) {
    p <- proposal.med.var(med, var)
    return(list(a = p$alpha * p$mean, b = p$alpha * (1 - p$mean), iters = p$iters))
}

## 3B: median, adaptive variance
proposal.params[[6]] <- function (med) {
    sd <- min(med, 1 - med, sqrt(0.1))
    p <- proposal.med.var(med, sd^2)
    return(list(a = p$alpha * p$mean, b = p$alpha * (1 - p$mean), iters = p$iters))
} 

## C: density maximization, fixed variance
proposal.params[[7]] <- function (center, var = 0.1) {
    roots <- uniroot.all(h.prime, lower = u.min(var), upper = u.max(var), center = center, var = var)
    ind <- sapply(roots, prop.dens, loc = center, var = var) %>% which.max()
    if (length(ind) > 1)  ind <- ind[1]
    
    u <- roots[ind]
    alpha <- u * (1 - u) / var - 1
    return(list(a = alpha * u, b = alpha * (1 - u)))
} 


## sampler function
sample.MH <- function (exp, initial = 0.25, S = 0.05, n.burn = 0, n.samps = 1e4, seed = 2024, a0 = 1, b0 = 1000, verbose = FALSE)  {
    set.seed(seed)
    
    samples <- numeric(n.samps)
    current <- initial
    counter <- 0
    
    for (j in 1:(n.burn + n.samps))  {
        if (exp %in% c(2, 5, 7))  {
            params.curr <- proposal.params[[exp]](current, S)
        }
        else  {
            params.curr <- proposal.params[[exp]](current)
        }

        if (verbose && j > n.burn && params.curr$a + params.curr$b == 1e-8)  {
            counter <- counter + 1
        } 
        
        prop <- rbeta(1, params.curr$a, params.curr$b)
        if (exp %in% c(2, 5, 7))  {
            params.prop <- proposal.params[[exp]](prop, S)
        }
        else  {
            params.prop <- proposal.params[[exp]](prop)
        }
        
        log.accept <- dbeta(prop, a0, b0, log = TRUE) - dbeta(current, a0, b0, log = TRUE) +
            dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
        # log.accept <- log(0.75 * dbeta(prop, 2, 5) + 0.25 * dbeta(prop, 10, 2)) - 
        #     log(0.75 * dbeta(current, 2, 5) + 0.25 * dbeta(current, 10, 2)) + 
        #     dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
        if (!is.nan(log.accept) && exp(log.accept) > runif(1))  {
            current <- prop
        } 
        if (j > n.burn) {
            samples[j - n.burn] <- current
        }
    }
    
    if (verbose)  {
        if (exp == 2 || exp == 5)  {
            sprintf("Concentration parameter underflow percentage in run %d: %.2f", i, counter / n.samps * 100) %>% print()
        }
    }
    
    return(samples)
}

## sampler function with 0/1 bit to handle numerical stability near 1
sample.MH.augmented <- function (exp, initial = 0.25, S = 0.05, n.burn = 0, n.samps = 1e4, seed = 2024, a0 = 1, b0 = 1000, verbose = FALSE)  {
    set.seed(seed)
    
    samples <- matrix(0L, nrow = n.samps, ncol = 2)
    current <- initial
    curr.bit <- (initial <= 0.5)
    counter <- 0
    
    for (j in 1:(n.burn + n.samps))  {
        ## Get parameters for Beta centered at current value
        if (exp %in% c(2, 5, 7))  {
            params.curr <- proposal.params[[exp]](current, S)
        }
        else  {
            params.curr <- proposal.params[[exp]](current)
        }
        
        if (verbose && j > n.burn && params.curr$a + params.curr$b == 1e-8)  {
            counter <- counter + 1
        } 
        
        ## Get proposal via probability integral transform; flip if necessary
        u <- runif(1)
        prop <- qbeta(u, params.curr$a, params.curr$b)
        if (prop <= 0.5)  {
            prop.bit <- curr.bit
        } else  {
            prop <- qbeta(1 - u, params.curr$b, params.curr$a)
            prop.bit <- 1 - curr.bit
        }
        
        ## Get parameters for Beta centered at proposal
        if (exp %in% c(2, 5, 7))  {
            params.prop <- proposal.params[[exp]](prop, S)
        }
        else  {
            params.prop <- proposal.params[[exp]](prop)
        }
        
        ## Compute acceptance probability, taking into account curr.bit and prop.bit
        if (curr.bit == 0 & prop.bit == 0)  {
            log.accept <- dbeta(prop, a0, b0, log = TRUE) - dbeta(current, a0, b0, log = TRUE) +
                dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
            # log.accept <- log(0.75 * dbeta(prop, 2, 5) + 0.25 * dbeta(prop, 10, 2)) - 
            #     log(0.75 * dbeta(current, 2, 5) + 0.25 * dbeta(current, 10, 2)) + 
            #     dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
        } else if (curr.bit == 0 & prop.bit == 1)  {
            log.accept <- dbeta(prop, b0, a0, log = TRUE) - dbeta(current, a0, b0, log = TRUE) +
                dbeta(current, params.prop$b, params.prop$a, log = TRUE) - dbeta(prop, params.curr$b, params.curr$a, log = TRUE)
            # log.accept <- log(0.75 * dbeta(prop, 5, 2) + 0.25 * dbeta(prop, 2, 10)) - 
            #     log(0.75 * dbeta(current, 2, 5) + 0.25 * dbeta(current, 10, 2)) + 
            #     dbeta(current, params.prop$b, params.prop$a, log = TRUE) - dbeta(prop, params.curr$b, params.curr$a, log = TRUE)
        } else if (curr.bit == 1 & prop.bit == 0)  {
            log.accept <- dbeta(prop, a0, b0, log = TRUE) - dbeta(current, b0, a0, log = TRUE) +
                dbeta(current, params.prop$b, params.prop$a, log = TRUE) - dbeta(prop, params.curr$b, params.curr$a, log = TRUE)
            # log.accept <- log(0.75 * dbeta(prop, 2, 5) + 0.25 * dbeta(prop, 10, 2)) - 
            #     log(0.75 * dbeta(current, 5, 2) + 0.25 * dbeta(current, 2, 10)) + 
            #     dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
        } else  {
            log.accept <- dbeta(prop, b0, a0, log = TRUE) - dbeta(current, b0, a0, log = TRUE) +
                dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
            # log.accept <- log(0.75 * dbeta(prop, 5, 2) + 0.25 * dbeta(prop, 2, 10)) - 
            #     log(0.75 * dbeta(current, 5, 2) + 0.25 * dbeta(current, 2, 10)) + 
            #     dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
        }
        
        ## accept/reject
        if (!is.nan(log.accept) && exp(log.accept) > runif(1))  {
            current <- prop
            curr.bit <- prop.bit
        } 
        if (j > n.burn) {
            samples[j - n.burn,] <- c(current, curr.bit)
        }
    }
    
    if (verbose)  {
        if (exp == 2 || exp == 5)  {
            sprintf("Concentration parameter underflow percentage in run %d: %.2f", i, counter / n.samps * 100) %>% print()
        }
    }
    
    return(samples)
}

