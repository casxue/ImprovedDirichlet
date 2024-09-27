# source("initialize_prior.R")

# library(Ternary)
library(parallel)
library(foreach)
# library(rootSolve)
library(rstan)
library(latex2exp)


## ACF and effective sample size, but not broken
ac.new <- function (vec, trunc = length(vec) - 1) {
    ac <- acf(vec, lag.max = trunc, plot = FALSE)[1:trunc][[1]]
    if (sum(is.nan(ac)) > 0)  {
        ac <- rep(1, length(ac))
    }
    return(ac)
}
ESS <- function (vec, M = length(vec), trunc = M - 1) {
    ac <- ac.new(vec, trunc)
    delta <- 1:trunc
    Teff <- M / (1 + 2 * sum((1 - delta / M) * ac))
    return(Teff)
}


######## Proposal distributions

proposal.params <- list()

## 1A: mean, fixed concentration
proposal.params[[1]] <- function (u, alpha = 5) {
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
proposal.params[[4]] <- function (med, alpha = 5) {
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

## 1C: density maximization, fixed concentration
proposal.params[[7]] <- function (center, alpha = 5) {
    return(dens.max.conc(center, alpha))
} 

## 2C: density maximization, fixed variance
proposal.params[[8]] <- function (center, var = 0.1) {
    return(dens.max.var(center, var))
} 

## 3C: maximum density, adaptive variance
proposal.params[[9]] <- function (center) {
    sd <- min(center, 1 - center, sqrt(0.1))
    return(dens.max.var(center, sd^2))
} 


## sampler function
sample.MH <- function (exp, initial = 0.25, alpha = 5, V = 0.1, n.burn = 0, n.samps = 1e4, seed = 2024, a0 = 1, b0 = 1000, is.mixture = FALSE, verbose = FALSE)  {
    set.seed(seed)
    
    samples <- numeric(n.samps)
    current <- initial
    counter <- 0
    
    for (j in 1:(n.burn + n.samps))  {
        ## compute proposal
        if (exp %% 3 == 1)  {
            params.curr <- proposal.params[[exp]](current, alpha)
        } else if (exp %% 3 == 2)  {
            params.curr <- proposal.params[[exp]](current, V)
        } else  {
            params.curr <- proposal.params[[exp]](current)
        }

        if (verbose && j > n.burn && params.curr$a + params.curr$b == 1e-8)  {
            counter <- counter + 1
        } 
        
        ## compute parameters of return probability
        prop <- rbeta(1, params.curr$a, params.curr$b)
        if (exp %% 3 == 1)  {
            params.prop <- proposal.params[[exp]](prop, alpha)
        } else if (exp %% 3 == 2)  {
            params.prop <- proposal.params[[exp]](prop, V)
        } else  {
            params.prop <- proposal.params[[exp]](prop)
        }
        
        ## compute acceptance ratio
        if (is.mixture)  {
            log.accept <- log(0.75 * dbeta(prop, 2, 5) + 0.25 * dbeta(prop, 10, 2)) - 
                log(0.75 * dbeta(current, 2, 5) + 0.25 * dbeta(current, 10, 2)) + 
                dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
        } else  {
            log.accept <- dbeta(prop, a0, b0, log = TRUE) - dbeta(current, a0, b0, log = TRUE) +
                dbeta(current, params.prop$a, params.prop$b, log = TRUE) - dbeta(prop, params.curr$a, params.curr$b, log = TRUE)
        }

        ## accept/reject
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
## NOT UPDATED YET
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

