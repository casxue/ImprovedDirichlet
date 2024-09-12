library(dplyr)
library(gtools)
library(LaplacesDemon)
library(scales)
library(parallel)
library(foreach)
library(rootSolve)
library(rstan)
library(latex2exp)


## concentration parameter corresponding to specified mean and variance for Beta, or for the smallest channel for Dirichlet
param.mean.variance <- function (u, delta) {
    if (length(u) == 1)  {
        return(u * (1 - u) / delta - 1)
    } else  {
        return(min(u) * (1 - min(u)) / delta - 1)
    }
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
    V <- u * (1 - u) / (alpha + 1)
    return(V)
}

## Identify mean & concentration for Beta distribution based on median & variance via coordinate descent
proposal.med.var <- function (c, V, precision = 3) {
    var.target <- V
    conc <- 2
    
    u <- arg.median.beta(c, conc)
    V <- u * (1 - u) / (conc + 1)
    
    if (V > var.target) {
        while (V > var.target) {
            conc <- conc * 2
            u <- arg.median.beta(c, conc)
            V <- u * (1 - u) / (conc + 1)
        }
        floor.log.conc <- log2(conc) - 1
    } else if (V < var.target) {
        while (V < var.target) {
            conc <- conc / 2
            u <- arg.median.beta(c, conc)
            V <- u * (1 - u) / (conc + 1)
        }
        floor.log.conc <- log2(conc)
    }
    
    M <- 10^precision
    log.conc <- mean(binsearch(variance.from.median.concentration, c(0, M), median = c, floor.log.conc = floor.log.conc, N = M, target = var.target)$where) / M + floor.log.conc
    conc <- 2^log.conc
    u <- arg.median.beta(c, conc)
    V <- u * (1 - u) / (conc + 1)
    
    return(list(mean = u, alpha = conc, var.actual = V))
}


## Minimum and maximum valid mean for a given variance, with small buffer
u.min <- function (V, eps = 1e-8)  {
    return((1 - sqrt(1 - 4 * V)) / 2 + eps)
}

u.max <- function (V, eps = 1e-8)  {
    return((1 + sqrt(1 - 4 * V)) / 2 - eps)
}


## Evaluate density at target point loc under a Beta with mean u and specified variance or concentration 
prop.dens <- function (u, loc = 0.001, V = 0, conc = 0, log = TRUE) {
    if (V > 0 && conc == 0) {
        alpha <- u * (1 - u) / V - 1
    } else if (V == 0 && conc > 0) {
        alpha <- conc
    } else {
        print("Need valid variance or concentration.")
        return(NA)
    }
    return(dbeta(loc, alpha * u, alpha * (1 - u), log = log))
}


## Objective function for density maximization with variance
h.prime <- function (u, center = 0.001, V = 0.01)  {
    return((log(center) - psigamma(u^2 * (1 - u) / V - u)) * (u * (2 - 3 * u) / V - 1) +
           (log(1 - center) - psigamma(u * (1 - u)^2 / V - (1 - u))) * ((3 * u^2 - 4 * u + 1) / V + 1) +
           psigamma(u * (1 - u) / V - 1) * (1 - 2 * u) / V)
}


## Objective function for density maximization with concentration parameter
g.prime <- function (u, center = 0.001, conc = 2)  {
    return(conc * log(center) - conc * log(1 - center) - psigamma(conc * u) * conc + psigamma(conc * (1 - u)) * conc)
}


## density maximization parameters, fixed variance
dens.max.var <- function (center, V = 0.01) {
    roots <- uniroot.all(h.prime, lower = u.min(V), upper = u.max(V), center = center, V = V)
    ind <- sapply(roots, prop.dens, loc = center, V = V) %>% which.max()
    if (length(ind) > 1)  ind <- ind[1]
    
    u <- roots[ind]
    alpha <- u * (1 - u) / V - 1
    return(list(a = alpha * u, b = alpha * (1 - u)))
} 


## density maximization parameters, fixed concentration
dens.max.conc <- function (center, conc = 1000, eps = 1e-8) {
    roots <- uniroot.all(g.prime, lower = eps, upper = 1 - eps, center = center, conc = conc)
    ind <- sapply(roots, prop.dens, loc = center, conc = conc) %>% which.max()
    if (length(ind) > 1)  ind <- ind[1]
    
    u <- roots[ind]
    return(list(a = conc * u, b = conc * (1 - u)))
} 


## density maximization parameters for Dirichlet
dens.max.dir <- function (center, V = 0.01) {
    i.star <- which.min(sapply(center, function (c) {min(c, 1 - c)}))[1]
    lst <- dens.max.var(center[i.star], V)
    conc <- lst$a + lst$b
    
    u <- numeric()
    for (i in 1:length(center)) {
        if (i == i.star)  {
            u <- c(u, lst$a / conc)
        } else  {
            u <- c(u, dens.max.conc(center[i], conc)$a / conc)
        }
    }
    return(u * conc)
} 



## Load COSMIC v2 signatures
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
cosmic <- read.csv("COSMIC_v2_SBS_GRCh37.txt", sep = "\t", header = T, row.names = 1)[channel.order,]


