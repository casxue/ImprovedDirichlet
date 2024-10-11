source("R/initialize_prior.R")
source("R/initialize_MH.R")
library(reshape)
library(HDInterval)


coverage <- function (theta0, method, n, cred = 0.95, alpha = 1, prior.loc = theta0, use.hpd = TRUE)  {
    if (method == "mean")  {
        a0 <- alpha * prior.loc
        b0 <- alpha * (1 - prior.loc)
    } else if (method == "median")  {
        params <- proposal.params[[4]](prior.loc, alpha)
        a0 <- params$a
        b0 <- params$b
    } else if (method == "maxdens")  {
        params <- proposal.params[[7]](prior.loc, alpha)
        a0 <- params$a
        b0 <- params$b
    } else  {
        stop("Invalid method")
    }
    
    if (use.hpd)  {
        ## Highest posterior density interval
        ints <- lapply(0:n, function (k) hdi(qbeta, cred, shape1 = a0 + k, shape2 = b0 + n - k))
        is.covered <- sapply(0:n, function (k) (theta0 > ints[[k + 1]][['lower']] & theta0 < ints[[k + 1]][['upper']]))
    } else  {
        ## Symmetric credible interval
        quantiles <- sapply(0:n, function (k) pbeta(theta0, a0 + k, b0 + n - k))
        is.covered <- quantiles > (1 - cred) / 2 & quantiles < (1 - (1 - cred) / 2)
    }
    
    weights <- sapply(0:n, function (k) dbinom(k, size = n, prob = theta0))
    return(sum(weights * is.covered))
}



min.log.c <- -6
max.log.c <- 0
cs <- 10^(seq(from = min.log.c, to = max.log.c, by = 0.05) %>% head(-1))

alpha <- 10
n <- 100

## Prior centered at true parameter
mean.coverage <- sapply(cs, function (theta) coverage(theta, "mean", n, alpha = alpha))
median.coverage <- sapply(cs, function (theta) coverage(theta, "median", n, alpha = alpha))
maxdens.coverage <- sapply(cs, function (theta) coverage(theta, "maxdens", n, alpha = alpha))

plot(cs, mean.coverage, type = "l", col = 2, lwd = 3, log = "x", ylim = c(0, 1), xaxt = "n",
     xlab = TeX("True parameter $\\theta_0$"), ylab = "Coverage", main = TeX("Prior centered at true parameter"),
     cex.lab = 1.5, cex.main = 1.5)
points(cs, median.coverage, type = "l", col = 3, lwd = 3)
points(cs, maxdens.coverage, type = "l", col = 4, lwd = 3)

axis(1, at = 10^(-6:0), labels = TeX(paste0("$10^{", -6:0, "}$")))
legend(1e-6, 0.8, bty = "n", x.intersp = 0.75, y.intersp = 1.5, lwd = 3, seg.len = 1.25, cex = 1.25, 
       legend = c("mean method", "median method", "maximum density"), col = 2:4)


## Prior centered at 1e-3
mean.coverage <- sapply(cs, function (theta) coverage(theta, "mean", n, alpha = alpha, prior.loc = 1e-3))
median.coverage <- sapply(cs, function (theta) coverage(theta, "median", n, alpha = alpha, prior.loc = 1e-3))
maxdens.coverage <- sapply(cs, function (theta) coverage(theta, "maxdens", n, alpha = alpha, prior.loc = 1e-3))

