setwd("/n/miller_lab/csxue/dirichlet/ImprovedDirichlet")

source("R/initialize_prior.R")
source("R/initialize_MH.R")

## Target distribution
a0 <- 0.5
b0 <- 0.5

## Proposal
alpha <- 5
V <- 0.1

## Initialization
N <- 100 # number of simulations
B <- 100 # number of burn-in samples
M <- 1e4 # length of each sampler
n.exp <- length(proposal.params)

## Sampling
print("Starting simulation") 
beta.samps <- foreach (i = 1:n.exp) %do% {
    foreach (n = 1:N) %do% {
        sample.MH(i, seed = n, n.burn = B, n.samps = M, alpha = alpha, V = V, a0 = a0, b0 = b0)
    } 
}
saveRDS(beta.samps, file = paste0("beta-samples_", a0, "-", b0, "_alpha-", alpha, "_var-", V, ".rds"))

## plotting
names <- paste0(rep(1:3, 3), rep(c("A", "B", "C"), each = 3))
beta.samps <- readRDS("beta-samples_1-1000_alpha-5_var-0.1.rds")
a0 <- 1
b0 <- 1000

## autocorrelation
lag.max <- 50
lags <- 0:lag.max
plot(0, type = "n", xlim = c(0, lag.max), ylim = c(0, 1), 
     xlab = "Lag", ylab = "Autocorrelation",
     cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.25)
for (j in c(1:3, 8))  {
    lapply(beta.samps[[j]], ac.new, trunc = 50) %>% 
        do.call(cbind, .) %>% 
        rowMeans() %>%
        c(1, .) %>% ## start with autocorrelation = 1 for lag of 0
        points(lags, ., type = "l", col = j, lwd = 3)
}
legend(7, 1, bty = "n", y.intersp = 1.25, 
       legend = c("I. max density at x; V = 0.1", 
                  TeX("II. mean = x; $\\alpha$ = 2"), 
                  "III. mean = x; V = 0.1", 
                  "IV. mean = x; adaptive variance"),
       col = c(8, 1:3), lwd = 3, seg.len = 1, cex = 1.5)

## K-S test for a Beta target
par(mfcol = c(3, 3))
for (k in 1:n.exp)  {
    suppressWarnings(
        Ds <- lapply(beta.samps[[k]], function (lst) {ks.test(lst, pbeta, a0, b0)$statistic}) %>% 
                unlist()
    )
    hist(Ds, main = names[k], breaks = 20, xlim = c(0, 1))
}
    
