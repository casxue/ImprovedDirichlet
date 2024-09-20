require(nnls)
require(dplyr)
require(LaplacesDemon)

##### Functions #####
reorder.cosmic <- function ()  {
    pyrs <- c("C", "T")
    nucs <- c("A", "C", "G", "T")
    
    our.order <- c()
    for (i in pyrs)  {
        for (j in setdiff(nucs, i))  {
            for (k in nucs)  {
                for (l in nucs)  {
                    our.order <- c(our.order, sprintf("%s[%s>%s]%s", k, i, j, l))
                }
            }
        }
    }
    
    cosmic.order <- c()
    for (i in nucs)  {
        for (j in pyrs)  {
            for (k in setdiff(nucs, j))  {
                for (l in nucs)  {
                    cosmic.order <- c(cosmic.order, sprintf("%s[%s>%s]%s", i, j, k, l))
                }
            }
        }
    }
    
    return(match(our.order, cosmic.order))
}

read.COSMIC.sigs <- function (version = "v2", build = "GRCh37")  {
    sigs.file <- file.path("data/", paste0("COSMIC_", version, "_SBS_", build, ".txt"))
    cosmic <- read.csv(sigs.file, sep = "\t", header = T, row.names = 1)[reorder.cosmic(),]
    return(as.matrix(cosmic))
}

loadings.nnls <- function (counts, sigs = read.COSMIC.sigs())  {
  theta <- apply(counts, 2, function (y) {nnls(sigs, y)$x})
  rownames(theta) <- colnames(sigs)
  return(theta) ## K x J
}

## prune loadings before generating synthetic data
## drop any signatures that contribute less than overall.threshold percent of the total mutation burden and 
##      less than individual.threshold percent of any one patient's mutation burden
trim.loadings <- function (loadings, overall.threshold = 2, individual.threshold = 10)  {
    K <- nrow(loadings)
    full.normalized.loadings <- loadings / sum(loadings)
    individual.normalized.loadings <- loadings / colSums(loadings)[col(loadings)]
    keep <- (rowSums(full.normalized.loadings) * 100 > overall.threshold) | 
            (apply(individual.normalized.loadings, 1, max) * 100 > individual.threshold)
    loadings[!keep,] <- 0
    return(loadings)
}

## simulate new counts from specified loadings and signatures matrices
## poisson - standard model
## contamination - add param % of a subject-specific contamination signature
## overdispersed - generate from negative binomial with mu = param 
## perturbed - generate new versions of each signature for each individual, calibrated so param is the expected cosine error wrt the COSMIC signature
## Defaults: contamination 2, overdispersed 2, perturbed 0.0025
simulate.counts <- function (loadings, sigs, mode = "poisson", param = 0, seed = 1)  {
    I <- nrow(sigs)
    J <- ncol(loadings)
    K <- ncol(sigs)

    set.seed(seed)
    if (mode == "poisson")  {
        means <- sigs %*% loadings
        new.counts <- means %>%
            sapply(function (y) {rpois(1, lambda = y)}) %>%
            matrix(nrow = I, ncol = J)
    } else if (mode == "overdispersed")  {
        if (param <= 1)  {
            print("Parameter must be >1 for overdispersed DGP")
            return()
        }
        means <- sigs %*% loadings
        new.counts <- means %>%
            sapply(function (y) {rnbinom(1, size = y, mu = param)}) %>%
            matrix(nrow = I, ncol = J)
    } else if (mode == "perturbed")  {
        if (param >= 1 || param <= 0)  {
            print("Parameter must be between 0 and 1 for perturbed DGP")
            return()
        }
        ## log (mean error / sparsity) = 3.6641 - 0.9820 log (concentration)
        ## sparsity = 1 / (I * ||sig||^2)
        sparsity <- apply(sigs, 2, function (sig) {1 / (I * sum(sig^2))})
        concentration <- exp( 3.6641 - log(param / sparsity) / 0.9820 )
        
        rand.sigs <- array(0L, dim = c(I, J, K))
        for (k in 1:K)  {
            rand.sigs[,,k] <- t(rdirichlet(J, concentration[k] * sigs[,k]))
        }
        
        new.counts <- matrix(0L, nrow = I, ncol = J)
        for (j in 1:J)  {
            means <- rand.sigs[,j,] %*% trimmed.loadings[,j]
            new.counts[,j] <- sapply(means, function (y) {rpois(1, lambda = y)})
        }
    } else if (mode == "contamination")  {
        if (param <= 0)  {
            print("Parameter must be >0 for contamination DGP")
            return()
        }
        error.sigs <- t(rdirichlet(J, rep(0.5, I)))
        error.loadings <- param * 0.01 * colSums(loadings)
        new.counts <- matrix(0L, nrow = I, ncol = J)
        for (j in 1:J)  {
            means <- cbind(sigs, error.sigs[,j]) %*% c(loadings[,j] * (1 - param * 0.01), error.loadings[j])
            new.counts[,j] <- sapply(means, function (y) {rpois(1, lambda = y)})
        }
    } else  {
        print("invalid mode")
        return()
    }
    colnames(new.counts) <- colnames(loadings)
    rownames(new.counts) <- rownames(sigs)
    return(new.counts)
}


# sigs <- read.COSMIC.sigs()  ## I x K
# counts <- read.csv("../../../BayesPowerSig/BPS-PCAWG-sigs/21_breast_cancers.mutations.txt", sep = "\t", row.names = 1)  ## I x J
# prelim.loadings <- loadings.nnls(counts, sigs)
# trimmed.loadings <- trim.loadings(prelim.loadings)
# simulate.counts(trimmed.loadings, sigs, "poisson")


