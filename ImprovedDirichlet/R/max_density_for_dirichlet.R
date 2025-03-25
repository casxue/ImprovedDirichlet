max_density_for_dirichlet <- function (c, theta, version, a.init = NULL,
                                       tol = 1e-8, max.iter = 100, step.size = 0.5, verbose = TRUE)  {
    ## input protection
    if (version == "concentration")  {
        if (theta <= 0)  {
            stop(paste("Concentration parameter must satisfy theta > 0, but input value is theta =", theta))
        }
    } else if (version == "kl_uniform")  {
        if (theta <= 0)  {
            stop(paste("KL divergence must satisfy theta > 0, but input value is theta =", theta))
        }
    } else if (version == "cosine_error")  {
        if (theta <= 0 || theta >= 1)  {
            stop(paste("Cosine error must satisfy 0 < theta < 1, but input value is theta =", theta))
        }
    } else  {
        stop(paste("Unknown argument for version.  Got", version, "but expected one of 'concentration', 'kl_uniform', or 'cosine_error'."))
    }

    K <- length(c)
    m <- which.min(c)

    ## initialization
    if (is.null(a.init))  {
        if (version == "kl_uniform" || version == "cosine_error")  {
            a.init <- 5 * (c + 1)
        } else  {
            a.init <- c
        }
    }

    for (attempt in 1:5)  {
        a <- a.init

        for (iter in 1:max.iter)  {
            s <- sum(a)

            ## gradient and Hessian of objective function
        	g <- digamma(a) - digamma(sum(a)) - log(c)
            H <- diag(trigamma(a)) - trigamma(sum(a))

            ## constraint function and its Jacobian
            if (version == "concentration")  {
                h <- s / theta - 1
                J <- rep(1 / theta, K)
            } else if (version == "kl_uniform")  {
                ent <- sum(lgamma(a)) - lgamma(s) + (s - K) * digamma(s) - sum((a - 1) * digamma(a))
                h <- - ent - lgamma(K) - theta
                J <- - (s - K) * trigamma(s) + (a - 1) * trigamma(a)
            } else if (version == "cosine_error")  {
                s2 <- sum(a^2)
                s3 <- sum(a^3)
                h <- log(s) - log(2) - log(s + 1) - log(s2) + log(s - s3 / s2) - log(theta)
                J <- 1 / s - 1 / (s + 1) - 2 * a / s2 + (1 - (s2 * 3 * a^2 - s3 * 2 * a) / s2^2) / (s - s3 / s2)
            }

            ## solve for update
            A <- cbind(H, J) %>% rbind(c(J, 0))
            y <- c(-g, -h)
            x <- solve(A, y)

            if (verbose)  {
                f <- sum(lgamma(a)) - lgamma(s) - sum((a - 1) * log(c))
                cat(sprintf("iter = %d \t f = %.10f \t h = %.10f \t a = %s \n", iter, f, h, paste(signif(a, 3), collapse = " ")))
            }

            ## update values
            a.old <- a
            a <- a.old + step.size * head(x, -1)

            ## fix non-positive entries
            violations <- which(a <= 0)
            a[violations] <- a.old[violations] / 2

            ## check for convergence
            if (sum(abs(a / a.old - 1)) + abs(h) < tol)  return(a)
        }
        if (verbose && log10(sum(a)) - log10(tol) < 16)  print("Sum of a values may be too large to handle in floating-point precision.  Can try reducing tol to force convergence.")

        step.size <- step.size / 5
        max.iter <- max.iter * 5
        if (verbose)  cat(sprintf("Failed to converge. Retrying with step.size = %.5f and max.iter = %d", step.size, max.iter))
    }
    stop("Failed to converge after all attempts.")
}
