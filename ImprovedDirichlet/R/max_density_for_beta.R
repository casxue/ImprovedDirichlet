max_density_for_beta <- function (c, v = NULL, alpha = NULL, a.init = NULL, b.init = NULL,
                                  tol = 1e-8, max.iter = 100, step.size = 0.5, verbose = FALSE)  {
    ## input protection
    if (is.null(v) && is.null(alpha))  {
        stop("Either target variance or concentration parameter must be specified")
    } else if (!is.null(v) && !is.null(alpha))  {
        stop("Cannot set target variance and concentration parameter simultaneously")
    }
    if (length(c) != 1 || c <= 0 || c >= 1)  {
        stop(paste("Target location must satisfy 0 < c < 1, but input value is c =", paste(c, collapse = " ")))
    }
    if (!is.null(v) && (v <= 0 || v >= 0.25))  {
        stop(paste("Variance must satisfy 0 < v < 0.25, but input value is v =", v))
    }
    if (!is.null(alpha) && alpha <= 0)  {
        stop(paste("Concentration parameter must satisfy theta > 0, but input value is alpha =", alpha))
    }

    for (attempt in 1:10)  {
        ## initialization
        a <- ifelse(is.null(a.init), c, a.init)
        b <- ifelse(is.null(b.init), 1 - c, b.init)


        for (iter in 1:max.iter)  {
            ## gradient of objective function
            df.a <- digamma(a) - digamma(a + b) - log(c)
            df.b <- digamma(b) - digamma(a + b) + log(1 - c)

            ## Hessian of objective function
                        df.aa <- trigamma(a) - trigamma(a + b)
            df.bb <- trigamma(b) - trigamma(a + b)
            df.ab <- - trigamma(a + b)
            H <- matrix(c(df.aa, df.ab, df.ab, df.bb), nrow = 2)

            ## constraint function and its Jacobian
            if (is.null(alpha))  {  ## variance constraint
                h <- log(a) + log(b) - 2 * log(a + b) - log(a + b + 1) - log(v)
                dh.a <- 1 / a - 2 / (a + b) - 1 / (a + b + 1)
                dh.b <- 1 / b - 2 / (a + b) - 1 / (a + b + 1)
            } else  {  ## concentration constraint
                h <- (a + b) / alpha - 1
                dh.a <- dh.b <- 1 / alpha
            }
            J <- c(dh.a, dh.b)

            ## solve for update
            A <- cbind(H, J) %>% rbind(c(J, 0))
            y <- c(-df.a, -df.b, -h)
            x <- solve(A, y)

            if (verbose)  {
                f <- - (a - 1) * log(c) - (b - 1) * log(1 - c) - lgamma(a + b) + lgamma(a) + lgamma(b)
                cat(sprintf("iter = %d \t a = %.10f \t b = %.10f \t f = %.10f \t h = %.10f\n", iter, a, b, f, h))
            }

            ## update values
            a.old <- a
            b.old <- b
            a <- a.old + step.size * x[1]
            b <- b.old + step.size * x[2]

            ## fix non-positive entries
            if (a < 0)  a <- a.old / 2
            if (b < 0)  b <- b.old / 2

            ## check for convergence
            if (abs(a / a.old - 1) + abs(b / b.old - 1) + abs(h) < tol)  return(c(a, b))
        }
        step.size <- step.size / 5
        max.iter <- max.iter * 5
        if (verbose)  cat(sprintf("Failed to converge. Retrying with step.size = %.5f and max.iter = %d", step.size, max.iter))
    }
    stop("Failed to converge after all attempts.")
}
