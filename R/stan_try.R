

library(rstan)



X <- sigs_all


# Compile the model
model_file <- "stan/simplex_regression.stan"  
fit <- stan(file = model_file, data = stan_data, iter = 2000, chains = 4)

# Check results
print(fit)



