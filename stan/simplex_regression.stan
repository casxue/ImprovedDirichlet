data {
    int<lower=1> S;                // Number of observations (samples)
    int<lower=1> I;                // Number of dimensions of y_s
    int<lower=1> J;                // Number of covariates
    matrix[I, J] X;                // Covariate matrix (I x J)
    simplex[I] y[S];               // Dirichlet-distributed observations
    vector[J] a;               // Prior parameter for Dirichlet distribution of pi
}

parameters {
    simplex[J] probs;                 // Covariate weights as a simplex
    real<lower=0> phi;             // Scale parameter for alpha
}

model {
    // Priors
    probs ~ dirichlet(a);  // Dirichlet prior for pi
    phi ~ gamma(1, 1);                // Gamma prior for phi
    vector[I] alpha;
    for (i in 1:I) {
      alpha[i] = phi * dot_product(probs, X[i, ]);  // alpha_i = phi * sum_j (pi_j * x_{i,j})
      }
    // Likelihood
    for (s in 1:S) {
        y[s] ~ dirichlet(alpha);   // Dirichlet likelihood for each observation
    }
}

