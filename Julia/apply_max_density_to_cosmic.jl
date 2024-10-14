# Apply the method max_density with kl to cosmic data, and save the output
# so we can read it in R

using Distributions
using SpecialFunctions
using LinearAlgebra
using Printf
using JLD2
using CSV
using DataFrames
using MAT

include("max_density.jl")

df = CSV.File("data/cosmic3.4.csv") |> DataFrame
X = Matrix(df[:, 2:end])

num_cols = size(X, 2)
num_ks = 500
#k_values = range(1.0, stop=300.0, length=num_ks)
## Create an empty matrix to store the results
#results_kl = zeros(Float64, size(X, 1), num_ks, size(X, 2))
#results_conc = zeros(Float64, size(X, 1), num_ks, size(X, 2))

# Evaluate the KL and the constrained concentration.
#for j in 1:size(X, 2)
#    @printf("signature = %i \n", j)
#    y = X[:, j]  # Get the column vector x
#    # Loop over k from 1 to 300 for values of theta
#    for (i, k) in enumerate(k_values)
#        # Constrain concentration
#        result = max_density_for_dirichlet(y, k, "concentration")
#        results_conc[:, i, j] = result
#        # Constrain KL
#        result = max_density_for_dirichlet(y, k, "kl_uniform")
#        results_kl[:, i, j] = result
#    end
#end

# Save the output
#matwrite("data/alpha_cosmic_concentration.mat", Dict("results_conc" => results_conc))
#matwrite("data/alpha_cosmic_KL.mat", Dict("results_kl" => results_kl))

# Evaluate the KL and the constrained concentration.
results_cosine = zeros(Float64, size(X, 1), num_ks, size(X, 2))
theta_values = range(0.001, 1, length = num_ks)
for j in 1:size(X, 2)
    @printf("signature = %i \n", j)
    y = X[:, j]  # Get the column vector x
    # Loop over k from values of theta
    for (i, theta) in enumerate(theta_values)
        # Constrain cosine error
        result = max_density_for_dirichlet(y, theta, "cosine_error")
        results_cosine[:, i, j] = result
    end
end

matwrite("data/alpha_cosmic_cosine.mat", Dict("results_cosine" => results_cosine))




