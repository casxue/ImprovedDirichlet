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
num_ks = 200
k_values = range(1.0, stop=200.0, length=num_ks)
# Create an empty matrix to store the results
results = zeros(Float64, size(X, 1), num_ks, size(X, 2))
for j in 1:size(X, 2)
    @printf("signature = %i \n", j)
    y = X[:, j]  # Get the column vector x
    # Loop over k from 1 to 200
    for (i, k) in enumerate(k_values)
        result = max_density_for_dirichlet(y, k, "kl_uniform")
        results[:, i, j] = result
    end
end

matwrite("data/alpha_kl_cosmic.mat", Dict("results" => results))



