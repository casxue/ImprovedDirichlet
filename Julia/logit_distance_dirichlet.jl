# Distribution of distance between logit transformed values

using Distributions
using PyPlot
using SpecialFunctions
using Statistics
include("max_density.jl")


logistic(x) = 1/(1+exp(-x))
logit(p) = log(p/(1-p))


# tag = "nonuniform-small-k"
tag = "nonuniform-large-k"


if tag=="nonuniform-small-k"
    c = [0.01, 0.1, 0.2, 0.3, 0.39]
    N = 10000000
    width = 1
    xmax = 500
elseif tag=="nonuniform-large-k"
    # c = ones(30); c = c/sum(c)
    c = (1:30); c = c/sum(c)
    N = 10000000
    width = 20
    xmax = 5000
else
    @assert(false, "Unknown tag")
end
v = 0.1
alpha = 1.0


# Max density method
a1 = max_density_for_dirichlet(c,v; alpha=alpha, verbose=false)
p1 = rand(Dirichlet(a1),N)
d1 = [sum(abs.(logit.(p1[:,i]) .- logit.(c))) for i=1:N]


# Mean method
a2 = alpha*c
# p2 = max.(rand(Dirichlet(a2),N), 1e-300)
p2 = rand(Dirichlet(a2),N)
d2 = [sum(abs.(logit.(p2[:,i]) .- logit.(c))) for i=1:N]

# d2 = zeros(N)
# for i=1:N
    # x = rand.(Gamma.(a2,1))
    # s = sum(x)
    # logit_p = [log(xi/(s-xi)) for xi in x]
    # if any(isinf.(logit_p)); println("x = $x"); println("s = $s"); end
    # d2[i] = sum(abs.(logit_p .- logit.(c)))
# end






close(1)
figure(1,figsize=(5.5,3.2)); clf(); grid(lw=0.2)
subplots_adjust(bottom=0.2)

h1 = fit(Histogram, d1, [(0:width:xmax); Inf])
h2 = fit(Histogram, d2, [(0:width:xmax); Inf])
edges = h1.edges[1]
f1 = h1.weights/(N*width)
f2 = h2.weights/(N*width)

# hist(d1; bins=[(0:xmax); Inf], density=true, histtype="step", label="max density")
# hist(d2; bins=[(0:xmax); Inf], density=true, histtype="step", label="mean method")

# hist(d1; bins=50, density=true, histtype="step", label="max density")
# hist(d2; bins=50, density=true, histtype="step", label="mean method")

xs = (edges[1:end-1] .+ edges[2:end])/2
plot(xs, f1, "b-", lw=2, label="max density")
plot(xs, f2, "r-", lw=2, label="mean method")
xlim(0,xmax*0.6)
xlabel("distance from target location")
ylabel("frequency")
# title("c = $c")
legend()
savefig("logit_distance-$tag.png",dpi=200)




nothing
