# Distribution of distance between logit transformed values

using Distributions
using PyPlot
using SpecialFunctions
include("max_density.jl")
include("median_matching.jl")


logsumexp(a,b) = (m=max(a,b); (m==Inf ? Inf : (m==-Inf ? -Inf : log(exp(a-m)+exp(b-m))+m)))


f(x,a,b) = pdf(Beta(a,b),x)

logistic(x) = 1/(1+exp(-x))

dlogistic(x) = exp(-x)/(1+exp(-x))^2

g(y,a,b,c) = (l=log(c/(1-c)); f(logistic(l-y),a,b)*dlogistic(l-y) + f(logistic(l+y),a,b)*dlogistic(l+y))

# log_f(x,a,b) = logpdf(Beta(a,b),x)


log_logistic(x) = (x > -700 ? -log1p(exp(-x)) : x)

log_1mlogistic(x) = (x < 700 ? -log1p(exp(x)) : -x)

log_f(d,a,b) = -logbeta(a,b) + (a-1)*log_logistic(d) + (b-1)*log_1mlogistic(d)

log_dlogistic(x) = -x + 2*log_logistic(x)

function log_g(y,a,b,c)
    l = log(c) - log(1-c)
    
    t1 = log_f(l-y,a,b) + log_dlogistic(l-y)
    t2 = log_f(l+y,a,b) + log_dlogistic(l+y)

    # t1 = log_f(logistic(l-y),a,b) + log_dlogistic(l-y)
    # t2 = log_f(logistic(l+y),a,b) + log_dlogistic(l+y)
    
    return logsumexp(t1,t2)
end
    


if false
# ys = 0.1:0.1:1000
ys = 1:1:100000


alpha = 0.1
c = 0.001 # 0.001
u = c
a = alpha*u
b = alpha*(1-u)


# g is the analytical formula for the density of Y = |logit(X) - logit(c)| where X~Beta(a,b).
# Compare with histograms from sampling

figure(1); clf(); grid(lw=0.2)
# n = 100000
# Y_samples = [abs(log(X/(1-X)) - log(c/(1-c))) for X in rand(Beta(a,b),n)]
# hist(Y_samples,bins=100,density=true)
plot(ys, exp.(log_g.(ys,a,b,c)), "m-")
# plot(ys, g.(ys,a,b,c), "r-")

end




# _____________________________________________________________________________
# Plots

cs = [0.001,0.2]
alphas = 10.0 .^ (1:-1:-1)

include_median = false

# ys = 1:1:100000
ys = 0.01:0.01:10

colors = ["g","r","y"]  #["y","r","g","b","c","v"]

close(1)

for (i_c,c) in enumerate(cs)

    figure(1, figsize=(5.5,3.2)); clf(); grid(lw=0.2)
    
    for (i_alpha,alpha) in enumerate(alphas)
        
        # Mean method
        u = c
        a1 = alpha*u
        b1 = alpha*(1-u)
        
        # Maximum density method
        # v = V(a1,b1)  # choose variance to match the mean method
        # a2,b2 = max_density_for_beta(c,v)
        a2,b2 = max_density_for_beta(c,0.1; alpha=alpha)
        println("alpha = $alpha    a2+b2 = $(a2+b2)")
        
        # Median matching method
        u = match_median(c,alpha)
        a3 = alpha*u
        b3 = alpha*(1-u)
        println("p = $(cdf(Beta(a3,b3),c))")
        
        plot(ys, exp.(log_g.(ys,a1,b1,c)), "$(colors[i_alpha])-", label="alpha=$alpha") #alpha=10^$(round(Int,log10(alpha)))")
        plot(ys, exp.(log_g.(ys,a2,b2,c)), "$(colors[i_alpha])--") #, label="max density, alpha=$alpha")
        if include_median
        plot(ys, exp.(log_g.(ys,a3,b3,c)), "$(colors[i_alpha]):") #, label="median method, alpha=$alpha")
        end
        
    end
    
    title("c = $c")
    # legend(fontsize=9, borderaxespad=0.1)
    legend()
    savefig("logit_distance-c=$c-with-median=$(include_median).png",dpi=200)
end


