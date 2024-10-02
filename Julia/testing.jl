# Testing the maximum density method

include("max_density.jl")

# _____________________________________________________________________________
# Example usage of Beta version a specific c and v value

c = 0.2  # target location
v = 0.05   # variance

a,b = max_density_for_beta(c,v)


V(a,b) = a*b/((a+b)^2 * (a+b+1))

println("v =      $v")
println("V(a,b) = ", V(a,b))
println("V(a,b)/v = ", V(a,b)/v)

println("u =      $(a/(a+b))")



using PyPlot
figure(1); clf(); grid(lw=0.2)
xs = 0.001:0.001:0.999
plot(xs, pdf.(Beta(a,b),xs), "b-")
yl,yu = ylim()
plot([c;c], [0;yu], "k--", lw=2)
ylim(yl,yu)


# _____________________________________________________________________________
# Example usage of Dirichlet version a specific c and v value

# c = 0.2*ones(3)  # target location
c = [0.001,0.1,0.2]  # target location
alpha = 1.0
v = 0.05   # variance

a = max_density_for_dirichlet(c,v; alpha=alpha)


# _____________________________________________________________________________
# Testing


function test_convergence(n_reps,max_v)
    @assert(max_v <= 0.25)
    n_failures = 0
    for rep = 1:n_reps
        if mod(rep,100)==0; println(rep); end
        
        c = rand()
        v = rand()*max_v
        
        a,b = max_density_for_beta(c,v; stepsize=0.1, maxiter=1000, verbose=false)
        
        fail = (abs(V(a,b)/v - 1) > 1e-6)
        
        n_failures += fail
    end
    
    return n_failures,n_reps
end



function test_convergence_alpha(n_reps)
    n_failures = 0
    for rep = 1:n_reps
        if mod(rep,100)==0; println(rep); end
        
        c = rand()
        alpha = exp(randn()*6 + 7)
        
        a,b = max_density_for_beta(c,0.1; stepsize=0.1, maxiter=1000, verbose=false, alpha=alpha)
        
        fail = (abs((a+b)/alpha - 1) > 1e-6)
        
        n_failures += fail
    end
    
    return n_failures,n_reps
end



function test_convergence_dirichlet(n_reps,k,max_v)
    @assert(max_v <= 0.25)
    n_failures = 0
    for rep = 1:n_reps
        if mod(rep,100)==0; println(rep); end

        c = rand(Dirichlet(k,0.5))
        v = rand()*max_v
        
        a = max_density_for_dirichlet(c,v; verbose=false)
        
        m = argmin(c)
        s = sum(a)
        fail = (abs(V(a[m],s-a[m])/v - 1) > 1e-6)
        
        n_failures += fail
    end
    
    return n_failures,n_reps
end



function test_convergence_dirichlet_alpha(n_reps,k)
    n_failures = 0
    for rep = 1:n_reps
        if mod(rep,100)==0; println(rep); end

        c = rand(Dirichlet(k,0.5))
        alpha = exp(randn()*6 + 7)
        
        a = max_density_for_dirichlet(c,0.1; verbose=false, alpha=alpha)
        
        fail = (abs(sum(a)/alpha - 1) > 1e-6)
        
        n_failures += fail
    end
    
    return n_failures,n_reps
end



# println("\nTesting with v specified ...")
# n_failures,n_reps = test_convergence(1000,0.25)
# println("Number of failures in testing: $n_failures / $n_reps")



# println("\nTesting with alpha specified ...")
# n_failures,n_reps = test_convergence_alpha(1000)
# println("Number of failures in testing: $n_failures / $n_reps")


println("\nTesting Dirichlet with variance specified ...")
n_failures,n_reps = test_convergence_dirichlet(1000,5,0.2)
println("Number of failures in testing: $n_failures / $n_reps")


println("\nTesting Dirichlet with alpha specified ...")
n_failures,n_reps = test_convergence_dirichlet_alpha(1000,5)
println("Number of failures in testing: $n_failures / $n_reps")




