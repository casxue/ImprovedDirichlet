# Density max procedure for selecting Beta parameters


using Distributions
using SpecialFunctions
using LinearAlgebra
using Printf



# _____________________________________________________________________________
# Density max procedure for selecting Beta parameters, using Newton's method for constrained optimization.
#
# Inputs:
#   c = target location between 0 and 1
#   v = variance between 0 and 0.25
#
# Optional inputs:
#   a_init = initial a for optimization algorithm
#   b_init = initial b for optimization algorithm
#   tol = convergence tolerance for optimization algorithm
#   alpha = concentration parameter  (If specified, this overrides the variance v.)
#
# Outputs:
#   a,b = Beta parameters
#
# Reference: This is based on the general algorithm for Newton's method with equality constraints described here:
#   https://www.cs.cmu.edu/~ggordon/10725-F12/scribes/10725_Lecture12.pdf
#
function max_density_for_beta(c, v; a_init=c, b_init=1-c, tol=1e-8, maxiter=100, stepsize=0.5, verbose=true, alpha=NaN)
    @assert(0 < c < 1, "Target location must satisfy 0 < c < 1, but input value is c = $c.")
    @assert(0 < v < 0.25, "Variance must satisfy 0 < v < 0.25, but input value is v = $v.")
    if !isnan(alpha); @assert(alpha > 0, "Concentration parameter must satisfy alpha > 0, but input value is alpha = $alpha."); end

    for attempt = 1:10
        a = a_init
        b = b_init
        a_old = a
        b_old = b
        
        for iter = 1:maxiter
            
            df_a = -log(c) - digamma(a+b) + digamma(a)
            df_b = -log(1-c) - digamma(a+b) + digamma(b)
            g = [df_a ; df_b]
            
            if isnan(alpha)
                h = log(a) + log(b) - 2*log(a+b) - log(a+b+1) - log(v)
                dh_a = 1/a - 2/(a+b) - 1/(a+b+1)
                dh_b = 1/b - 2/(a+b) - 1/(a+b+1)
            else
                h = (a + b)/alpha - 1
                dh_a = 1/alpha
                dh_b = 1/alpha
            end
            
            J = [dh_a dh_b]
            
            df_aa = -trigamma(a+b) + trigamma(a)
            df_bb = -trigamma(a+b) + trigamma(b)
            df_ab = -trigamma(a+b)
            
            H = [df_aa df_ab ; df_ab df_bb] 
            
            A = [H J' ; J 0]
            y = [-g ; -h]
            x = A\y
            
            if verbose
                f = -(a-1)*log(c) - (b-1)*log(1-c) - lgamma(a+b) + lgamma(a) + lgamma(b)  # only used for printing the value of the objective function
                @printf("iter = %i     a = %.10f      b = %.10f     f = %.10f    h = %.10f\n", iter, a, b, f, h)
            end
            
            a_old = a
            b_old = b
            
            a = a_old + stepsize*x[1]
            b = b_old + stepsize*x[2]
            
            if a<0; a = a_old/2; end
            if b<0; b = b_old/2; end
            
            if abs(a/a_old - 1) + abs(b/b_old - 1) + abs(h) < tol; return a,b; end
        end
        
        if verbose; println("Failed to converge. Retrying with stepsize=$stepsize and maxiter=$maxiter."); end
        stepsize = stepsize/5
        maxiter = maxiter*5
    end
    
    @assert(false, "Failed to converge after all attempts.")
    
    # return NaN,NaN
end


lgamma_stirling(x) = 0.5*log(2*pi/x) + x*(log(x) - 1)

lbeta_stirling(a) = (k=length(a); s=sum(a); 0.5*(k-1)*log(2*pi) + 0.5*log(s) - 0.5*sum(log.(a)) + sum(a.*log.(a/s)))



# _____________________________________________________________________________
# Density max procedure for selecting Dirichlet parameters, using Newton's method for constrained optimization.
#
# Inputs:
#   c = (length k vector) target location such that sum(c) = 1 and all(c > 0)
#   theta = interpretation depends on version
#   version = meaning of theta, can be either:
#       "concentration": theta = concentration parameter
#       "kl_uniform":    theta = Kullback-Leibler divergence from the uniform distribution on the simplex
#       "cosine_error":  theta = mean cosine error between a sample from the Dirichlet and its mean
#
# Optional inputs:
#   a_init = initial a for optimization algorithm
#   tol = convergence tolerance for optimization algorithm
#
# Outputs:
#   a = (length k vector) Dirichlet parameters
#
# Reference: This is based on the general algorithm for Newton's method with equality constraints described here:
#   https://www.cs.cmu.edu/~ggordon/10725-F12/scribes/10725_Lecture12.pdf
#
function max_density_for_dirichlet(c, theta, version; a_init=Float64[], tol=1e-8, maxiter=100, stepsize=0.5, verbose=false)
    @assert((abs(sum(c)-1)<1e-10) && all(c.>0), "Target location must satisfy sum(c)=1 and all(c>0), but input value is c = $c.")
    if version=="concentration"
        @assert(theta > 0, "Concentration parameter must satisfy theta > 0, but input value is theta = $theta.")
    elseif version=="kl_uniform"
        @assert(theta > 0, "KL divergence must satisfy theta > 0, but input value is theta = $theta.")
    elseif version=="cosine_error"
        #@assert(0 < theta < 1, "Cosine error must satisfy 0 < theta < 1, but input value is theta = $theta.")
    else
        @assert(false, "Unknown argument for version.  Got $version but expected one of 'concentration', 'kl_uniform', or 'cosine_error'.")
    end

    k = length(c)
    m = argmin(c)
    
    if isempty(a_init)
        if version=="kl_uniform"
            # a_init=10*exp(theta/2)*c  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            # a_init=10*exp(theta/2)*(c + ones(k))/2
            a_init=10*(c + ones(k))/2
            # a_init=10*c
        elseif version=="cosine_error"
            a_init=10*(c + ones(k))/2
        else
            a_init = c
        end
    end
    
    for attempt = 1:5
        a = a_init
        a_old = a
        
        for iter = 1:maxiter
            
            # gradient of objective function
            g = digamma.(a) .- digamma(sum(a)) .- log.(c)
            
            # Hessian of objective function
            H = diagm(trigamma.(a)) .- trigamma(sum(a))
            
            if version == "concentration"
                h = sum(a)/theta - 1        # Constraint function
                J = ones(1,k)*(1/theta)     # Jacobian of constraint function
                
            elseif version == "kl_uniform"
                s = sum(a)
                ent = sum(lgamma.(a)) - lgamma(s) + (s - k)*digamma(s) - sum((a .- 1).*digamma.(a))
                # ent = lbeta_stirling(a) + (s - k)*digamma(s) - sum((a .- 1).*digamma.(a))
                h = -ent - lgamma(k) - theta
                J = -(s-k)*trigamma(s) .+ (a.-1).*trigamma.(a)
                J = vec(J)'
            
            elseif version == "cosine_error"
                s = sum(a)
                s2 = sum(a.^2)
                s3 = sum(a.^3)
                h = -log(2) + log(s) - log(s+1) - log(s2) + log(s - s3/s2) - log(theta)
                J = 1/s .- 1/(s+1) .- 2*a/s2 .+ (1 .- (s2*3*(a.^2) .- s3*2*a)/s2^2) / (s - s3/s2)
                J = vec(J)'
                
                # h = log(sum(a.*c)) - log(s) - log(sum(c.*c)) - log(1-theta)
                # J = c/sum(a.*c) .- 1/s
                # h = log(sum(a.*c)) - 0.5*log(sum(a.*a)) - 0.5*log(sum(c.*c)) - log(1-theta)
                # J = c/sum(a.*c) .- a/sum(a.*a)
                # h = (1 - sum((a/s).*c)/sum(c.*c))/theta - 1
                # J = -(s*c)/sum(c.*c) .- sum((-a/s^2).*c)/sum(c.*c)
                # h = log(sum(a.*c)) - log(sum(c.*c)) - log(1-theta)
                # J = c/sum(a.*c)
                # J = vec(J)'
            end
            
            A = [H J' ; J 0]
            y = [-g ; -h]
            x = A\y
            
            if verbose
                f = sum(lgamma.(a)) - lgamma(sum(a)) - sum((a.-1).*log.(c))  # only used for printing the value of the objective function
                @printf("iter = %i     f = %.10f    h = %.10f\n", iter, f, h)
                println("     a = ",round.(a; digits=3))
                if isnan(f); @assert(false); end
                # readline()
            end
            
            a_old = a
            
            a = a_old + stepsize*x[1:end-1]
            
            violations = (a .< 0)
            a[violations] = a_old[violations]/2
            
            if sum(abs.(a./a_old .- 1)) + abs(h) < tol; return a; end
        end
        
        if verbose && (log10(sum(a)) - log10(tol) < 16); println("Sum of a values may be too large to handle in floating-point precision.  Can try reducing tol to force convergence."); end
        if verbose; println("Failed to converge. Retrying with stepsize=$stepsize and maxiter=$maxiter."); end
        stepsize = stepsize/5
        maxiter = maxiter*5
    end
    
    
    @assert(false, "Failed to converge after all attempts.")
    
    # return NaN,NaN
end




