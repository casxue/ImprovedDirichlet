
library(JuliaCall)
# Setup Julia
julia <- julia_setup(installJulia = TRUE)

# Install packages
julia_install_package_if_needed("Distributions")
julia_install_package_if_needed("SpecialFunctions")
julia_install_package_if_needed("LinearAlgebra")
julia_install_package_if_needed("Printf")

# Call the file
julia_source("Julia/max_density.jl")

# Call the functions
max_density_for_beta <- function(c, v){
  out <- julia_call("max_density_for_beta", c, v)
  return(out)
}



