function beta(rho::Vector{Float64})
    r = (
         rho.^2, # beta
         2*rho,  # beta prime
         2*ones(size(rho)) # beta second
        )
    return r
end
