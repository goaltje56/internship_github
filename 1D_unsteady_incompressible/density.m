function rho = density(NPI, relax_rho, rho, p, P_atm, T)
    for I = 1:NPI+1
        if I == 1
            rho(I) = (1-relax_rho)*rho(I) + relax_rho*(p(I+1)+P_atm)/(287*T(I))
        else
            rho(I) = (1-relax_rho)*rho(I) + relax_rho*(p(I)+P_atm)/(287*T(I))
        end
    end
end