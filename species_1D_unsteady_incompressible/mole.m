function [X_k, rho, Gamma] = mole(NPI, n, Y_k, rho_k, MW, rho, Gamma, Gamma_k)
global Runiv
    for j = 1:n
        for I = 1:NPI+2
            if j == 1
                N(1,I) = 0;
                rho(1,I) = 0;
                Gamma(1,I) =0;
            end
            N_k(j,I) = Y_k(j,I)/MW(j);
            N(1,I)   = N(1,I)+N_k(j,I);
        end
    end
    
    for j = 1:n
        for I = 1:NPI+2
            X_k(j,I) = N_k(j,I)/N(1,I);
            rho(1,I) = rho(1,I) + rho_k(j,I)*X_k(j,I);
            Gamma(1,I) = Gamma(1,I) + Gamma_k(j)*X_k(j,I);
        end
    end
end
