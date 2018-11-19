function [X_k,  D_k, rho_real, rho2_real] = mole(NPI, n, Y_k, Y2_k, MW,  D, rho_s)
global Runiv
    for j = 1:n
        for I = 1:NPI+2
            
            if j == 1
                N(1,I) = 0;
            end
            
            D_k(j,I) = 0;
            N_k(j,I) = Y_k(j,I)/MW(j);
            N(1,I)   = N(1,I)+N_k(j,I);
            
        end
    end
    
    for j = 1:n
        for I = 1:NPI+2
            X_k(j,I) = N_k(j,I)/N(1,I);
        end
    end

    for j = 1:n
        for i = 1:NPI+1
            D_k(j,i) = D(j,:)*X_k(:,i);
            rho_real(1,i) = rho_s*Y_k(:,i);
            rho2_real(1,i) = rho_s*Y2_k(:,i);            
        end
    end
end
