function [rho_k D_k Y_k p_k M rho_mix rho_old f_old Gamma_mix] = species(NPI, n, Y, rho, p, D, MW, Gamma)
global Runiv
M = 0;
N = 0;

    for j = 1:n
        
        for I = 1:NPI+2
            p_k (j,I)   = p(j);         % pressure of species k
            Y_k(j,I)    = Y(j);         % mass of species k
            N_k(j,I)    = Y_k(j,I)/MW(j);   % moles of species k
            MW_k(j,I)   = MW(j);        % molar weight of species k
            rho_k(j,I)  = rho(j);       % density of species k
            rho_mix(1,I)= 0;            % average density of mixture
            Gamma_mix(1,I) =0;
            Gamma_k(j, I) = Gamma(j);
            D_k(j,I)    = 0;         % Diffusivity of Species k
        end
        
        N = N + N_k(j,1);
    end
    
    for j = 1:n
        
        for I = 1:NPI+2 
            X_k(j,I) = N_k(j,I)/N;
            rho_mix(1,I) = rho_mix(1,I)+rho_k(j,I)*X_k(j,I);
            Gamma_mix(1,I) = Gamma_mix(1,I) + Gamma_k(j,I)*X_k(j,I);

        end
        
    end

    for j = 1:n
        for i = 1:NPI+2
            D_k(j,i) = D(j,:)*X_k(:,i);
        end
    end

    rho_old = rho_mix;
    f_old   = Y_k;
end