function [rho_k Y_k D_k p_k m_k M rho_mix rho_old f_old] = species(NPI, n, m, rho, p, D, MW)
global Runiv
M = 0;
N = 0;
    for j = 1:n
        
        for I = 1:NPI+2
            p_k (j,I)   = p(j);         % pressure of species k
            D_k(j,I)    = D(j);         % Diffusivity of Species k
            m_k(j,I)    = m(j);         % mass of species k
            N_k(j,I)    = m(j)/MW(j);   % moles of species k
            MW_k(j,I)   = MW(j);        % molar weight of species k
            rho_k(j,I)  = rho(j);       % density of species k
            rho_mix(1,I)= 0;            % average density of mixture
        end
        
        M = M + m_k(j,1);
        N = N + N_k(j,1);
    end
    
    for j = 1:n
        
        for I = 1:NPI+2 
            Y_k(j,I) = m_k(j,I)/M;
            X_k(j,I) = N_k(j,I)/N;
            rho_mix(1,I) = rho_mix(1,I)+rho_k(j,I)*X_k(j,I);
        end
        
    end
    rho_old = rho_mix;
    f_old   = Y_k;
end