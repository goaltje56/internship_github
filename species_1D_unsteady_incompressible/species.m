function [rho_k Y_k D_k p_k m_k M] = species(NPI, n, m, rho, p, D, MW)

M = 0;
    for j = 1:n
        
        for I = 1:NPI+2
            rho_k(j,I)= rho(j);
            p_k (j,I) = p(j);
            D_k(j,I)  = D(j);
            m_k(j,I)  = m(j);
            MW_k(j,I) = MW(j);
        end
        
        M = M + m_k(j,1);
    end
    
    for j = 1:n
        
        for I = 1:NPI+1 
            Y_k(j,I) = m_k(j,I)/M;
        end
        
    end
    
end