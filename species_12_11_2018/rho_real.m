function [m_in m_out m_sink] = rho_real(NPI, n, Y_k, rho_s, MW,  f_old)
global Runiv
    for j = 1:n
        for I = 1:NPI+2
            if j == 1
                N(1,I) = 0;
                Nold(1,I) = 0;
                rho_TRUE(1,I) = 0;
                rho_TRUEold(1,I) = 0;
                
            end
            N_k(j,I) = Y_k(j,I)/MW(j);
            N_kold(j,I) = f_old(j,I)/MW(j);
            N(1,I)   = N(1,I)+N_k(j,I);
            Nold(1,I) = Nold(1,I)+N_kold(j,I);
        end
    end
    
    for j = 1:n
        for I = 1:NPI+2
            X_k(j,I) = N_k(j,I)/N(1,I);
            X_kold(j,I) = N_kold(j,I)/Nold(1,I);
            rho_TRUE(1,I) = rho_TRUE(1,I) + rho_s(j)*X_k(j,I);
            rho_TRUEold(1,I) = rho_TRUEold(1,I) + rho_s(j)*X_kold(j,I);
        end
    end
    
    m_in = rho_TRUEold;
    m_out = rho_TRUE;
    m_sink = rho_TRUEold.*(f_old(1,:)-Y_k(1,:));       %% note, this is case specific!!
end
