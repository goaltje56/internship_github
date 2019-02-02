function [Y_k, Y2_k, rho_real, rho2_real]  = bound(NPI,Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real)
        
    Y_k(:,1) = Y_in;
    Y2_k(:,1) = Y2_in;
    Y2_k(:,2) = Y2_k(:,3);
    
    rho_real(NPI+2) = rho_real(NPI+1);
    rho2_real(NPI+2) = rho2_real(NPI+1);
    
    Y_k(1, NPI+2) = Y_k(1, NPI+1);
    Y_k(2, NPI+2) = Y_k(2, NPI+1);
    Y_k(3, NPI+2) = Y_k(3, NPI+1);
    Y_k(4, NPI+2) = Y_k(4, NPI+1);

    Y2_k(1, NPI+2) = Y2_k(1, NPI+1);
    Y2_k(2, NPI+2) = Y2_k(2, NPI+1);
    Y2_k(3, NPI+2) = Y2_k(3, NPI+1);
    Y2_k(4, NPI+2) = Y2_k(4, NPI+1);
end