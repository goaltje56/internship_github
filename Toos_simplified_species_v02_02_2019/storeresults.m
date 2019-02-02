function [rho_old rho2_old f_old f2_old] = storeresults(NPI, rho, rho2, Y_k, Y2_k, rho_old, rho2_old, f_old, f2_old, n)
    
    for I = 2:NPI+2
        rho_old(I) = rho(I);
        rho2_old(I) = rho2(I);
    end

    for j = 1:n
        for I = 2:NPI+2
            f_old(j,I) = Y_k(j,I);
            f2_old(j,I) = Y2_k(j,I);
        end
    end
end