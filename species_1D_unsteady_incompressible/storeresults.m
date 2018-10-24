function [u_old pc_old rho_old f_old] = storeresults(NPI, u, pc, rho, Y_k, u_old, pc_old, rho_old, f_old, n)
    for I = 3:NPI+1
        i = I;
        u_old(i) = u(i);
    end
    
    for I = 2:NPI+1
        pc_old(I) = pc(I);
%         T_old(I)  = T(I);
        rho_old(I) = rho(I);
    end

    for j = 1:n
        for I = 2:NPI+1
            f_old(j,I) = Y_k(j,I);
        end
    end
end