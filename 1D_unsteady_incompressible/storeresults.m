function [u_old pc_old T_old rho_old p_old] = storeresults(NPI, u, pc, p, T, rho, u_old, pc_old, T_old, rho_old, p_old)
    for I = 3:NPI+1
        i = I;
        u_old(i) = u(i);
    end
    
    for I = 2:NPI+1
        pc_old(I) = pc(I);
        p_old(I) = p(I);
        T_old(I)  = T(I);
        rho_old(I) = rho(I);
    end
end