function [u_old pc_old T_old] = storeresults(NPI, u, pc, T, u_old, pc_old, T_old)
    for I = 3:NPI+1
        i = I;
        u_old(i) = u(i);
    end
    
    for I = 2:NPI+1
        pc_old(I) = pc(I);
        T_old(I)  = T(I);
    end
end