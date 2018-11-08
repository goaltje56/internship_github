function [u Y_k p] = bound(NPI,rho,x,x_u,A,u, u_in, X_in, Y_k, p)
%     T(1) = T(2);                     % temperature at inlet
    u(2) = u_in;                    % velocity at inlet
    
%     Y_k(1,1) = 0.8;
%     Y_k(2,1) = 0.1;
%     Y_k(3,1) = 0.05;
%     Y_k(4,1) = 0.05;

    Y_k = MoleToMass(X_in, Y_k);

    F_u = conv(NPI, rho, x, x_u, u);
    m_in = F_u(2)*A; 
    m_out = F_u(NPI+1)*A;
    
    u(NPI+2) = u(NPI+1)*m_in/m_out; % velocity at outlet
%     T(NPI+2) = T(NPI+1);
%     p(NPI+2) = p(NPI+1);
%     pc(NPI+2) = pc(NPI+1);

%     T(NPI+2) = 0;
    Y_k(1, NPI+2) = Y_k(1, NPI+1);
    Y_k(2, NPI+2) = Y_k(2, NPI+1);
    Y_k(3, NPI+2) = Y_k(3, NPI+1);
    Y_k(4, NPI+2) = Y_k(4, NPI+1);

end