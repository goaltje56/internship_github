function [u T Y_k m_in m_out p] = bound(NPI,rho,x,x_u,A,u, u_in, T, Y_k, p)
    T(1) = 273;                     % temperature at inlet
    u(2) = u_in;                    % velocity at inlet
%     u_pseudo(2) = u_in;                    % velocity at inlet

    Y_k(1,1) = 0;
    Y_k(2,1) = 1;
    F_u = conv(NPI, rho, x, x_u, u);
    m_in = F_u(2)*A; 
    m_out = F_u(NPI+1)*A;
    
    u(NPI+2) = u(NPI+1)*m_in/m_out; % velocity at outlet
%     u_pseudo(NPI+2) = u_pseudo(NPI+1)*m_in/m_out; % velocity at outlet

%     T(NPI+2) = T(NPI+1);
%     p(NPI+1) = p(NPI);

    T(NPI+2) = 273;
    Y_k(1, NPI+2) = Y_k(1, NPI+1);
    Y_k(2, NPI+2) = Y_k(2, NPI+1);

end