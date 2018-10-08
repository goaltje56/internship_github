function [u T m_in m_out] = bound(NPI,rho,x,x_u,A,u, u_in, T)
    T(1) = 373;                     % temperature at inlet
    u(2) = u_in;                    % velocity at inlet
    p(2) = 200000;
    F_u = conv(NPI, rho, x, x_u, u);
    m_in = F_u(2)*A(1); 
    m_out = F_u(NPI+1)*A(end);
    
    u(NPI+2) = u(NPI+1)*m_in/m_out; % velocity at outlet
%     T(NPI+1) = 200;
    T(NPI+2) = T(NPI+1);
end