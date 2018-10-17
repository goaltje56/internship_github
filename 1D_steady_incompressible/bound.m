function [u T m_in m_out p] = bound(NPI,rho,x,x_u,A,u, u_in, T, p)
    T(1) = 273;                     % temperature at inlet
    u(2) = u_in;                    % velocity at inlet
    
    F_u = conv(NPI, rho, x, x_u, u);
    m_in = F_u(2)*A; 
    m_out = F_u(NPI+1)*A;
%     p(NPI+1) = p(NPI);
    u(NPI+2) = u(NPI+1)*m_in/m_out; % velocity at outlet
%     T(NPI+1) = T(NPI);
end