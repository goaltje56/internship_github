function [u u_guess T m_in m_out p] = bound(NPI,rho,x,x_u,A,u, u_in, T, u_guess, p)
    T(1) = 373;                     % temperature at inlet
    u(2) = u_in;                    % velocity at inlet
    u_guess(2) = u_in;
    
    F_u = conv(NPI, rho, x, x_u, u);
    m_in = F_u(2)*A; 
    m_out = F_u(NPI+1)*A;
    
    u(NPI+2) = u(NPI+1)*m_in/m_out; % velocity at outlet
    u_guess(NPI+2) = u_guess(NPI+1)*m_in/m_out; % velocity at outlet

%     pc(NPI+1) = 0;
%     p(NPI+1)  = 100;
    T(NPI+1) = T(NPI);
end