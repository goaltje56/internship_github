function [m_in, m_out] = globcont(NPI, rho,x, x_u, u, A)
F_u = conv(NPI, rho, x, x_u, u);
 
 m_in = F_u(2)*A;
 m_out = F_u(NPI+1)*A;
 
end