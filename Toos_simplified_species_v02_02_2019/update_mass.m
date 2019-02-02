function [m_sink, m_out, m_in, m_system] = update_mass(NPI, mass, Y_sink,Y_k, u_in)
m_in = sum(mass);

   for i = 2:NPI+1
       m_sink(:,i)      = m_in*Y_sink(:,i)/u_in;
       m_out            = m_in + sum(m_sink(:,i));
       m_in             = m_out;
   end
  
   m_system    = m_in + m_sink;
   m_in = sum(mass);


end