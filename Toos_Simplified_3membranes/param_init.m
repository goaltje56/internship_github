function [u, u2, rho, d_u, b, SP, Su, relax_rho, relax_f, Dt, u_old] = param_init(NPI, u_in);
%     m_in    = 1;
%     m_out   = 1;
    Dt = 0.05;
    for I=1:NPI+2
        i = I;
        u(i)    = u_in;         % Velocity in x-direction
        u2(i)   = 0;
        rho(I)  = 2;            % Density
        d_u(i)  = 0;            % Variable d[i] to calculate pc
        b(I)    = 0;            % The general constant
        SP(I)   = 0;            % Source term
        Su(I)   = 0;            % Source term
        u_old(i)= u(i);         % Velocity in x-direction old timestep
        rho_old(i)= rho(i);         % Temperature old timestep

    end

    relax_rho = 0.5;
    relax_f = 1;              % relaxation for temperature

end
