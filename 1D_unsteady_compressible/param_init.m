function [u, p, pc, T, rho, mu, Cp, Gamma, d_u, b, SP, Su, relax_u, relax_pc, relax_T, relax_rho, Dt, u_old, T_old, pc_old] = param_init(NPI, u_in);
    m_in    = 1;
    m_out   = 1;
    Dt = 0.1;
    for I=1:NPI+2
        i = I;
        u(i)    = 1;         % Velocity in x-direction
        pc(I)   = 0;            % Pressure correction 
        p(I)    = 2;            % Relative pressure
        T(I)    = 273;          % Temperature
        rho(I)  = 1;            % Density
        mu(I)   = 2.5*10^(-5);  % Viscosity
        Cp(I)   = 1013;         % Heat capacity [J/(kg*K)]
        Gamma(I)= 0.5; %0.0315/Cp(I); % Thermal conductivity
        d_u(i)  = 0;            % Variable d[i] to calculate pc
        b(I)    = 0;            % The general constant
        SP(I)   = 0;            % Source term
        Su(I)   = 0;            % Source term
        u_old(i)= u(i);         % Velocity in x-direction old timestep
        pc_old(i)= u(i);        % Pressure correction old timestep
        T_old(i)= u(i);         % Temperature old timestep
    end
    u(NPI+1)   = 0.5*u_in;         % Velocity in x-direction
    relax_u  = 0.2;             % relaxation for velocity
    relax_pc = 1.1 - relax_u;   % relaxation for pressure
    relax_T = 1;              % relaxation for temperature
    relax_rho = 0.2;
end
