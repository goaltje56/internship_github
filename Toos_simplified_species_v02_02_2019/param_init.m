function [u, u2, relax_f, Dt] = param_init(NPI, u_in);
    Dt = 0.5;
    for I=1:NPI+2
        i = I;
        u(i)    = u_in;         % Velocity in x-direction
        u2(i)   = 0;
    end
    
    relax_f = 1;              % relaxation for temperature
end