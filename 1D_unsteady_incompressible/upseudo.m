function [u_guess, d_u] = upseudo(NPI, rho, x, x_u, u, A, relax_u, d_u,mu, Dt, Dx, u_guess)
    Istart_u = 3;
    F_u = conv(NPI, rho, x, x_u, u);
    u_fric = fric_u(NPI, mu, x, x_u, u);
    Dh = 0.001;
    for I = 3:NPI+1
        i = I;
           
        Fw = ((F_u(i)+F_u(i-1))/2)*A;
        Fe = ((F_u(i)+F_u(i+1))/2)*A;
        
        % transport by diffusion eq 5.8b
        Dw = (mu(I-1)/(x_u(i)-x_u(i-1)))*A;
        De = (mu(I)/(x_u(i+1)-x_u(i)))*A;
        
        % Friction
        Mw = ((u_fric(i)+u_fric(i-1))/2)/Dh;
        Me = ((u_fric(i)+u_fric(i+1))/2)/Dh;
        
        % coefficients (hybrid differencing scheme)
        aW(i) = max([Fw Dw+Fw/2 0]);
        aE(i) = max([-Fe De-Fe/2 0]);
        b(i)  = (Mw + Me) /2;
        
        aPold = 0.5*(rho(I-1)+rho(I))*Dx/Dt;
        % without time dependent terms 
        aP(i) = aW(i)+aE(i)+ Fe - Fw + aPold;
            
        % pressure correction 
        d_u(i) = A*relax_u/aP(i);
        
        % guessed velocity
        u(i) = (aW(i)+aE(i)+b(i)) /aP(i);
        
    end
end