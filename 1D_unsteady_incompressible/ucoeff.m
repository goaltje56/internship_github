function [aP aE aW b d_u Istart_u u T] = ucoeff(NPI, rho, x, x_u, u, p, A, relax_u, d_u,mu, u_in, T, Dt, u_old, Dx)
    Istart_u = 3;
    [u T m_in m_out] = bound(NPI,rho,x,x_u,A,u, u_in, T);
    F_u = conv(NPI, rho, x, x_u, u);
    fric_u = u_fric(NPI, mu, x, x_u, u);
    Dh = 0.001;
    for I = 3:NPI+1
        i = I;
           
        % convective mass flux eq 5.8a
%         Fw = (((rho(I-1)+rho(I))*u(i)/2)+((rho(I-1)+rho(I-2))*u(i-1)/2))*A/2;  % rho*u at west of cell face
%         Fe = (((rho(I+1)+rho(I))*u(i+1)/2)+((rho(I)+rho(I-1))*u(i)/2))*A/2;    % rho*u at east of cell face
        
        Fw = ((F_u(i)+F_u(i-1))/2)*A;
        Fe = ((F_u(i)+F_u(i+1))/2)*A;
        
        % transport by diffusion eq 5.8b
        Dw = (mu(I-1)/(x_u(i)-x_u(i-1)))*A;
        De = (mu(I)/(x_u(i+1)-x_u(i)))*A;

        Mw = ((fric_u(i)+fric_u(i-1))/2)/Dh;
        Me = ((fric_u(i)+fric_u(i+1))/2)/Dh;
        
        % coefficients (hybrid differencing scheme)
        aW(i) = max([Fw Dw+Fw/2 0]);
        aE(i) = max([-Fe De-Fe/2 0]);
        b(i)  = (Me+Mw)/2;
        
        aPold = 0.5*(rho(I-1)+rho(I))*Dx/Dt;
        % without time dependent terms 
        aP(i) = aW(i)+aE(i)+ Fe - Fw + aPold;
            
        % pressure correction 
        d_u(i) = A*relax_u/aP(i);
        
        % putting integrated pressure gradient in RHS
        % to solve with TDMA alogrithm
        b(i) = b(i) + (p(I-1) - p(I))*A +aPold*u_old(i);
        
        % relaxation to aP and put last term on right side into source term
        aP(i) = aP(i)/relax_u;
        b(i)  = b(i) + (1-relax_u)*aP(i)*u(i);
        
    end
end