function [aP aE aW b d_u Istart_u u T] = ucoeff(NPI, rho, x, x_u, u, p, A, relax_u, d_u,mu, u_in, T, Dt, u_old, Dx)
    Istart_u = 3;
    [u T m_in m_out] = bound(NPI, rho, x, x_u, A, u, u_in, T);
    F_u = conv(NPI, rho, x, x_u, u);
    fric_u = u_source(NPI, mu, x, x_u, u);
    
    Dh = 0.0001;     % the Hydraulic diameter for pipe friction
    for I = 3:NPI+1
        i = I;
           
        % convective mass flux eq 5.8a
  
        Fw = ((F_u(i)+F_u(i-1))/2)*A(I);
        Fe = ((F_u(i)+F_u(i+1))/2)*A(I);
        
        Mw = ((fric_u(i)+fric_u(i-1))/2)/Dh;
        Me = ((fric_u(i)+fric_u(i+1))/2)/Dh;
        
        % transport by diffusion eq 5.8b
        Dw = (mu(I-1)/(x_u(i)-x_u(i-1)))*A(I);
        De = (mu(I)/(x_u(i+1)-x_u(i)))*A(I);
              
        % coefficients (hybrid differencing scheme)
        aW(I) = max([Fw Dw+Fw/2 0]);
        aE(I) = max([-Fe De-Fe/2 0]);
        b(I)  = (Mw + Me)/2;
        
        aPold = 0.5*(rho(I-1)+rho(I))*Dx/Dt;
        % without time dependent terms 
        aP(I) = aW(I)+aE(I)+ Fe - Fw + aPold;
            
        % pressure correction 
        d_u(i) = A(I)*relax_u/aP(I);
        
        % putting integrated pressure gradient in RHS
        % to solve with TDMA alogrithm
        b(i) = b(i) + (p(I-1) - p(I))*A(I) +aPold*u_old(i);
        
        % relaxation to aP and put last term on right side into source term
        aP(I) = aP(i)/relax_u;
        b(I)  = b(I) + (1-relax_u)*aP(I)*u(i);
        
    end
end