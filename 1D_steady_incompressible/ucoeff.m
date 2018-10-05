function [aP aE aW b d_u Istart_u u] = ucoeff(NPI, rho, x, x_u, u, p, A, relax_u, d_u,mu, u_in)
    Istart_u = 3;
    [u m_in m_out] = bound(NPI,rho,x,x_u,A,u, u_in);
    F_u = conv(NPI, rho, x, x_u, u);

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
        
%         if I == 3
%             aE(i) = max([-Fe De-Fe/2 0]);
%             aW(i) = 0;
%             aP(i) = aE(i) + aW(i) + Fe - Fw;         % but what if u<0 ??        
%             b(i)  = Fw*u(1);
                  
%         if I == NPI + 1
%             aE(i) = 0;
%             aW(i) = max([Fw Dw+Fw/2 0]);
%             aP(i) = Fe;         % but what if u<0 ??        
%             b(i)  = 0;
%             
%         else 
       
        % coefficients (hybrid differencing scheme)
        aW(i) = max([Fw Dw+Fw/2 0]);
        aE(i) = max([-Fe De-Fe/2 0]);
        b(i)  = 0;
        aP(i) = aW(i)+aE(i)+ Fe - Fw;
            
%         end
        % pressure correction 
        d_u(i) = A*relax_u/aP(i);
        
        % putting integrated pressure gradient in RHS
        % to solve with TDMA alogrithm
        b(i) = b(i) + (p(I-1) - p(I))*A ;
        
        % relaxation to aP and put last term on right side into source term
        aP(i) = aP(i)/relax_u;
        b(i)  = b(i) + (1-relax_u)*aP(i)*u(i);
        
    end
end