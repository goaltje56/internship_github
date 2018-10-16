function [aE aW aP b Istart_F] = Fcoeff(NPI, rho, A, x, x_u, u, Y_k, D, relax_f, Dt, f_old, Dx);
    Istart_F = 2;
    F_u = conv(NPI, rho, x, x_u, u);
    
    for I=Istart_F:NPI+1
        i = I;
        
        Fw = F_u(i)*A;
        Fe = F_u(i+1)*A;
                
        % diffusivity D, at the INTERFACE is calculated 
        % with the use of a harmonic mean
        Dw = ((rho(I-1)*D(I-1)*rho(I)*D(I))/(rho(I-1)*D(I-1)*(x(I)-x_u(i)) + rho(I)*D(I)*(x_u(i)-x(I-1))))*A;
        De = ((rho(I)*D(I)*rho(I+1)*D(I+1))/(rho(I)*D(I)*(x(I+1)-x_u(i+1)) + rho(I+1)*D(I+1)*(x_u(i+1)-x(I))))*A;
        
        % source terms
        SP(I) = 0;
        Su(I) = 0;
        
        % the coefficients (with hybrid differencing scheme)
        aW(I) = max([Fw,  Dw+Fw/2, 0]);
        aE(I) = max([-Fe, De-Fe/2, 0]);
        aPold = rho(I)*Dx/Dt;
        
        aP(I) = aW(I) + aE(I) + Fe - Fw - SP(I) +aPold;
        b(I) = Su(I) + aPold*f_old(I);

        aP(I) = aP(I)/relax_f;              
       
        b(I) = b(I) + (1-relax_f)*aP(I)*Y_k(I);
    end
    
end