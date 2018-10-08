function [aE aW aP b Istart_T] = Tcoeff(NPI, rho, A, x, x_u, u, T, Gamma, relax_T, Dt, T_old);
    Istart_T = 2;
    F_u = conv(NPI, rho, x, x_u, u);
    
    for I=Istart_T:NPI+1
        i = I;
        
        Fw = F_u(i)*A(I);
        Fe = F_u(i+1)*A(I);
                
        % conductivity, Gamma, at the INTERFACE is calculated 
        % with the use of a harmonic mean
        Dw = ((Gamma(I-1)*Gamma(I))/(Gamma(I-1)*(x(I)-x_u(i)) + Gamma(I)*(x_u(i)-x(I-1))))*A(I);
        De = ((Gamma(I)*Gamma(I+1))/(Gamma(I)*(x(I+1)-x_u(i+1)) + Gamma(I+1)*(x_u(i+1)-x(I))))*A(I);
        
        % source terms
        SP(I) = 0;
        Su(I) = 0;
        
        % the coefficients (with hybrid differencing scheme)
        aW(I) = max([Fw,  Dw+Fw/2, 0]);
        aE(I) = max([-Fe, De-Fe/2, 0]);
        aPold = rho(I)*A/Dt;
        
        aP(I) = aW(I) + aE(I) + Fe - Fw - SP(I) +aPold;
        b(I) = Su(I) + aPold*T_old(I);

        aP(I) = aP(I)/relax_T;
        
       
       
        b(I) = b(I) + (1-relax_T)*aP(I)*T(I);
    end
    
end