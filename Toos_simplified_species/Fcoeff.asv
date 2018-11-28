function [aE aW aP b Istart_F Y_sink] = Fcoeff(NPI, rho, A, x, x_u, u, Y_k, Y2_k,T, rho_real, rho2_real, P_k, D, relax_f, Dt, f_old, Dx, rho_old, sink, Pp, Pr, MW1, MW2, w);
    Istart_F = 2;
    F_u = conv(NPI, rho, x, x_u, u);
    global Runi
    
    for I=Istart_F:NPI+1
        i = I;
        
        Fw = F_u(i)*A;
        Fe = F_u(i+1)*A;
                
        % diffusivity D, at the INTERFACE is calculated 
        % with the use of a harmonic mean
%         Dw = ((rho(I-1)*D(I-1)*rho(I)*D(I))/(rho(I-1)*D(I-1)*(x(I)-x_u(i)) + rho(I)*D(I)*(x_u(i)-x(I-1))))*A;
%         De = ((rho(I)*D(I)*rho(I+1)*D(I+1))/(rho(I)*D(I)*(x(I+1)-x_u(i+1)) + rho(I+1)*D(I+1)*(x_u(i+1)-x(I))))*A;
%         
        Dw = ((rho(I-1)*D(I-1)*rho(I)*D(I))/(rho(I-1)*D(I-1)*(x(I)-x_u(i)) + rho(I)*D(I)*(x_u(i)-x(I-1))))*A;
        De = ((rho(I)*D(I)*rho(I+1)*D(I+1))/(rho(I)*D(I)*(x(I+1)-x_u(i+1)) + rho(I+1)*D(I+1)*(x_u(i+1)-x(I))))*A;
   
        
        Pe = Fw/Dw;
        % source terms
        if sink == 1
            SP(I) = 0;
            Su(I) = -P_k(I)*(Pr*MW1(I)*Y_k(I)-Pp*MW2(I)*Y2_k(I))*Dx*w*Dt;      % sink term depends on mass fraction
            Y_sink(I) = -P_k(I)*(Pr*MW1(I)*Y_k(I)-Pp*MW2(I)*Y2_k(I))*Dx*w*Dt;
        else 
        SP(I) = 0;
        Su(I) = 0;
        Y_sink(I) = 0;
        end
        
        aPold = rho_old(I)*Dx/Dt;
        
        if I == 2
            aW(i) = 0;
            aE(i) = max([-Fe De-Fe/2 0]);
            aP(i) = max([Fw Dw+Fw/2 0]) + aE(I) + Fe - Fw + aPold;
            Su(i) =  Su(i) + max([Fw Dw+Fw/2 0])*Y_k(I-1);
        elseif I == NPI +1
            aW(I) = max([Fw,  Dw+Fw/2, 0]);
            aE(I) = 0;
            Su(i)  = Su(i) + max([-Fe De-Fe/2 0])*Y_k(i+1); 
            aP(I) = aW(I) + max([-Fe De-Fe/2 0]) + aE(I) + Fe - Fw - SP(I) +aPold;
        else 
        % the coefficients (with hybrid differencing scheme)
        aW(I) = max([Fw,  Dw+Fw/2, 0]);
        aE(I) = max([-Fe, De-Fe/2, 0]);
        
        aP(I) = aW(I) + aE(I) + Fe - Fw - SP(I) +aPold;
        end
        b(I) = Su(I) + aPold*f_old(I);

        aP(I) = aP(I)/relax_f;              
       
        b(I) = b(I) + (1-relax_f)*aP(I)*Y_k(I);
    end
    
end