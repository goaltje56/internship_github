function [aE aW aP b Istart_pc pc] = pccoeff(NPI, rho, A, x, x_u, u, d_u, pc, rho_old, Dx, Dt)
    Istart_pc = 2;

    F_u = conv(NPI, rho, x, x_u, u);
    for I = Istart_pc:NPI
        i = I;
        
        % see eq. 6.32 
        
        SP(I) = 0;
        Su(I) = 0;
        
        % convective mass flux
        Fw = ((F_u(i)+F_u(i-1))/2)*A;
        Fe = ((F_u(i)+F_u(i+1))/2)*A;
        
        % the coefficients
        aE(I) = (rho(I+1)+rho(I))*d_u(i+1)*A/2;
        aW(I) = (rho(I)+rho(I-1))*d_u(i)*A/2;
        
        aP(I) = aE(I) + aW(I) - SP(I);
        b(I) = F_u(i)*A-F_u(i+1)*A + (rho_old(I)-rho(I))*A*Dx/Dt;
%         end
        
        pc(I) = 0;
        
    end
    
end