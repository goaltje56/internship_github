function [aE aW aP b Istart_pc pc] = pccoeff(NPI, rho, A, x, x_u, u, d_u, pc)
    Istart_pc = 2;
%     pc(NPI+2) = 0; 
%     pc(1) = 0;
    F_u = conv(NPI, rho, x, x_u, u);
    for I = Istart_pc:NPI+1
        i = I;
        
%         Fw = (((rho(I-1)+rho(I))*u(i)/2)+((rho(I-1)+rho(I-2))*u(i-1)/2))*A/2;  % rho*u at west of cell face
%         Fe = (((rho(I+1)+rho(I))*u(i+1)/2)+((rho(I)+rho(I-1))*u(i)/2))*A/2;    % rho*u at east of cell face
        
%         Fw = (rho(I-1)*(u(i-1)+u(i))/2)*A;  % rho*u at west of cell face
%         Fe = (rho(I)*(u(i)+u(i+1))/2)*A;    % rho*u at east of cell face
               
        % see eq. 6.32 
        
        SP(I) = 0;
        Su(I) = 0;

%         if I == 2;
%             aW(i) = 0;
%             aE(I) = (rho(I+1)+rho(I))*d_u(i+1)*A/2;
%             aP(I) = aE(I) + aW(I) - SP(I);
%            
%         else
%         if I == NPI+2;
%             aE(i) = 0;
% %             pc(I+1) = pc(I);
%             aW(I) = (rho(I-1)+rho(I))*((d_u(i)+d_u(i-1))/2)*A/2;
%             aP(I) = aW(I) - SP(I);
%             b(I) = 0.5*(F_u(i-1)*A-F_u(i)*A);
%         else
            
        Fw = ((F_u(i)+F_u(i-1))/2)*A;
        Fe = ((F_u(i)+F_u(i+1))/2)*A;
        % the coefficients
        aE(I) = (rho(I+1)+rho(I))*d_u(i+1)*A/2;
        aW(I) = (rho(I)+rho(I-1))*d_u(i)*A/2;
        
        aP(I) = aE(I) + aW(I) - SP(I);
        b(I) = F_u(i)*A-F_u(i+1)*A;
%         end
        
        pc(I) = 0;
        
    end
    
end