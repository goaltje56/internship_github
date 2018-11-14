function F_u = conv(NPI, rho, x, x_u, u)
    for I=2:NPI+2 % start at 3 beceause I-2 is outside domain for i=1
        i = I;
        F_u(i) = ((rho(I-1)*(x(I)-x_u(i)))+(rho(I)*(x_u(i)-x(I-1))))*u(i)/ (x(I)-x(I-1)); 
%         if I == 2
%             F_u(i) = ((((rho(I)+rho(I-1))/2)*u(i))+(((rho(I-1)))*u(i-1)))/2             
%         else
%             F_u(i) = ((((rho(I)+rho(I-1))/2)*u(i))+(((rho(I-1)+rho(I-2))/2)*u(i-1)))/2 
%         end
        %% density exact average velocity
%         F_u(i) = ((u(i)+u(i-1))/2)*rho(I); 
    end
end