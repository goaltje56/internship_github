function F_u = conv(NPI, rho, x, x_u, u)

    for I=2:NPI+2 
        i = I;
        F_u(i) = ((rho(I-1)*(x(I)-x_u(i)))+(rho(I)*(x_u(i)-x(I-1))))*u(i)/ (x(I)-x(I-1)); 
    end
    
end