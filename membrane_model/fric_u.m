function F_u = fric_u(NPI, mu, x, x_u, u)

    for I=2:NPI+2 % start at 3 beceause I-2 is outside domain for i=1
        i = I;
        F_u(i) = -((mu(I-1)*(x(I)-x_u(i)))+(mu(I)*(x_u(i)-x(I-1))))*u(i)/ (x(I)-x(I-1)); 
    end
    
end