function fric_u = u_fric(NPI, mu, x, x_u, u)

    for I=2:NPI+2 
        i = I;
        fric_u(i) = -((mu(I-1)*(x(I)-x_u(i)))+(mu(I)*(x_u(i)-x(I-1))))*u(i)/ (x(I)-x(I-1)); 
    end
    
end