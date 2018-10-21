function [dudx ] = derivatives(NPI,u, x, x_u)
    for I = 2:NPI+1
        dudx(I) = (u(i+1)-u(i))/(x_u(i+1)-x_u(i));
end