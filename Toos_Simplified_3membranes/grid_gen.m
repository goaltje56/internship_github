function [NPI, x, x_u] = grid_gen(Dx, XMAX)
    % Length of the element
    NPI  = XMAX/Dx;
    x(1) = 0;
    x(2) = 0.5*Dx;
    
    % Length variable for the scalar points in the x direction
    for I=3:NPI+1
        x(I) = x(I-1)+Dx;
    end
    x(NPI+2) = x(NPI+1) + 0.5*Dx;
%     x(NPI+3) = x(NPI+1) + Dx;
    
    % Length variable for the velocity components u[i] in the x direction
    % Used to make staggered grid
    x_u(1) = 0;
    x_u(2) = 0;
    for i=3:NPI+2
        x_u(i) = x_u(i-1)+Dx;
    end   
%     x_u(NPI+1) = x_u(NPI)+Dx; 
end