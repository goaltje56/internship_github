function [x] = self_mass(M, x, sink, Xr, X, n)
% Mr    = x(1)
% Mp    = x(2)
% XAp   = x(3)
% XBp   = x(4)
% XCp   = x(5)
% XDp   = x(6)

%% these equations depend on the boundary conditions!!
Mr_new =0;

while abs(x(1) - Mr_new) > 0.0001
    A = ((Xr-x(3:n)) ~= 0);
    Mr_new = x(1);
    
    for i = 1:n      
        if A(i) == true
            x(1) = M*(X(i)-x(i+2))/(Xr(i)-x(i+2));
        end
    end

    
    x(3:n+2) = (M*X- x(1)*Xr)/(M-x(1));
    x(3:n+2) = x(3:n+2).*sink;
    x(3:n+2) = x(3:n+2)/sum(x(3:n+2));
    
end

Mp = M - x(1);
% x(1) = Mr;
x(2) = Mp;

end