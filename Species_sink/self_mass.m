function [x] = self_mass(M, x, sink, Xr, X, n)
% Mr    = x(1)
% Mp    = x(2)
% XAp   = x(3)
% XBp   = x(4)
% XCp   = x(5)
% XDp   = x(6)

% M X = Mr Xr + Mp Xp
% Mr  = M (X - Xp)/ (Xr - Xp)
% Xp = (M X - Mr Xr) / (M-Mr)

%% these equations depend on the boundary conditions!!
Mr_new =0;

if n == 2
    x(1) = M * (X(1)-X(2)-1)/ (Xr(1)-Xr(2)-1);
    x(n+2) = (x(1)*Xr(1)-M*X(1))/(M-x(1))+1;
    x(n+1) = 1 - x(n+2);
else 

while abs(x(1) - Mr_new) > 0.0001
    A = ((Xr-x(3:n+2)') ~= 0);
    Mr_new = x(1);
    i = 1;
    if A(i) == true
            x(1) = M*(X(i)-x(i+2))/(Xr(i)-x(i+2));
    else
            i = i+1;
    end
    while x(1) > M
        x(1) = x(1)/1.2;
    end
    counter = 0;
    for j = 1:n
        if j ~= i
            x(j+2) = ((M*X(j)- x(1)*Xr(j))/(M-x(1)))*sink(j);
            counter = counter + x(j+2);
        end
    end
    x(3:n+2) = x(3:n+2).*sink;
    x(i+2)   = 1 - counter;
%     x(3:n+2) = x(3:n+2)/sum(x(3:n+2));
end

end
    x(2) =  M - x(1);
end