function [x Xr] = self_mass(M, x, sink, Xr, X, Y_sink, n)
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

x(2) = -sum(sum(Y_sink(1:end)));
x(1) = M - x(2);

for j = 1:n
    x(j+2) = 1;
    for i = 1:n
        if i ~= j
            x(j+2) = x(j+2) + sum(Y_sink(i,1:end))/x(2);
        end
    end
end

for i = 1:n
    Xr(i) = (M*X(i) -x(2)*x(i+2))/x(1);
end
a= M*X(1);
b= x(2)*x(3);
c = 0;
end