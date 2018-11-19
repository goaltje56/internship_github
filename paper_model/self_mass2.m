function [x Y2] = self_mass2(M, M2, x, sink, Y_sink, Y2_k, n)
% Mr        = x(1)
% Msink     = x(2)
% MAsink    = x(3)
% MBsink    = x(4)
% MCsink    = x(5)
% MDsink    = x(6)
% Mp        = x(7)

% M X = Mr Xr + Mp Xp
% Mr  = M (X - Xp)/ (Xr - Xp)
% Xp = (M X - Mr Xr) / (M-Mr)

%% these equations depend on the boundary conditions!!
Mr_new =0;

x(2) = -sum(sum(Y_sink(1:end)));            % x(2) = positive if sink
x(7) = M2 - x(2);
x(1) = M - x(2);

for j = 1:n
    x(j+2) = 0;
    x(j+2) = x(j+2) - sum(Y_sink(j,1:end));
end

for i = 1:n
    Y2(i) = (M2*Y2_k(i) + x(i+2))/x(7);
end

end