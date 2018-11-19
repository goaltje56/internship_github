function [x Y2] = self_mass2(M, M2, x, sink, Y_sink, Y2_k, MW, n)
% Mr        = x(1)      TOTAL mass at retenate
% Msink     = x(2)      TOTAL mass that permeates
% Mp        = x(3)      TOTAL mass at permeate side
% MAsink    = x(4)      Permeating mass of species A
% MBsink    = x(5)      Permeating mass of species B
% MCsink    = x(6)      Permeating mass of species C
% MDsink    = x(7)      Permeating mass of species D
% YAsink    = x(8)      Mass fraction of permeated species
% YBsink    = x(9)
% YCsink    = x(10)
% YDsink    = x(11)

% M X = Mr Xr + Mp Xp
% Mr  = M (X - Xp)/ (Xr - Xp)
% Xp = (M X - Mr Xr) / (M-Mr)

%% these equations depend on the boundary conditions!!
Mr_new =0;

x(2) = -sum(sum(Y_sink(1:end)));            % x(2) = positive if sink
x(3) = M2 + x(2);
x(1) = M - x(2);

for j = 1:n
    x(j+3) = 0;
    x(j+3) = x(j+3) - sum(Y_sink(j,1:end)); % permeate mass
end

for i = 1:n
    Y2(i) = (M2*Y2_k(i) + x(i+3))/x(3);     % mass per species
end

for j = 1:n
    x(j+3+n) = x(j+3)/sum(x(4:(3+n)));      % permeate mass fraction
end                                         % (only the source/sink species)

end