function [x theta] = stagecut(M, Y_sink)
% Mr        = x(1)      TOTAL mass at retenate
% Msink     = x(2)      TOTAL mass that permeates

%% these equations depend on the boundary conditions!!

x(2) = -sum(sum(Y_sink(1:end)));            % x(2) = positive if sink
x(1) = M - x(2);
theta = x(2)/M;

end