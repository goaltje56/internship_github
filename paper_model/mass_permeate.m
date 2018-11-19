function [Yp Mp] = mass_permeate(Y_in, Y_k, M, rho_real, u_in)
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
if M-rho_real*u_in ~= 0
    Yp = (M*Y_in - rho_real*u_in*Y_k)/ (M-rho_real*u_in);
    Mp = (M-rho_real*u_in);
else 

    Yp = 0;
    Mp = 0;
end
end