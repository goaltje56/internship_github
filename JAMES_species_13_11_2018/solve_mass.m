function [F] =  solve_mass(x, M, XA, XB, XC, XD, XAr, XBr, XCr, XDr, sink)
% Mr    = x(1)
% Mp    = x(2)
% XAp   = x(3)
% XBp   = x(4)
% XCp   = x(5)
% XDp   = x(6)

%% these equations depend on the boundary conditions!!

F(1) = M*XA - (x(1)*XAr + x(2)*x(3)*sink(1));
F(2) = M*XB - (x(1)*XBr + x(2)*x(4)*sink(2));
F(3) = M*XC - (x(1)*XCr + x(2)*x(5)*sink(3));
F(4) = M*XD - (x(1)*XDr + x(2)*x(6)*sink(4));

F(5) = M - ( x(1) + x(2) );
F(6) = 1 - (x(3)*sink(1) + x(4)*sink(2) + x(5)*sink(3) + x(6)*sink(4));

end