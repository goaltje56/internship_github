function [F] =  solve_mass2(x, M, XA, XB, XAr, XBr)
% Mr    = x(1)
% Mp    = x(2)
% XAp   = x(3)
% XBp   = x(4)
% XCp   = x(5)
% XDp   = x(6)

%% these equations depend on the boundary conditions!!

F(1) = M*XA - (x(1)*XAr + x(2)*x(3));
F(2) = M*XB - (x(1)*XBr);
% F(3) = M*XC - (x(1)*XCr);
% F(4) = M*XD - (x(1)*XDr);

F(3) = M - ( x(1) + x(2) );
% F(6) = 1 - (x(3) + x(4) + x(5) + x(6));


end