%% To correct the pressure and the volocity see eq 6.24, 6.25
  
function b_p = pcorr(NPI, rho,A, u, aP_u, aW_u, aE_u, x, x_u)
    for I=3:NPI+1
        i = I;
        b_p(i) = (rho(I-1)*(x(I)-x_u(i)))+(rho(I)*(x_u(i)-x(I-1)))/(x(I)-x(I-1)); 
        b_p(i) = (b_p(i)*A*u(i)*(aW_u(i)+aE_u(i)))/aP_u(i);
    end
end