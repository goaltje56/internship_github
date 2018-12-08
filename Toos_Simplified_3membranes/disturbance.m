function Y_in = disturbance(dEND, dSTART, Dt, Y_in, Y_final)
x = [0 : Dt : dEND-dSTART];

Y_in =  Y_final+(Y_in-Y_final)*exp(-0.5*x);
end 
