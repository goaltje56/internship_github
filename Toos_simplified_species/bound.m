function Y_k  = bound(NPI,Y_in, Y_k)
        
    Y_k(:,1) = Y_in;

    Y_k(1, NPI+2) = Y_k(1, NPI+1);
    Y_k(2, NPI+2) = Y_k(2, NPI+1);
    Y_k(3, NPI+2) = Y_k(3, NPI+1);
    Y_k(4, NPI+2) = Y_k(4, NPI+1);

end