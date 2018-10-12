function Y_k = species_bound(NPI, n, Y_k)

for j = 1:n    
    for I = 1:NPI+2
        if j ==1
            Y_total(1,I) = 0;
        end
        Y_total(1,I) = Y_total(1,I) + Y_k(j,I);        
    end
end

for j = 1:n   
    for I = 1:NPI+2      
        Y_k(j,I) = Y_k(j,I)/Y_total(1,I);
    end    
end

end