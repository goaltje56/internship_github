function Y_k = species_bound(NPI, n, Y_k)

Y_k = Y_k./sum(Y_k);
% for j = 1:n    
%     for I = 2:NPI+1
%         if j ==1
%             Y_total(1,I) = 0;
%         end
%         Y_total(1,I) = Y_total(1,I) + Y_k(j,I);        
%     end
% end
% 
% for j = 1:n   
%     for I = 2:NPI+1      
%         Y_k(j,I) = Y_k(j,I)/Y_total(1,I);
%     end    
% end

end