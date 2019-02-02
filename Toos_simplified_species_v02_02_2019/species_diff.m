function D_k = species_diff(NPI, T, iSp1, iSp2, Prop, n)
global Sp
% D_k = zeros(n,n);

for i = 1:n
    for j = 1:n
        D_k(i,j) = DiffProp(T,Sp(iSp1(i)),iSp2(j),'Diffusivity');
    end


end 
    
end
