function Y_k = MoleToMass(X_in, Y_k)
    
global Sp

% Define some species
iO2  = find(strcmp({Sp.Name},'O2'));
iCO2 = find(strcmp({Sp.Name},'CO2'));
iN2  = find(strcmp({Sp.Name},'N2'));
iAr  = find(strcmp({Sp.Name},'AR'));

iAll = [iO2 iCO2 iN2 iAr];
n           = length(iAll);                % number of species [-]


MW   = [Sp(iAll).Mass];      % Molar weight of species [gr/mol]

X_k = X_in/(sum(X_in));

    for i = 1:n
        Y_k(i,1) = X_k(i)*MW(i)/(MW*X_k);
    end

end