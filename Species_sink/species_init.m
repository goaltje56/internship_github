function [mass, moles, Y_k, X_k, iAll, MW, D, sink, n] = species_init(NPI)
MechanismFile = 'fuels.trot';

global Sp

%% input moles
moles = [100; 50; 50; 50];
sink  = [1 0 1 0];

% Define some species
iO2  = find(strcmp({Sp.Name},'O2'));
iCO2 = find(strcmp({Sp.Name},'CO2'));
iN2  = find(strcmp({Sp.Name},'N2'));
iAr  = find(strcmp({Sp.Name},'AR'));

iAll = [iO2 iCO2 iN2 iAr];
n           = length(iAll);                % number of species [-]


MW   = [Sp(iAll).Mass];%[31.9988 44.01 28.0134 39.9480];      % Molar weight of species [gr/mol]

X_k = moles/(sum(moles));
m_total = MW*moles;

    for i = 1:n
        Y_k(i,:) = X_k(i)*MW(i)/(MW*X_k);
    end
mass = Y_k*m_total;
D       = species_diff(NPI, 300, iAll, iAll, 'Diffusivity', n); % Diffusivity of species [m^2/s]
    
end
