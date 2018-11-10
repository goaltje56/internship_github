function [mass, moles, rho_s, Y_k, X_k, Y_in, X_in, iAll, MW, D, D_k, P_k, f_old, sink, n] = species_init(NPI, massflow)

%% Read species data
% Make these variable global
global El Sp Patm

MechanismFile = 'fuels.trot';

% Read the reaction mechanism
[El, Sp] = ReadTrotDat(MechanismFile);

%% Define some species
iO2  = find(strcmp({Sp.Name},'O2'));
iCO2 = find(strcmp({Sp.Name},'CO2'));
iN2  = find(strcmp({Sp.Name},'N2'));
iAr  = find(strcmp({Sp.Name},'AR'));

iAll  = [iO2 iCO2 iN2 iAr];

MW   = [Sp(iAll).Mass];                  % Molar weight of species [gr/mol]
n = length(iAll);

%% mole flow rate in system and sink term
if massflow == 1
    mass = [4000;   2000;   1000;    1000 ];
    moles = (mass./MW');
else
    moles = [100;   50;   50;    50 ];
%     moles = [20.78; 79.22];
end

    X_in  = moles/sum(moles);
 
    
sink  = [1   0   1  0];
% P_n   = [7.155; 1.255]*10^(-9);
P_n   = [7.155; 3.16; 1.255; 0  ]*10^(-9);  % Permeability of species
rho_s = [1.429 1.98 1.2504 1.784];          % 'Real' density of species [kg/m^3]
% rho_s = [1.429 1.2504];          % 'Real' density of species [kg/m^3]



    for i = 1:n
        Y(i,:) = X_in(i)*MW(i)/(MW*X_in);
    end
   Y_in  = Y(:,1);
m_total = MW*moles;                         % total mass    
mass = Y*m_total;                           % total mass per species.
D       = species_diff(NPI, 300, iAll, iAll, 'Diffusivity', n); % Diffusivity of species [m^2/s]
  
    for j = 1:n
        
        for I = 1:NPI+2
            p_k (j,I)   = Patm;             % pressure of species k
            Y_k(j,I)    = Y(j);             % mass of species k
            MW_k(j,I)   = MW(j);            % molar weight of species k
            P_k(j,I)    = P_n(j);           % Permeability
            X_k(j,I)    = X_in(j);
        end

    end

    for j = 1:n
        for i = 1:NPI+2
            D_k(j,i) = D(j,:)*Y_k(:,i);    % Diffusion is in terms of mass 
        end                                % so use massfraction!!
    end
    
    f_old = Y_k;
end
