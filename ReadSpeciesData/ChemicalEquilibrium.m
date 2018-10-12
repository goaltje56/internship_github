%% Read species data
% Make these variable global
global Runiv El Sp;

% MechanismFile = 'H2-CO.trot';
MechanismFile = 'fuels.trot';

% Read the reaction mechanism
[El, Sp] = ReadTrotDat(MechanismFile);
% Number of species
Nsp = length(Sp);
% Universal gas constant
Runiv = 8.314462175;

%% Plot some info
fprintf('\n##  %-16s %12s %12s\n','Species name','H0','cp @1200K');
fprintf('%-20s %12s %12s\n','','(J/mol)','(J/mol.K)');
for i = 1:Nsp
    fprintf('%.2i  %-16s %12.4e %12.4e\n', i, Sp(i).Name, ...
            ThermProp(298.15,Sp(i),'Enthalpy','Mole'), ...
            ThermProp(1200.0,Sp(i),'Cp','Mole'));
end
fprintf('\n');

%% Loop over equivalence ratio
phi = 0.5:0.02:1.5;

% Allocate output arrays
Tequil = zeros(length(phi),1);
Xequil = zeros(length(phi),Nsp);

for i = 1:length(phi)

% Reactants
% Composition [mole fractions]
Xreac = zeros(1,Nsp);
Xreac(strcmp({Sp.Name},'O2')) = 0.21;
Xreac(strcmp({Sp.Name},'N2')) = 0.79;

% Stoichiometric coefficient for hydrogen !! USE H2-CO.trot MECHANISM !!
% nu_mole = 0.5;
% Xreac(strcmp({Sp.Name},'H2')) = 0.21 / nu_mole * phi(i);

% Stoichiometric coefficient for methane
nu_mole = 2.0;
Xreac(strcmp({Sp.Name},'CH4')) = 0.21 / nu_mole * phi(i);

% Stoichiometric coefficient for propane
% nu_mole = 5.0;
% Xreac(strcmp({Sp.Name},'C3H8')) = 0.21 / nu_mole * phi(i);

% Stoichiometric coefficient for n-heptane
% nu_mole = 11.0;
% Xreac(strcmp({Sp.Name},'N-C7H16')) = 0.21 / nu_mole * phi(i);

% Normalize
Xreac = Xreac / sum(Xreac);

% Temperature [K]
Treac = 298.15;
% Pressure [Pa]
p = 1.01325e5;

% Call the equilibrium solver
if i==1
    % First time without initial guess
    [Xeq, Teq] = ChemEquil(Xreac,Treac,p);
    fprintf('\nComputing .');
else
    % Use previous result as initial guess
    [Xeq, Teq] = ChemEquil(Xreac,Treac,p,Xequil(i-1,:),Tequil(i-1));
    fprintf('.');
end

% Save results
Tequil(i) = Teq;
Xequil(i,:) = Xeq;

end % Loop over phi
fprintf('\n');

%% Plotting
figure(1)

subplot(1,2,1)
plot(phi,Tequil,'-')
xlabel('Equivalence ratio')
ylabel('Temperature [K]')

subplot(1,2,2)
iO2  = find(strcmp({Sp.Name},'O2'));
iCO2 = find(strcmp({Sp.Name},'CO2'));
iCO  = find(strcmp({Sp.Name},'CO'));
iH2O = find(strcmp({Sp.Name},'H2O'));
iH2  = find(strcmp({Sp.Name},'H2'));
plot(phi,Xequil(:,[iO2 iCO2 iCO iH2O iH2]),'-')
xlabel('Equivalence ratio')
ylabel('Species mole fractions')
legend(Sp([iO2 iCO2 iCO iH2O iH2]).Name)

%%
figure(2)
iNO  = find(strcmp({Sp.Name},'NO'));
plot(phi,Xequil(:,[iNO])*1e6,'-')
xlabel('Equivalence ratio')
ylabel('NO (ppm)')

%%
% set(gcf,'PaperPosition',[0 0 25 10])
% print -depsc equil.eps
% print -dpng -r300 equil.png

