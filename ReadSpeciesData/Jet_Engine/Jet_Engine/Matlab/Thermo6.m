clear all;close all;clc;
%%
warning off
abspath_to_generalfolder='C:\Users\s137280\Documents\tue leerjaar 2\Thermodynamics\Jet_Engine\Matlab\General'; % absolute reference to General folder
addpath(abspath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Some easy units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;
%% Used by Nasa pols
global Runiv pref
Runiv=8.314472;
pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Nasa is ready
TA=300; % Initial temperature
cFuel='Gasoline';           % Proper indices to database
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});

%% Determine compositions Yi of initial and final mixture and do your stuff
SpS=Sp(iSp);                                                                %% Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];

%% determine Aircomp and Fuelcomp                                                  
Xair = [0 0.21 0 0 0.79]  ;                                      %% mole fraction
Mair = Xair*Mi';
Yair = Xair.*Mi/Mair;
Yfuel = [1 0 0 0 0];

%% 
x = Sp(iSp(1)).Elcomp(3);                           %% line 33
y = Sp(iSp(1)).Elcomp(2);                           %% line 34
a = (x+(y/4));                                      %% line 35

MO2 = Sp(iSp(2)).Mass;
MN2= Sp(iSp(5)).Mass;
MF= Sp(iSp(1)).Mass;
MCO2= Sp(iSp(3)).Mass;
MH2O= Sp(iSp(4)).Mass;

MA=a*(MO2+(Xair(5)/Xair(2))*MN2);                   %% line 43
phi = [0.125:0.025:1];
%%AFstoi = (11.035*Mair/0.21)/Mi(1);
AF = (MA/MF)./phi;                                  %% line 46

%% determine compositions mix AFstoi at beginning
YF = 1./(1+AF);                                     %% line 49
YA = 1-YF;                                          %% line 50
YO_before = 1/(1+(Xair(5)*MN2)/(Xair(2)*MO2)).*YA;    %% line 51    
YN2 = 1/(1+(Xair(2)*MO2)/(Xair(5)*MN2)).*YA;        %% line 52

%% determine compostions mix AFstoi at end
YCO2 = (YF./MF).*x.*MCO2;                           %% line 55
YH2O = (YF./MF).*y/2.*MH2O;                         %% line 56
YO2after = 1-YCO2-YH2O-YN2;                         %% line 57

h1 = YF.*HNasa(TA,Sp(55))+YO_before.*HNasa(TA,Sp(4))+YN2.*HNasa(TA,Sp(48)); %% line 59

%% h2 for every temperature and mixture
i=1;
for T=200:1:3000
    for j = 1:length(AF)
        h2(i,j) = YCO2(j)*HNasa(T,Sp(16))+YH2O(j)*HNasa(T,Sp(6))+YN2(j)*HNasa(T,Sp(48))+YO2after(j)*HNasa(T,Sp(4)); %% line 65
    end
    i=i+1;
end

T_i = [200:1:3000]; 
for i=1:length(AF)                                                                % Compute properties for all species for TR temperature range
    T_ad(i) = interp1(h2(:,i),T_i,h1(i));           %%line 72
end

figure(1);
plot(phi,T_ad);
grid on;
xlabel('Equivalence ratio [-]');
ylabel('Temperature [K]');
title('Adiabatic temperature v.s. Equivalence ratio');

figure(2);
plot(AF,T_ad);
grid on;
xlabel('Air Fuel ratio [-]');
ylabel('Temperature [K]');
title('Adiabatic flame temperature v.s. Air Fuel ratio');

fprintf('%6s %9.2f %9.2f\n','Adiabtic flame temperature at AF =',T_ad(36))
fprintf('%6s %9.2f %9.2f\n','Adiabtic flame temperature at 8*AF =',T_ad(1))
