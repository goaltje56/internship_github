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
%%TA=300; % Initial temperature
cFuel='Gasoline';           % Proper indices to database
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});

%% Determine compositions Yi of initial and final mixture and do your stuff
SpS=Sp(iSp);                                                                %% Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];

AFstoi = (11.035/0.21)*28.84/(12*7.76+1.008*13.1);                           %% mass of air for stoi devided by mass of fuel

%% determine Aircomp and Fuelcomp                                                  
Xair = [0 0.21 0 0 0.79]  ;                                      %% mole fraction
Mair = Xair*Mi';
Yair = Xair.*Mi/Mair;
Yfuel = [1 0 0 0 0];

%% determine total compositions at beginning
Yair1 = AFstoi.*Yair;
Y1 = (Yfuel+Yair1)./(1+sum(Yair1));

for i=1:NSp
    hi(i)    = HNasa(Tref,SpS(i));
end

%% total enthalpy before
hfuel = Yfuel*hi';
hair = Yair*hi';
%%h1 =  Yfuel*hi' + Yair*hi';
h1 = Y1*hi';

%% determine comp after reaction
N_total = [0 0 7.76 6.55 ((11.035*0.79)/0.21)]   ;          %%total amount of mole of every species
Ntotal = sum(N_total);
X2 = N_total./sum(N_total);
M2 = X2*Mi';                                                
Y2 = X2.*Mi/M2;

TR = [200:1:3000]; 
for i=1:NSp                                                                 % Compute properties for all species for TR temperature range
    hia(:,i) = HNasa(TR,SpS(i));
end
                    
h2 = h1;                        %% enthalpy at the end is the same
h2_a = Y2*hia';

T2 = interp1(h2_a ,TR ,h1)
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
end
h2check = Y2*hi2';

%% the amount of air for AF 0.2-1
j = [0.2:0.001:1];
for i=1:NSp
Yair1i(:,i) = j.*AFstoi.*Yair(i);
end 



