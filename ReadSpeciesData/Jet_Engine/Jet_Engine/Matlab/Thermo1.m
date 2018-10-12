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

%% determine comp before reaction                                                   
X1 = [1 11.035 0 0 ((11.035/0.21)*0.79)]  ;                                      %% mole fraction
M1 = X1.*Mi;
Y1 = M1./Mi;

%% determine comp after reaction
X2 = [0 0 7.76 6.55 ((11.035/0.21)*0.79)]   ;                                                                
M2 = X2.*Mi;
Y2 = M2./Mi;

%% get rid of the inf

for i=1:NSp
    if Y1(:,i) == inf
        Y1(:,i)=0;
    end
    if Y2(:,i) == inf
        Y2(:,i)=0;
    end
end

TR = [200:1:5000]; 
for i=1:NSp                                                                 % Compute properties for all species for TR temperature range
    hia(:,i) = HNasa(TR,SpS(i));
end

for i=1:NSp
    hi(i)    = HNasa(Tref,SpS(i));
end

h1 = Y1*hi';                    %% the enthalpy at the beginning
h2 = h1;                        %% enthalpy at the end is the same
h2_a = Y2*hia';

T2 = interp1(h2_a ,TR ,h1)
% T=[0:50:2000]
% for i =1:NSp
%     Hdummy(T) = HNasa(T(i),SpS(i))
% end
%i = 0;
%for T = 50:50:3000
 %   i = i + 1;
  %  Tdummy(i) = T;
   % ho_dummy(i) = HNasa(T,Sp(4));
%end
%interp1(ho_dummy,Tdummy,2*10E-6)

