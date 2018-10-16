clear all;close all;clc;
%% 
warning off
abspath_to_generalfolder='C:\Users\s137280\Documents\tue leerjaar 2\Thermodynamics\Jet_Engine\Matlab\General'; % absolute reference to General folder
addpath(abspath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Used by Nasa pols should not be changed
global Runiv Pref
Runiv=8.314472;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Given conditions. Take the ones from your specific case                  
v1=200;Tamb=250;P3overP2=7;Pamb=55*kPa;mfurate=0.68*kg/s;AF=71.25;             % These are the ones from the book
cFuel='Gasoline';
%% Select all species
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air comp
Xair = [0 0.21 0 0 0.79];                                                   % Order is important
MAir = Xair*Mi';
Yair = Xair.*Mi/MAir;
%% Fuel comp
Yfuel = [1 0 0 0 0];                                                        % Only fuel
%% Range of enthalpies/thermal part of entropy of species
TR = [200:1:3000];
for i=1:NSp                                                                 % Compute properties for all species for TR temperature range
    hia(:,i) = HNasa(TR,SpS(i));
    sia(:,i) = SNasa(TR,SpS(i));
end
hair_a= Yair*hia';                                                          % enthalpy of air for range of T
sair_a= Yair*sia';                                                          % thermal part of entropy of air for range of T
%% [1-2] Diffusor
sPart = 'Diffusor';
T1 = Tamb;
P1 = Pamb;
Rg = Runiv/MAir;
for i=1:NSp
    hi(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi';
v2=0;
h2 = h1+0.5*v1^2-0.5*v2^2;
T2 = interp1(hair_a,TR,h2);                                                 % Not exactly correct but nearly. Why??? can also do a search
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';                                                        % Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;
Pr = exp(lnPr);
P2 = P1*Pr;
S1  = s1thermal - Rg*log(P1/Pref);
S2  = s2thermal - Rg*log(P2/Pref);
% Print to screen
fprintf('\nStage %12s\n      | %9i %9i\n',sPart,1,2);
fprintf('%6s| %9.2f %9.2f\n','Temp',T1,T2);
fprintf('%6s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
%% Here starts your part

%% determining fractions of the mix at the beginning.
mair1= AF*Yair;
Nair1 = mair1./Mi;
Nfuel = Yfuel./Mi;
Xair1 = Nair1./sum(Nair1);
mtotal= mair1+Yfuel;
Ymix1 = mtotal/(AF+1);
Nmix1 = mtotal./Mi;                          %%amount of mole per specie at the beginning
Xmix1 = Nmix1./sum(Nmix1);

%%determining fractions of the mix at the end
NCO2_H2O=Nmix1(1).*[0 0 7.76 6.55 0];
MCO2_H2O = NCO2_H2O.*Mi;
NO_used = Nmix1(1).*[0 11.035 0 0 0];  %% the amount of mole O2 necessary for 1 kg combustion
Nair_after = Nair1-NO_used;
N_after = Nair_after+NCO2_H2O;
Xafter=N_after./sum(N_after);
Mafter = Xafter*Mi';
Yafter = Xafter.*Mi./sum(Mafter);

%%[2-3 compressor]
P3 = P2*P3overP2;
cp_airT2= CpNasa(T2,SpS(2))*Yair(2)+CpNasa(T2,SpS(5))*Yair(5);
cv_airT2= CvNasa(T2,SpS(2))*Yair(2)+CvNasa(T2,SpS(5))*Yair(5);
gamma_air=cp_airT2/cv_airT2;
T3 = T2*(P3/P2)^((gamma_air-1)/gamma_air);
for i=1:NSp
    hl(i)    = HNasa(T3,SpS(i));
end
h3 = Yair*hl';
%% Wincomp= mflow *(h3-h2)
mairflow=mfurate*AF;
Wincomp = mairflow*(h3-h2);
v3 = 0;

%%[3-4 combustor]

P4=P3; %% the pressure is constant
%% T4 is adiabaticflame temperature
h4 = Yafter*hia';
T4 = interp1(h4,TR,h3);
for i=1:NSp
    hl2(i)    = HNasa(T4,SpS(i));
end

h4check = Yafter*hl2';

v4 = 0; %% behoud van massa

%% [4-5 Turbine]

h5 = h2+h4check-h3;          %%conservation of energy
v5 = 0;

h_5 = Yafter*hia';
T5 = interp1(h_5,TR,h5);
for i=1:NSp
    hl3(i)    = HNasa(T5,SpS(i));
end



