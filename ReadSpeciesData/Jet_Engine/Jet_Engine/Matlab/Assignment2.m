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
sair_a= Yair*sia';                                                          % thermal part os entropy of air for range of T
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
%%[2-3 compressor]
P3 = P2*P3overP2;
gamma_air=1.4;
T3 = T2*(P3/P2)^((gamma_air-1)/gamma_air);
%% Wc = mflow(h3-h2) and Wc/mflow = Cp(T3-T2)= 1/2 v3^2
cp3


