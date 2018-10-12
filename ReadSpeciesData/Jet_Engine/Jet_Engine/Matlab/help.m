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
v1=200;Tamb=250;P3overP2=9;Pamb=55*kPa;mfurate=0.68*kg/s;AF=102.78;             % groupsdata
cFuel='CH4';
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
    hi1(i)    = HNasa(T1,SpS(i));
end
h1 = Yair*hi1';
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

%% [2-3] compressor
sPart='Compressor';
P3=P2*P3overP2;
v3=0;
s_druk = Rg*log(P3/Pref);           

s3thermal=S2+s_druk;
T3=interp1(sair_a,TR,s3thermal);

for i=1:NSp
    h_3(i)    = HNasa(T3,SpS(i));
end
h3 = Yair*h_3';

%%[3-4] combustor
sPart='Combustor';
P4=P3;
v4=0;

MCH4=SpS(1).Mass;                  
Mair=0.21*SpS(2).Mass+0.79*SpS(5).Mass;

Yi3=[0.0096 0.2298 0 0 0.7606];    
Yi4=[0 0.1912 0.0265 0.0217 0.7606 ];    

h4a = Yi3*h_3';
h4 = Yi4*hia';
T4=interp1(h4,TR,h4a)


for i=1:NSp
    h_4(i)    = HNasa(T4,SpS(i));
end
h4b = Yi4*h_4';    


for i=1:NSp
si4(i) = SNasa(T4,SpS(i));
end
s4thermal = Yi4*si4';
s4thermal = sum(s4thermal);
S4 = s4thermal - Rg*log(P4/Pref) ;

 %% [4-5] Turbine
sPart = 'Turbine';
v5=v4
h5 = h4b + h2 - h4a;
h_5 = Yi4*hia';

T5=interp1(h_5,TR,h5)

for i=1:NSp
    hi5(i)  =HNasa(T5,SpS(i));
    si4(i)  =SNasa(T4,SpS(i));
    si5(i)  =SNasa(T5,SpS(i));
end
h5check=Yi4*hi5';
s5thermal=Yi4*si5';

lnPr_a=(s5thermal-s4thermal)/Rg;
Pr_a = exp(lnPr_a);
P5=P4*Pr_a

S5  = s5thermal - Rg*log(P5/Pref);

%% [5-6] Nozzle
sPart = 'Nozzle';
P6 = Pamb
S6 = S5;
s6thermal= S5+ (Rg*log(P6/Pref));
s_6= Yi4*sia';

T6 = interp1(s_6 , TR, s6thermal);

for i=1:NSp
    hi6(i)    = HNasa(T6,SpS(i));
end

h6=Yi4*hi6';
v6 = sqrt(2*(h5check-h6))