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
%% [2-3] compressor
sPart='Compressor';
P3=P2*P3overP2;
gamma= 1.4; %%dit heb ik uit het boek gehaald en weet niet of het altijd geldt
T3=T2*P3overP2^((gamma-1)/gamma);
v3=0;
fprintf('\nStage %12s\n      %9i\n',sPart,3);
fprintf('%6s| %9.2f\n','Temp',T3);
fprintf('%6s| %9.2f  [kPa]\n','Press',P3/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v3);

%%[3-4] combustor
sPart='Combustor';
P4=P3;
for i=1:5
     hi3(i)=HNasa(T3,SpS(i));                   
end

for i=1:length(TR)                                
        for k=1:5
            hi4(k,i)=HNasa(TR(i),SpS(k));
        end      
    end

Yi3=[0.0096 0.2298 0 0 0.7606];    
Yi4=[0 0.1912 0.0265 0.0217 0.7606 ];

h3=sum(hi3.*Yi3);
h4=hi4'*Yi4';

T4=interp1(h4, TR, h3);

v4=0;

fprintf('\nStage %12s\n      %9i\n',sPart,4);
fprintf('%6s| %9.2f\n','Temp',T4);
fprintf('%6s| %9.2f  [kPa]\n','Press',P4/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v4);

%%Turbine [4-5]
sPart='Turbine';
h5 = h4 + h2 - h3;
s5 = s4;
v5=0;

for i=1:NSp
    hi5(i) = HNasa(T5,SpS(i));
    si4(i) = SNasa(T4,SpS(i));
    si5(i) = SNasa(T5,SpS(i));
end
h5check=Yi4*hi5';
s4thermal=Yi4*si4';
s5thermal=Yi4*si5';
lnPr=(s5thermal-s4thermal)/Rg;
Pr=exp(lnPr);
P5=P4*Pr;
smix_a=Yi4*sia';

fprintf('\nStage %12s\n      %9i\n',sPart,5);
fprintf('%6s| %9.2f\n','Temp',T5);
fprintf('%6s| %9.2f  [kPa]\n','Press',P5/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v5);

%%Nozzle[5-6]
Spart='Nozzle';
sPart = 'Nozzle';
P6 = Pamb;

for i=1:NSp
    si6(i)=SNasa(T5,SpS(i));
end
S6=sair_a-Rg*log(P6/Pref);
T6=interp1(S6,TR,s5);

for i=1:NSp
    hi6(i)  =HNasa(T6,SpS(i));
end

h6=0;
for i=1:5
    h6=h6+Yi6(i)*hi6(i);
end
v6 = sqrt(2*(h5check-h6));

fprintf('\nStage %12s\n      %9i\n',sPart,6);
fprintf('%6s| %9.2f\n','Temp',T6);
fprintf('%6s| %9.2f  [kPa]\n','Press',P6/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v6);




