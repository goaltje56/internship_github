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
v1=200;Tamb=250;P3overP2=7;Pamb=55*kPa;mfurate=0.68*kg/s;AF=71.25;             % These are the ones from group 138
cFuel='Gasoline';
%% Select all species
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air comp
Xair = [0 0.21 0 0 0.79];                        %%line 25                           % Order is important
MAir = Xair*Mi';
Yair = Xair.*Mi/MAir;                            %%line 27
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
v2=0;                           %% line 47
h2 = h1+0.5*v1^2-0.5*v2^2;      %% line 48
T2 = interp1(hair_a,TR,h2);     %% line 49                                           % Not exactly correct but nearly. Why??? can also do a search
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
    si1(i)    = SNasa(T1,SpS(i));
    si2(i)    = SNasa(T2,SpS(i));
end
h2check = Yair*hi2';             %% line 55                                % Why do I do compute this h2check value? Any ideas?
s1thermal = Yair*si1';
s2thermal = Yair*si2';
lnPr = (s2thermal-s1thermal)/Rg;        
Pr = exp(lnPr);
P2 = P1*Pr;                         %% line 60
S1  = s1thermal - Rg*log(P1/Pref);
S2  = s2thermal - Rg*log(P2/Pref);

% Print to screen
fprintf('\nStage %12s\n      | %9i %9i\n',sPart,1,2);
fprintf('%6s| %9.2f %9.2f\n','Temp',T1,T2);
fprintf('%6s| %9.2f %9.2f  [kPa]\n','Press',P1/kPa,P2/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v1,v2);
%% Here starts your part
%% determining fractions of the mix at the beginning.

x = Sp(iSp(1)).Elcomp(3);                   %% line 72        
y = Sp(iSp(1)).Elcomp(2);                   %% line 73        
a = (x+(y/4));                              %% line 74        

MO2 = Sp(iSp(2)).Mass;                      %% line 76
MN2= Sp(iSp(5)).Mass;
MF= Sp(iSp(1)).Mass;
MCO2= Sp(iSp(3)).Mass;
MH2O= Sp(iSp(4)).Mass;                      %% line 80

YF = 1./(1+AF);                             %% line 82        
YA = 1-YF;                                  %% line 83
YO_before = 1/(1+Xair(5)*MN2/(Xair(2)*MO2)).*YA;  %% line 84      
YN2 = 1/(1+(Xair(2)*MO2)/(Xair(5)*MN2)).*YA;      %% line 85  

Ymix1 = [YF YO_before 0 0 YN2];             %% line 87
Nmix1 = Ymix1./Mi;
Xmix1 = Nmix1./sum(Nmix1);
Mmix1 = Xmix1*Mi';

Xair = [0 0.21 0 0 0.79];                                                   % Order is important
MAir = Xair*Mi';
Yair = Xair.*Mi/MAir;

%%determining fractions of the mix at the end
YCO2 = (YF./MF).*x.*MCO2;                           %% line 97
YH2O = (YF./MF).*(y/2).*MH2O;                       %% line 98
YO2after = 1-YCO2-YH2O-YN2;                         %% line 99

Yafter = [0 YO2after YCO2 YH2O YN2];                %% line 100
Nafter = Yafter./Mi;
Xafter = Nafter./sum(Nafter);
Mafter = Xafter*Mi';

%%[2-3 compressor]
sPart='Compressor';
P3 = P2*P3overP2;                                   %% line 108
s_druk3 = Rg*log(P3/Pref);           %%s3 can be found since S2=S3 and P3 is known
s3thermal = S2+s_druk3;

T3 = interp1(sair_a,TR,s3thermal);  %% line 112 because we now s3thermal we can find T3.
v3=0;                                %%conservation of energy line 113

for i=1:NSp
    hl(i)    = HNasa(T3,SpS(i));
end
h3 = Yair*hl';                      %%line 118 finding h3


mairflow=mfurate*AF;                %% line 121
Wincomp = mairflow*(h3-h2);         %% line 122

fprintf('\nStage %12s\n      %9i\n',sPart,3);
fprintf('%6s| %9.2f\n','Temp',T3);
fprintf('%6s| %9.2f  [kPa]\n','Press',P3/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v3);


%%[3-4 combustor]
sPart='Combustor';
%% h3 = h4                          line 132
P4=P3; %% line 133 the pressure is constant

%% T4 is adiabaticflame temperature
h_4= Ymix1*hl';             %% line 136
h4 = Yafter*hia';           %% line 137
T4 = interp1(h4,TR,h_4);    %% line 138
for i=1:NSp
    hl2(i)    = HNasa(T4,SpS(i)); %% line 140
end
h4check = Yafter*hl2';

for i=1:NSp
si4(i) = SNasa(T4,SpS(i));
end
s4thermal = Yafter*si4';
s4thermal = sum(s4thermal);

%% Rg2 = Runiv/Mafter

S4 = s4thermal - Rg*log(P4/Pref) ;      %%line 152

v4 = 0; %% line 154 energy conversation
fprintf('\nStage %12s\n      %9i\n',sPart,4);
fprintf('%6s| %9.2f\n','Temp',T4);
fprintf('%6s| %9.2f  [kPa]\n','Press',P4/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v4);

%% [4-5 Turbine]
sPart = 'Turbine';

h5 = h4check-(mairflow*(h3-h2)/(mairflow+mfurate));    %% line 163 conservation of energy
v5 = 0;                                                %% line 164 

h_5 = Yafter*hia';                                     %% line 166
T5 = interp1(h_5,TR,h5);                               %% line 167
for i=1:NSp
    hi5(i)    = HNasa(T5,SpS(i));
end
h5check = Yafter*hi5';

for i=1:NSp
    si5(i)    = SNasa(T5,SpS(i));
end
s5thermal = Yafter*si5';

lnPr2 = (s5thermal-s4thermal)/Rg;
Pr2 = exp(lnPr2);
P5 = P4*Pr2;                    %%line 180
S4  = s4thermal - Rg*log(P4/Pref);
S5  = s5thermal - Rg*log(P5/Pref);
fprintf('\nStage %12s\n      %9i\n',sPart,5);
fprintf('%6s| %9.2f\n','Temp',T5);
fprintf('%6s| %9.2f  [kPa]\n','Press',P5/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v5);

%%[5-6 Nozzle]
sPart = 'Nozzle';

P6 = Pamb;                      %%line 191
s_druk6 = Rg*log(P6/Pref);
s6thermal = S5+s_druk6;
s_6= Yafter*sia';

T6 = interp1(s_6,TR,s6thermal); %%line 196
for i=1:NSp
    hi6(i)    = HNasa(T6,SpS(i));
end
h6 = Yafter*hi6';               %% line 200
v6 = (2*(h5-h6))^0.5;           %% line 201

fprintf('\nStage %12s\n      %9i\n',sPart,6);
fprintf('%6s| %9.2f\n','Temp',T6);
fprintf('%6s| %9.2f  [kPa]\n','Press',P6/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v6);