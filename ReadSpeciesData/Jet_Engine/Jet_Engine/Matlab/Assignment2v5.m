clear all;close all;clc;
%% 
warning off
abspath_to_generalfolder='C:\Users\s149640\Documents\Thermodynamica\Jet engine\Matlab\General'; % absolute reference to General folder
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
v1=200;Tamb=300;P3overP2=9;Pamb=100*kPa;mfurate=0.58*kg/s;AF=204.42;             % These are the ones from the book
cFuel='H2';
%% Select all species
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
NSp = length(SpS);
Mi = [SpS.Mass];
%% Air comp
Xair = [0 0.21 0 0 0.79]                                                    % Order is important
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
%% Diffusor [1-2]
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
fprintf('\n');
%% Compressor [2-3]
sPart='Compressor';

v3=0;
P3=P2*P3overP2;
S3=S2;

s3thermal=S3+Rg*log(P3/Pref);
T3=interp1(sair_a, TR, s3thermal);

mairrate=AF*mfurate;

for j=1:Nsp
    hiair3(i)=HNasa(T3,SpS(i));
end

hair3=Yair*hiair3';
% W1=mairrate*(h2-hair3);
% Wdot blz 523

% Print to screen
fprintf('\nStage %12s\n      %9i\n',sPart,3);
fprintf('%6s| %9.2f\n','Temp',T3);
fprintf('%6s| %9.2f  [kPa]\n','Press',P3/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v3);
fprintf('\n');
%% Combuster [3-4]
sPart='Combustor';

NH2=1;

Nair1=((NH2*Mi(1))*AF)/(0.21*Mi(2)+0.79*Mi(5));
NO21=Nair1;
NN21=3.76*NO21;

NH2O2=NH2;
NO22=Nair1-0.5*NH2;
NN22=3.76*NO21;

Xmix1=[NH2 NO21 0 0 NN21];
Xmix2=[0 NO22 0 NH2O2 NN22];
Mmix1=Xmix1.*Mi;
Mmix2=Xmix2.*Mi;


for j=1:NSp
    Ymix1(j)=Mmix1(j)/(sum(Mmix1));
    Ymix2(j)=Mmix2(j)/(sum(Mmix2));
end

h_4 = Ymix1*hiair3';
h4 = Ymix2*hia';
T4 = interp1(h4,TR,h_4);

P4=P3;
v4=0;

for i=1:NSp
    hi2(i)    = HNasa(T4,SpS(i));
end
h4check = Ymix2*hi2';

for i=1:NSp
si4(i) = SNasa(T4,SpS(i));
end
s4thermal = Ymix2*si4';
s4thermal = sum(s4thermal);

S4 = s4thermal - Rg*log(P4/Pref) ;

v4 = 0; %% energy conversation

% Print to screen
fprintf('\nStage %12s\n      %9i\n',sPart,4);
fprintf('%6s| %9.2f\n','Temp',T4);
fprintf('%6s| %9.2f  [kPa]\n','Press',P4/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v4);
fprintf('\n');
%% Turbine [4-5]
sPart='Turbine';
S5=S4;
v5=v4;

h5=h4check+((h2-hair3)*mairrate)/(mfurate+mairrate);
%h5=sum(h5);
h_5 = Ymix2*hia';
T5 = interp1(h_5,TR,h5);

for i=1:NSp
    si4(i)=SNasa(T4,SpS(i));
    hi4(i)=HNasa(T4,SpS(i));
    si5(i)=SNasa(T5,SpS(i));
    hi5(i)=HNasa(T5,SpS(i));
end

s4thermal=Ymix2*si4';
s5thermal=Ymix2*si5';
h4=Ymix2*hi4';
h5check=Ymix2*hi5';

lnPr=((s5thermal-s4thermal)/Rg);
Pr=exp(lnPr);
P5=P4*Pr;
%P5=Pref*exp((-S4+s5thermal)/Rg);
S4  = s4thermal - Rg*log(P4/Pref);
S5  = s5thermal - Rg*log(P5/Pref);

% Print to screen
fprintf('\nStage %12s\n      %9i\n',sPart,5);
fprintf('%6s| %9.2f\n','Temp',T5);
fprintf('%6s| %9.2f  [kPa]\n','Press',P5/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v5);
fprintf('\n');
%% Nozzle [5-6]
sPart='Nozzle';
%S6=S5;

S_6=Ymix2*sia';

P6=Pamb;
s6thermal=S5+Rg*log(P6/Pref);
T6=interp1(S_6,TR,s6thermal);
%T6=interp1(smix4,TR,s6thermal);

for j=1:NSp
    hi6(i)=HNasa(T6,SpS(i));
end
h6=Ymix2*hi6';
v6=sqrt(2*(h5check-h6));

% Print to screen
fprintf('\nStage %12s\n      %9i\n',sPart,6);
fprintf('%6s| %9.2f\n','Temp',T6);
fprintf('%6s| %9.2f  [kPa]\n','Press',P6/kPa);
fprintf('%6s| %9.2f %9.2f  [m/s]\n','v',v6);
fprintf('\n');