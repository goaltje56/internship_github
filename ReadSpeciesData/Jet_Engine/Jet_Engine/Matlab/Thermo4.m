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


%% determine Aircomp and Fuelcomp                                                  
Xair = [0 0.21 0 0 0.79]  ;                                      %% mole fraction
Mair = Xair*Mi';
Yair = Xair.*Mi/Mair;
Yfuel = [1 0 0 0 0];

x = Sp(iSp(1)).Elcomp(3);
y = Sp(iSp(1)).Elcomp(2);
a = (x+(y/4))

MO = Sp(iSp(1)).Mass;
MN2= Sp(iSp(4)).Mass;
MF= Sp(iSp(7)).Mass;
MCO2= Sp(iSp(2)).Mass;
MH2O= Sp(iSp(3)).Mass;

MA=a*(MO2+Xair(5)/Xair(2))*MN2;
phi = [0.125:0.025:1]
AFstoi = (11.035*Mair/0.21)/Mi(1);
AF = 8*AFstoi;

%% determine compositions mix AFstoi at beginning
Mairin= AF.*Yair;



Ninitial=[1 11.035 0 0 (11.035*0.79/0.21)];
Nintot = sum(Ninitial);
Xinitial= Ninitial./Nintot;
Minitial=Xinitial*Mi';
Yinitial=Xinitial.*Mi/Minitial;

%% determine compostions mix AFstoi at end
Xend = [0 0 7.76 6.55 (11.035*0.79/0.21)];
Mend = Xend*Mi';
Yend = Xend.*Mi/Mend;

for i=1:NSp
    hi(i)    = HNasa(Tref,SpS(i));
end

%% h for AFstoi
hinitial= Yinitial*hi';

%% determine composition mix 8*AFstoi at beginning
AF2 = 8*AFstoi;
Mair2 = Mi(1)*AF2;
Mtot = Mair2+Mi(1);
Mtotair2 = Mair2.*Yair;
Mtot2 = Mtotair2;
Mtot2(1)=Mi(1);
Yin = Mtot2./Mtot;
Nin = Yin.*Mtot./Mi;
Nintot = sum(Nin);
Xin = Nin./Nintot;

%% determine composition mix 8*AFstoi after
Nend=Nin-Ninitial;
Nend(3)= Xend(3);
Nend(4)= Xend(4);
Nendtot=sum(Nend);
Xend2= Nend./Nendtot;
Mend2= Xend2*Mi';
Yend2 = Xend2.*Mi/Mend2;

%% h for 8*AFstoi
hinitial2=Yin*hi';

TR = [200:1:3000]; 
for i=1:NSp                                                                 % Compute properties for all species for TR temperature range
    hia(:,i) = HNasa(TR,SpS(i));
end

%%determining T AFstoi
hend = hinitial;                        %% enthalpy at the end is the same
hend_a = Yend*hia';


T2 = interp1(hend_a ,TR ,hinitial)
for i=1:NSp
    hi2(i)    = HNasa(T2,SpS(i));
end
h2check = Yend*hi2';

%%determining T 8*AFstoi
hend2 = hinitial2;
hend_2 = Yend2*hia';

T3 = interp1(hend_2 ,TR ,hinitial2)
for i=1:NSp
    hiend(i) = HNasa(T3,SpS(i));
end
hendcheck = Yend2*hiend';

%% the amount of air for AF 0.2-1
%plot(TR,hend_a)
%hold on
for i=1:2801
hplot(:,i)=hinitial;
end
%plot(TR,hplot)
Tloop = [T3:2:T2];
Eqloop1 = [1:0.00797:8];
for i=1:879
Eqloop(i)=Eqloop1(880-i);
end


AFloop1 = [AFstoi:0.1135:AF2];
for i=1:879
AFloop(i)=AFloop1(880-i);
end

plot(Tloop,AFloop)

%%plot(Tloop,Eqloop);


