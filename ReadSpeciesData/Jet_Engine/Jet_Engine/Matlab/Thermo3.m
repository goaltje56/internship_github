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
Y3 = Y2;

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
j = [0.2:0.0002857:1];              %%the step to vary phi
for i=1:NSp
Mair1i(:,i) = j.*AFstoi.*Yair(i);       %%total mass of air for 0.2-1 AFR
end 
p1=size(Mair1i);
p=p1(1,1);
for n=1:p
    for i=1:NSp
        Yfueli(n,i)=Yfuel(i);           %%all mass fraction of gasoline for 0.2-1 AFR
    end
end

for n=1:p
for i=1:NSp
    Mtotali(n,i)=Yfueli(n,i)+Mair1i(n,i);    %%total of all mass per specie
end
end
Mtotal=sum(Mtotali,2);                       %% mass of the mix

for n=1:p
for i=1:NSp
Ytotal1(n,i)=Mtotali(n,i)./Mtotal(n);           %% fraction of the mix
end
end

for n=1:p
for i=1:NSp
htotal(n,i) = Ytotal1(n,i).*hi(i);              %%the enthalpy of every specie 
end
end

htotal1=sum(htotal,2);

htotal2=htotal1;

%% determine fractions after reaction
k=ones(size(j));
k=k-j;
k=transpose(k);

%% massa verhouding zuurstof staat tot CO2 en H2O
Mo=[0 0 341.44 118.0048 0];
Mo1= Mo./353.12;

%% molverhouding na bij volledige verbranding
Ntotaal1= [0 0 7.76 6.55 11.035*0.21/0.79];
Ntotaal = sum(Ntotaal1);
Xtotaal = Ntotaal1./Ntotaal;
Mna1= Xtotaal*Mi';




for i=1:NSp
Mafteri(:,i) = j.*AFstoi.*Mo1(i);       %%total mass of CO2  for 0.2-1 AFR
end


Mfo= Mtotali;
Mfo(:,3:5)=0;
Mfo(:,2)=0;

for n =1:p
  for i=1:NSp
Madd(n,i) = Mfo(n,i).*k(n);
    end
end

Mafter=Madd+Mafteri;                    %% total mass of every specie
Mtotal2=sum(Mafter,2);                  %% total mass after

for n=1:p
for i=1:NSp
Yafteri(n,i)=Mafter(n,i)./Mtotal2(n);           %% fraction of the mix
end
end

for n=1:p
for i=1:NSp
hafteri(n,i) = Yafteri(n,i)*hia(n,i)';
end
end

hafter=sum(hafteri,2);


for n=1:2801
T2(n) = interp1(hafter,TR ,htotal1(n));
end

for n=1:2801
for i=1:NSp
    hi2a(n,i)    = HNasa(T2(n),SpS(i));
end
end

for n=1:2801
    for i=1:NSp
h2checka(n,i) = Yafteri(n,i).*hi2a(n,i);
    end
end


h2checkb = sum(h2checka,2);







