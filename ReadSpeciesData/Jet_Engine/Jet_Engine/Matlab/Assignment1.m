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
TA=300; % Initial temperature
cFuel={'Gasoline'};
iSp = myfind({Sp.Name},{'O2','CO2','H2O','N2','O'});                            % Proper indices to database
%% Determine compositions Yi of initial and final mixture and do your stuff

HNasa(Tref,Sp(55))
