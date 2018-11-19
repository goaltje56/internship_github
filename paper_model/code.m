% This program solves steady convection-diffustion problems
% using the simple algorithm described in ch. 6.4 in 
% "Computational Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. 
% Symbols and equations cited are from this reference unless
% mentioned otherwise

% Reference: Computational Fluid Dynamics, H.K. Versteeg and W. 
%            Malalasekera, Longman Group Ltd, 1995
% -------------------------------------------------------------------------

clear all;
close all;
clc;

% -------------------------------------------------------------------------
%% set timer to indicate the computation time
timerVal = tic;

%% initializing
global Patm Runi      
Patm        = 101325;           % athmosphesric pressure [Pa]
Runi        = 8.314;

NPI         = 100;              % number of grid cells in x-direction [-] 
XMAX        = 1;                % length of the domain [m]
% u_in        = 1.5;              % inflow velocity [m/s]
Q_in        = 300;              % flow rate [mol/s]
T           = 298;              % temperature
A           = 1;                % area of one cell [m^2]
Total_time  = 100;              % total simulation time [s]
P_k         = [30.78 5.7]*10^(-10);
alpha       = [P_k(1)/P_k(2) P_k(2)/P_k(1)] ;
X0          = [0.205 0.795];       % initial guess for Mr, Mp and Y_{1:n}
Pr          = 790.8*10^3;
Pp          = 101.3*10^3;
gamma       = Pp/Pr;
% a           = 1 - alpha;
% b           = -1+alpha+(1/gamma)+ (X0/gamma)*(alpha-1);
% c           = -alpha*X0/gamma;
for i = 1:2
    a           = 1 - alpha(i);
    b           = -1+alpha(i)+(1/gamma)+ (X0(i)/gamma)*(alpha(i)-1);
    c           = -alpha(i)*X0(i)/gamma;
    Xp_in(i)       = (-b+sqrt(b^2-4*a*c))/(2*a);
end
%% grid generation
[Dx, x, x_u] = grid_gen(NPI,XMAX);   % create staggered grid




