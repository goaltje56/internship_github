% This script is made to postprocess results 
% regarding the stability of the steady_incompressible1D code.
% Different inlet velocities are used and the 
% minimum required gridpoints are determined
% as well as the number of iterations required for convergence up 
% to 4 decimals to the exact solution and the computation time. 

clear all;
close all;
clc;

u    = [1  3   5   10];                 % velocity at inlet
NPI  = [7  43  70  173];                % minimum required grid points
iter = [47 383 607 1448];               % iterations till convergence
CPU  = [0.1105 1.2724 6.3346 151.9467]; % computation time

figure(1)
hold on
% plot(u,NPI,'bo','LineWidth',2)
% plot(u,iter,'ro','LineWidth',2)
% plot(u,CPU,'ko','LineWidth',2)
% plot(NPI,iter,'ro','LineWidth',2)
plot(NPI,CPU,'ko','LineWidth',2) %'Yscale','log','Xscale','log',
% plot(iter,CPU,'ko','LineWidth',2)
legend('CPU','Location','NorthWest')
set(gca,'Yscale','log','Xscale','log','box','on', 'LineWidth', 2, 'FontSize', 15)
grid on
xlabel('NPI [-] ','LineWidth', 2)
axis([0 700 0 300]);