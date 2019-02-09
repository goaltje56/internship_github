% This program solves unsteady 1D gas separation problems for n species
% under isobaric, isothermal and frictionless circumstances.
% Symbols and equations are cited from "Computational Fluid Dynamics" 
% by H.K. Versteeg and W. Malalasekera unless
% mentioned otherwise

% Reference: Computational Fluid Dynamics, H.K. Versteeg and W. 
%            Malalasekera, Longman Group Ltd, 1995
% Author   : Catharina van Gool 
% -------------------------------------------------------------------------

clear all;
close all;
clc;

% -------------------------------------------------------------------------
%% set path to folder where values must be stored, clean the old file
path_Results1 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_simplified_species_v02_02_2019\results\retenate.txt';
path_Results2 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_simplified_species_v02_02_2019\results\permeate.txt';

% add name taggs to created file and close file again.
test = fopen(path_Results1,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s \n', 'Time','Position', 'u_Position', 'X1', 'X2', 'X3','X4');
fclose(test);

% add name taggs to created file and close file again.
test = fopen(path_Results2,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', 'X2_1', 'X2_2', 'X2_3', 'X2_4', 'stagecut');
fclose(test);
% -------------------------------------------------------------------------
%% set timer to indicate the computation time
timerVal = tic;

%% initializing
global Patm Runiv      
Patm        = 101325;       % athmosphesric pressure                [Pa]
Runiv        = 8.314;       % universal gas constant                [J/mol.K]

NPI         = 200;          % number of grid cells in x-direction   [-] 
XMAX        = 1;            % length of the domain                  [m]
u_in        = 1;            % inflow velocity                       [m/s]
T           = 298;          % temperature                           [K]
A           = 1;            % area of one cell                      [m^2]
Total_time  = 200;          % total simulation time                 [s]
Pr          = 400*10^3;     % pressure retenate                     [Pa]
Pp          = 40*10^3;      % pressure permeate                     [Pa]
w           = 1;            % width of cell domain                  [m]

% species properties some values have to be set manually!!
[rho_s, rho,Gamma, Gamma_k, MW1, MW2, Y_k, X_k, Y_in, X_in, Y2_k, X2_k, Y2_in, X2_in, iAll, MW, rho_real, rho2_real, rho_old, rho2_old, D, D_k, D2_k, P_k, f_old, f2_old, sink, n] = species_init(NPI);

% make a vector with initial values for all non-specie dependent parameters 
[u, u2, relax_f, Dt] = param_init(NPI, u_in);

%% grid generation
[Dx, x, x_u] = grid_gen(NPI,XMAX);      % create staggered grid

store_times = [0:0.5:Total_time];       % define sample points to save data
ii          = 1;

%--------------------------------------------------------------------------
%% The main calculation part
for time = 0:Dt:Total_time

     Y_in = species_bound(NPI, n, Y_in);
    [Y_k Y2_k rho_real rho2_real] = bound(NPI, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);                          % Apply boundary condtions
    
    %% Species retenate
    for i = 1:n
        [aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Istart_f, Y_sink(i,:)] = Fcoeff(NPI, rho_real, A, x, x_u, u, Y_k(i,:), Y2_k(i,:),T, rho_real, rho2_real, P_k(i,:), D_k(i,:), relax_f, Dt, f_old(i,:), Dx, rho_old, sink(i), Pp, Pr, MW1, MW2, w);
        Y_k(i,:) = solve_eq(NPI,aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Y_k(i,:), 2);
    end
% make sure the mass frations won't become negative    
    for j = 1:n
        for i = 1:NPI+1
            if Y_k(j,i) <0 %|| sum(-Y_sink(j,1:i))>rho(1)*u_in
                Y_k(j,i) = 0;
                Y_sink(j,i) = 0;
            end               
        end
    end
    
    Y_k = species_bound(NPI, n, Y_k);
    
    [Y_k, Y2_k,rho_real, rho2_real ] = bound(NPI, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);  
    
    % determine new value for X_k MW rho_real and D_k 
    [X_k, X2_k, Gamma_k, D_k, rho_real, rho2_real, MW1, MW2]          = mole(NPI, n,Gamma, Y_k, Y2_k, MW, D, rho_s);
    
    % Due to assumption rho = constant:
    u2 = -sum(Y_sink)./rho_real(2:end);
    u2(1)= 0;     % no inlet velocity

    u2(end+1) = u2(end);
    
    
    %% Species Permeate
    for i = 1:n
        [aE_f2(i,:), aW_f2(i,:), aP_f2(i,:), b_f2(i,:), Istart_f, Y2_sink(i,:)] = F2coeff(NPI, rho,rho_s, A, x, x_u, u2, Y_k(i,:), Y2_k(i,:),T, rho_real, rho2_real, P_k(i,:), D_k(i,:), relax_f, Dt, f2_old(i,:), Dx, rho2_old, sink(i), Y_sink(i,:), sum(Y_sink), n);
        Y2_k(i,:) = solve_eq(NPI,aE_f2(i,:), aW_f2(i,:), aP_f2(i,:), b_f2(i,:), Y2_k(i,:), 2);
    end    
    
    Y2_k = species_bound(NPI, n, Y2_k);
    [Y_k, Y2_k,rho_real, rho2_real ]        = bound(NPI, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);  
    
    %% determine stage cut
    for i = 2:NPI+2
        [x_dummy(i,:), theta(i)] = stagecut((rho_real(1)*u_in), Y_sink(:,1:i-1));
    end

   [X_k, X2_k, Gamma_k, D_k, rho_real, rho2_real, MW1, MW2]          = mole(NPI, n,Gamma, Y_k, Y2_k, MW, D, rho_s);

    % store results of this run as old results for next iteration
    [rho_old, rho2_old, f_old, f2_old] = storeresults(NPI, rho_real, rho2_real, Y_k, Y2_k, rho_old, rho2_old, f_old, f2_old, n);
    
    % determine new value for X_k rho_real and D_k 
    [X_k, X2_k, Gamma_k, D_k, rho_real, rho2_real]          = mole(NPI, n, Gamma, Y_k, Y2_k, MW, D, rho_s);
    [Y_k, Y2_k,rho_real, rho2_real ]        = bound(NPI, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);  

%     store data at different time steps
    if time == store_times(ii)
        time_x = time*ones(1,length(u));
        file1 = fopen(path_Results1,'a');
        fprintf(file1,'%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f \n',[time_x; x; x_u; X_k(1,:); X_k(2,:); X_k(3,:); X_k(4,:)]);    
        fprintf(file1,'\n');        
        fclose(file1);
        
        file2 = fopen(path_Results2,'a');
        fprintf(file2,'%-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f\n',[X2_k(1,:); X2_k(2,:); X2_k(3,:); X2_k(4,:); theta]);    
        fprintf(file2,'\n');        
        fclose(file2);        
        
    ii = ii +1;
    end

end

elapsedTime = toc(timerVal)


figure(2)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Velocity [m/s] ','LineWidth', 2)
axis([0 XMAX+Dx 0 2]);
p2 = plot(x_u(2:NPI+2),u(2:NPI+2),'r','LineWidth',2);
legend(p2,'Velocity','Location','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(3)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction retenate [-] ','LineWidth', 2)
axis([0 XMAX+Dx 0 1]);
p3 = plot(x(1:NPI+1),Y_k(1,1:NPI+1),'r','LineWidth',2);
p4 = plot(x(1:NPI+1),Y_k(2,1:NPI+1),'b','LineWidth',2);
p5 = plot(x(1:NPI+1),Y_k(3,1:NPI+1),'k','LineWidth',2);
p6 = plot(x(1:NPI+1),Y_k(4,1:NPI+1),'c','LineWidth',2);
legend([p3, p4, p5, p6],'O_2','CO2','N_2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(4)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction permeate[-] ','LineWidth', 2)
axis([0 XMAX+Dx 0 2]);
p3 = plot(x(2:NPI+1),Y2_k(1,2:NPI+1) ,'r','LineWidth',2);
p4 = plot(x(2:NPI+1),Y2_k(2,2:NPI+1) ,'b','LineWidth',2);
p5 = plot(x(2:NPI+1),Y2_k(3,2:NPI+1) ,'k','LineWidth',2);
p6 = plot(x(2:NPI+1),Y2_k(4,2:NPI+1) ,'c','LineWidth',2);
legend([p3, p4, p5, p6],'O_2','CO2','N_2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

% figure(6)
% hold on
% grid on
% xlabel('Geometric position [m] ','LineWidth', 2)
% ylabel('Mass flow[-] ','LineWidth', 2)
% % axis([0 XMAX+Dx 0 2]);
% p3 = plot(x(2:NPI+1),x_dummy(2:NPI+1,1),'r','LineWidth',2);
% % p4 = plot(x(1:NPI+1),-m_out(100,2:NPI+1),'b','LineWidth',2);
% p5 = plot(x(2:NPI+1), x_dummy(2:NPI+1,3),'k','LineWidth',2);
% legend([p3, p5],'m_{retenate}','m_{permeate}','NorthEast')
% set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(8)
hold on
grid on
xlabel('Stage cut, \Theta [-] ','LineWidth', 2)
ylabel('X_{permeate} [-] ','LineWidth', 2)
% axis([0 0.3 0 1]);
ylim([0 1])
p11 = plot(theta(2:NPI+2), X2_k(1,2:NPI+2),'r','LineWidth',2);
p12 = plot(theta(2:NPI+2), X2_k(2,2:NPI+2),'b','LineWidth',2);
p13 = plot(theta(2:NPI+2), X2_k(3,2:NPI+2),'k','LineWidth',2);
p14 = plot(theta(2:NPI+2), X2_k(4,2:NPI+2),'c','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% legend([p11, p12, p13, p14],'O_2','CO2','H_2O','AR','Location','Best')

figure(9)
hold on
grid on
xlabel('Stage cut, \Theta [-] ','LineWidth', 2)
ylabel('X_{retentate} [-] ','LineWidth', 2)
ylim([0 1])
% axis([0 0.3 0 1]);
p11 = plot(theta(2:NPI+2), X_k(1,2:NPI+2),'r','LineWidth',2);
p12 = plot(theta(2:NPI+2), X_k(2,2:NPI+2),'b','LineWidth',2);
p13 = plot(theta(2:NPI+2), X_k(3,2:NPI+2),'k','LineWidth',2);
p14 = plot(theta(2:NPI+2), X_k(4,2:NPI+2),'c','LineWidth',2);
set(gca,'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% legend([p11, p12, p13, p14],'O_2','CO2','H_2O','AR','Location','Best')

% Mp = Y_sink(2,:)+Y_sink(4,:);

figure(10)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mole fraction retenate [-] ','LineWidth', 2)
axis([0 XMAX+Dx 0 1]);
p1 = plot(x(1:NPI+1),X_k(1,1:NPI+1),'r','LineWidth',2);
p2 = plot(x(1:NPI+1),X_k(2,1:NPI+1),'b','LineWidth',2);
p3 = plot(x(1:NPI+1),X_k(3,1:NPI+1),'k','LineWidth',2);
p4 = plot(x(1:NPI+1),X_k(4,1:NPI+1),'c','LineWidth',2);
title('Membrane 1')
legend([p1, p2, p3, p4],'O2','CO2','H2O','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(11)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mole fraction permeate [-] ','LineWidth', 2)
axis([0 XMAX+Dx 0 1]);
p1 = plot(x(1:NPI+1),X2_k(1,1:NPI+1),'r','LineWidth',2);
p2 = plot(x(1:NPI+1),X2_k(2,1:NPI+1),'b','LineWidth',2);
p3 = plot(x(1:NPI+1),X2_k(3,1:NPI+1),'k','LineWidth',2);
p4 = plot(x(1:NPI+1),X2_k(4,1:NPI+1),'c','LineWidth',2);
title('Membrane 1')
legend([p1, p2, p3, p4],'O2','CO2','H2O','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)