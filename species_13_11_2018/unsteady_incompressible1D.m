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
%% set path to folder where values must be stored, clean the old file
path_Results1 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\species_12_11_2018\results\output.txt';
path_Results2 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\species_12_11_2018\results\diff.txt';
path_Results3 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\species_12_11_2018\results\x_dummy.txt';

% add name taggs to created file and close file again.
test = fopen(path_Results1,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', 'Time','Position', 'u_Position','density','velocity','Temperature', 'Pressure', 'species1', 'species2', 'species3','species4');
fclose(test);

% add name taggs to created file and close file again.
test = fopen(path_Results2,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n','Dspecies1', 'Dspecies2', 'Dspecies3','Dspecies4','X1', 'X2', 'X3','X4');
fclose(test);

% add name taggs to created file and close file again.
test = fopen(path_Results3,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n','Mr', 'Mp', 'Yp1', 'Yp2', 'Yp3', 'Yp4', 'stagecut');
fclose(test);
% -------------------------------------------------------------------------
%% set timer to indicate the computation time
timerVal = tic;

%% initializing
global Patm       
Patm        = 101325;           % athmosphesric pressure [Pa]

NPI         = 100;              % number of grid cells in x-direction [-] 
XMAX        = 1;                % length of the domain [m]
u_in        = 0.0015;           % inflow velocity [m/s]
A           = 1;                % area of one cell [m^2]
Total_time  = 200;              % total simulation time [s]
x0 = [10 10 0 0 0 0];       % initial guess for Mr, Mp and Y_{1:n}
% massflow    = 0;                % if 1 massflow else moleflow
Pr          = 691*10^3;
Pp          = 103*103;

% species properties some values have to be set manually!!
[rho_s, Y_k, X_k, Y_in, X_in, Y2_k, X2_k, Y2_in, X2_in, iAll, MW, MW_mix, MW2_mix, D, D_k, D2_k, P_k, f_old, f2_old, sink, n] = species_init(NPI);

% make a vector with initial values for all non-specie dependent parameters 
[u, rho, p, pc, T, mu, Cp, d_u, b, SP, Su, relax_u, relax_pc, relax_T, relax_rho, relax_f, Dt, u_old, rho_old, T_old, pc_old]...
 = param_init(NPI, u_in);

%% grid generation
[Dx, x, x_u] = grid_gen(NPI,XMAX);   % create staggered grid

store_times = 0:10:Total_time;      % define sample points to save data
ii          = 1;

%--------------------------------------------------------------------------
%% The main calculation part
for time = 0:Dt:Total_time
    for z =1:100
    [u, Y_k, p] = bound(NPI,rho,x,x_u,A,u, u_in, Y_in, Y_k, p);                          % Apply boundary condtions
    
    % momentum
    [aP_u, aE_u, aW_u, b_u, d_u, Istart_u, u]   = ...
    ucoeff(NPI, rho, x, x_u, u, p, A, relax_u, d_u, mu, u_in, Dt, u_old, Dx, rho_old);          % Determine coefficients
    u                                           = solve_eq(NPI, aE_u, aW_u, aP_u, b_u, u, 3);   % Solve equation
    [u, Y_k, p]                    = bound(NPI,rho,x,x_u,A,u, u_in, Y_in, Y_k, p);       % Apply boundary condtions

    % pressure correction (modified form of continuity equation)
    [aE_pc, aW_pc, aP_pc, b_pc, Istart_pc, pc]  = pccoeff(NPI, rho, A, x, x_u, u, d_u, pc, rho_old, Dx, Dt);    % Determine coefficients
    pc                                          = solve_eq(NPI-1, aE_pc, aW_pc, aP_pc, b_pc, pc, 2);            % Solve equation

    % correction for pressure and velocity
    [p u pc] = velcorr(NPI, pc, p, u, relax_pc, d_u);                                           % Pressure & velocity correction

%     %% Temperature
%     [aE_T, aW_T, aP_T, b_T, Istart_T, aPold_T] = Tcoeff(NPI, rho, A, x, x_u, u, T, Gamma, relax_T, Dt, T_old, Dx);
%     [TR, r_T] = GS_solve(NPI+1,T, aW_T, aE_T, aP_T, b_T, 10^(-6));
%  	  T = solve_eq(NPI,aE_T, aW_T, aP_T, b_T, T, 2);
%     
    %% Species 
    for i = 1:n
        [aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Istart_f, Y_sink(i,:)] = Fcoeff(NPI, rho, rho_s(i), A, x, x_u, u, Y_k(i,:), Y2_k(i,:), Pr, Pp, P_k(i,:), D_k(i,:), relax_f, Dt, f_old(i,:), Dx, rho_old, sink(i), i);
        Y_k(i,:) = solve_eq(NPI,aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Y_k(i,:), 2);
    end
    Y_k = species_bound(NPI, n, Y_k);
    
    [u, Y_k, p] = bound(NPI,rho,x,x_u,A,u, u_in, Y_in, Y_k, p);  

    end
    
    for i = 2:NPI+1
        [x_dummy(i,:), Xr(i,:)] = self_mass((u(2)*rho(2)), x0, sink, Y_k(:,i), Y_k(:,1), Y_sink(:,1:i), n);
    end
    Xr_new = Xr';
    Y_k(:,2:end-1) = Xr_new(:,2:end);
    [u, Y_k, p] = bound(NPI,rho,x,x_u,A,u, u_in, Y_in, Y_k, p);  

    % store results of this run as old results for next iteration
    [u_old, pc_old, rho_old, f_old] = storeresults(NPI, u, pc, rho, Y_k, u_old, pc_old, rho_old, f_old, n);
    
    % determine new value for X-k and D_k 
    [X_k, D_k]          = mole(NPI, n, Y_k, MW, D);
    
%     store data at different time steps
    if time == store_times(ii)
        time_x = time*ones(1,length(u));
        test = fopen(path_Results1,'a');
        fprintf(test,'%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n',[time_x; x; x_u; rho; u; T; p; Y_k(1,:); Y_k(2,:); Y_k(3,:); Y_k(4,:)]);    
        fprintf(test,'\n');        
        fclose(test);

        file2 = fopen(path_Results2,'a');
        fprintf(file2,'%-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f\n',[D_k(1,:); D_k(2,:); D_k(3,:); D_k(4,:); X_k(1,:); X_k(2,:); X_k(3,:); X_k(4,:)]);    
        fprintf(file2,'\n');        
        fclose(file2);
        
        file3 = fopen(path_Results3,'a');
        fprintf(file3,'%-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f\n',[x_dummy(:,1)'; x_dummy(:,2)'; x_dummy(:,3)'; x_dummy(:,4)'; x_dummy(:,5)'; x_dummy(:,6)'; x_dummy(:,2)'/(rho(5)*u(5))]);    
        fprintf(file3,'\n');        
        fclose(file3);        
        
    ii = ii +1;
    end

end

elapsedTime = toc(timerVal)

figure(1)
hold on
% plot(x,ones(1,length(x)),'bo','LineWidth',2)
% plot(x_u,1.2*ones(1,length(x_u)),'ro','LineWidth',2)
% legend('cell center','cell face','Location','NorthWest')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
axis([0 XMAX+Dx 0 1]);
p1 = plot(x(2:NPI+1),p(2:NPI+1),'b','LineWidth',2)
% plot(x(1:NPI+1),T(1:NPI+1),'k','LineWidth',2)
% plot(x_u(2:NPI+2),u(2:NPI+2),'sr','LineWidth',2);
% plot(x(1:NPI+1),T2(1:NPI+1),'r','LineWidth',2)

% plot(x(1:NPI+2),pc(1:NPI+2),'sb','LineWidth',2)
% plot(x(2:NPI+1),rho(2:NPI+1),':c','LineWidth',2)
% plot(x(2:NPI+1),d_u(2:NPI+1),':k','LineWidth',2)
% legend('P','u','P_c','\rho','d_u','Location','SouthWest')
legend(p1,'P','Location','NorthEast')

for i = 2:NPI+1
mdot(i) = rho(i)*u(i)*A;
end

figure(2)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
axis([0 XMAX+Dx 0 0.002]);
p2 = plot(x_u(2:NPI+2),u(2:NPI+2),'r','LineWidth',2);
legend(p2,'Velocity','Location','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(3)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction [-] ','LineWidth', 2)
axis([0 XMAX+Dx 0 2]);
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
ylabel('Mass fraction [-] ','LineWidth', 2)
axis([0 XMAX+Dx 0 2]);
p3 = plot(x(2:NPI+1),Xr(2:NPI+1,1),'r','LineWidth',2);
p4 = plot(x(2:NPI+1),Xr(2:NPI+1,2),'b','LineWidth',2);
p5 = plot(x(2:NPI+1),Xr(2:NPI+1,3),'k','LineWidth',2);
p6 = plot(x(2:NPI+1),Xr(2:NPI+1,4),'c','LineWidth',2);
legend([p3, p4, p5, p6],'O_2','CO2','N_2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(6)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass flow[-] ','LineWidth', 2)
% axis([0 XMAX+Dx 0 2]);
p3 = plot(x(2:NPI+1),x_dummy(2:NPI+1,1),'r','LineWidth',2);
% p4 = plot(x(1:NPI+1),-m_out(100,2:NPI+1),'b','LineWidth',2);
p5 = plot(x(2:NPI+1), x_dummy(2:NPI+1,2),'k','LineWidth',2);
legend([p3, p5],'m_{retenate}','m_{permeate}','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(7)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Normalized Permeate[-] ','LineWidth', 2)
% axis([0 XMAX+Dx 0 2]);
% p3 = plot(x(2:NPI+1),x_dummy(2:NPI+1,1),'r','LineWidth',2);
% p4 = plot(x(1:NPI+1),-m_out(100,2:NPI+1),'b','LineWidth',2);
p5 = plot(x(2:NPI+1), x_dummy(2:NPI+1,2)/(rho(5)*u(5)),'k','LineWidth',2);
legend([p3, p5],'m_{retenate}','m_{permeate}','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(8)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Stage cut[-] ','LineWidth', 2)
% axis([0 XMAX+Dx 0 2]);
% p3 = plot(x(2:NPI+1),x_dummy(2:NPI+1,1),'r','LineWidth',2);
% p4 = plot(x(1:NPI+1),-m_out(100,2:NPI+1),'b','LineWidth',2);
p5 = plot(x_dummy(2:NPI+1,2)/(rho(5)*u(5)),x_dummy(2:NPI+1,5) ,'k','LineWidth',2);
% legend([p3, p5],'m_{retenate}','m_{permeate}','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
