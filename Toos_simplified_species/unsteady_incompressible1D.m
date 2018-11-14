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
path_Results1 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_simplified_species\results\output.txt';
path_Results2 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_simplified_species\results\diff.txt';
path_Results3 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_simplified_species\results\x_dummy.txt';

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
global Patm Runi      
Patm        = 101325;           % athmosphesric pressure [Pa]
Runi        = 8.314;

NPI         = 100;              % number of grid cells in x-direction [-] 
XMAX        = 1;                % length of the domain [m]
u_in        = 0.0025;           % inflow velocity [m/s]
T           = 278;              % temperature
A           = 1;                % area of one cell [m^2]
Total_time  = 200;              % total simulation time [s]
x0 = [0 0 0 0 0 0];       % initial guess for Mr, Mp and Y_{1:n}
Pr          = 691*10^3;
Pp          = 103*103;

% species properties some values have to be set manually!!
[rho_s, Y_k, X_k, Y_in, X_in, Y2_k, X2_k, Y2_in, X2_in, iAll, MW, rho_real, rho2_real, rho_old, rho2_old, D, D_k, D2_k, P_k, f_old, f2_old, sink, n] = species_init(NPI);

% make a vector with initial values for all non-specie dependent parameters 
[u, d_u, b, SP, Su, relax_rho, relax_f, Dt, u_old] = param_init(NPI, u_in);

%% grid generation
[Dx, x, x_u] = grid_gen(NPI,XMAX);   % create staggered grid

store_times = 0:10:Total_time;      % define sample points to save data
ii          = 1;

%--------------------------------------------------------------------------
%% The main calculation part
for time = 0:Dt:Total_time
    for z =1:100
    Y_k = bound(NPI, Y_in, Y_k);                          % Apply boundary condtions
    
    %% Species 
    for i = 1:n
        [aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Istart_f, Y_sink(i,:)] = Fcoeff(NPI, rho_real, rho_s(i), A, x, x_u, u, Y_k(i,:), Y2_k(i,:),T, rho_real, rho2_real, P_k(i,:), D_k(i,:), relax_f, Dt, f_old(i,:), Dx, rho_old, sink(i), i);
        Y_k(i,:) = solve_eq(NPI,aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Y_k(i,:), 2);
    end
    Y_k = species_bound(NPI, n, Y_k);
    
    Y_k = bound(NPI, Y_in, Y_k);  
    
    % determine new value for X_k rho_real and D_k 
    [X_k, D_k, rho_real, rho2_real]          = mole(NPI, n, Y_k, Y2_k, MW, D, rho_s);
    
    for i = 2:NPI+1
        [x_dummy(i,:), Xr(i,:)] = self_mass((rho_real(1)*u_in), x0, sink, Y_k(:,i), Y_k(:,1), Y_sink(:,1:i), n);
    end

%     for i = 2:NPI+1
%         [Yp(i,:) Mp(i,:)] = mass_permeate(Y_k(:,1), Y_k(:,i), rho_real(1)*u_in, rho_real(i), u_in)
%     end
   
    % correct the mass fractions at retenate side
    Xr_new = Xr';
    Y_k(:,2:end-1) = Xr_new(:,2:end);
    Y_k = bound(NPI,Y_in, Y_k);  
    Y_k = species_bound(NPI, n, Y_k);

    %% Permeate Species 
    for i = 1:n
        [aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Istart_f, Y_source(i,:)] = F2coeff(NPI, rho2_real, rho_s(i), A, x, x_u, u, Y_k(i,:), Y2_k(i,:),T, rho_real, rho2_real, P_k(i,:), D_k(i,:), relax_f, Dt, f2_old(i,:), Dx, rho2_old, sink(i), Y_sink(i,:), i);
        Y2_k(i,:) = solve_eq(NPI,aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Y2_k(i,:), 2);
    end
    Y2_k = species_bound(NPI, n, Y2_k);
    end
    % store results of this run as old results for next iteration
    [rho_old, rho2_old, f_old, f2_old] = storeresults(NPI, rho_real, rho2_real, Y_k, Y2_k, rho_old, rho2_old, f_old, f2_old, n);
    
    % determine new value for X_k rho_real and D_k 
    [X_k, D_k, rho_real, rho2_real]          = mole(NPI, n, Y_k, Y2_k, MW, D, rho_s);
    
% %     store data at different time steps
%     if time == store_times(ii)
%         time_x = time*ones(1,length(u));
%         test = fopen(path_Results1,'a');
%         fprintf(test,'%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n',[time_x; x; x_u; rho; u; T; p; Y_k(1,:); Y_k(2,:); Y_k(3,:); Y_k(4,:)]);    
%         fprintf(test,'\n');        
%         fclose(test);
% 
%         file2 = fopen(path_Results2,'a');
%         fprintf(file2,'%-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f\n',[D_k(1,:); D_k(2,:); D_k(3,:); D_k(4,:); X_k(1,:); X_k(2,:); X_k(3,:); X_k(4,:)]);    
%         fprintf(file2,'\n');        
%         fclose(file2);
%         
%         file3 = fopen(path_Results3,'a');
%         fprintf(file3,'%-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f\n',[x_dummy(:,1)'; x_dummy(:,2)'; x_dummy(:,3)'; x_dummy(:,4)'; x_dummy(:,5)'; x_dummy(:,6)'; x_dummy(:,2)'/(rho(5)*u(5))]);    
%         fprintf(file3,'\n');        
%         fclose(file3);        
%         
%     ii = ii +1;
%     end

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
p5 = plot(x_dummy(2:NPI+1,2)/(rho_real(1)*u_in),x_dummy(2:NPI+1,5) ,'k','LineWidth',2);
% legend([p3, p5],'m_{retenate}','m_{permeate}','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
