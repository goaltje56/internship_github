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
path_Results1 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_Simplified_3membranes\results\output.txt';
path_Results2 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_Simplified_3membranes\results\diff.txt';
path_Results3 = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_Simplified_3membranes\results\x_dummy.txt';

% add name taggs to created file and close file again.
test = fopen(path_Results1,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n', 'Time','Position', 'u_Position', 'X_k1', 'X_k2', 'X2_k1', 'X2_k2', 'stagecut1');
fclose(test);

% add name taggs to created file and close file again.
test = fopen(path_Results2,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n','Position2', 'u_Position2', 'X3_k1', 'X3_k2', 'X4_k1', 'X4_k2', 'stagecut2');
fclose(test);

% add name taggs to created file and close file again.
test = fopen(path_Results3,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-12s %-12s %-12s\n','Position3', 'u_Position3', 'X5_k1', 'X5_k2', 'X6_k1', 'X6_k2', 'stagecut3');
fclose(test);
% -------------------------------------------------------------------------
%% set timer to indicate the computation time
timerVal = tic;

%% initializing
global Patm Runiv      
Patm        = 101325;           % athmosphesric pressure [Pa]
Runiv        = 8.314;

Dx          = 0.005;              % number of grid cells in x-direction [-] 
L1        = 0.88;                % length of the domain [m]
L2        = 0.265;
L3        = 0.665
% XMAX        = 0.71;                % length of the domain [m]
u_in        = 1;              % inflow velocity [m/s]
T           = 298;              % temperature
A           = 1;                % area of one cell [m^2]
Total_time  = 200;              % total simulation time [s]
Pr          = 400*10^3;
Pp          = 40*10^3;
w           = 1;
dSTART      = 50;
dEND        = 150;
z = 0;

%% grid generation
[NPI1, x1, x_u1] = grid_gen(Dx,L1);   % create staggered grid for membrane 1
[NPI2, x2, x_u2] = grid_gen(Dx,L2);   % create staggered grid for membrane 2
[NPI3, x3, x_u3] = grid_gen(Dx,L3);   % create staggered grid for membrane 3

% species properties some values have to be set manually!!
[rho_s, MW3, MW4, Y3_k, X3_k, Y3_in, X3_in, Y4_k, X4_k, Y4_in, X4_in, iAll, MW, rho3_real, rho4_real, rho3_old, rho4_old, D, D3_k, D4_k, P_k, f3_old, f4_old, sink, n] = species_init(NPI2);
[rho_s, MW5, MW6, Y5_k, X5_k, Y5_in, X5_in, Y6_k, X6_k, Y6_in, X6_in, iAll, MW, rho5_real, rho6_real, rho5_old, rho6_old, D, D5_k, D6_k, P_k, f5_old, f6_old, sink, n] = species_init(NPI3);
[rho_s, MW1, MW2, Y_k, X_k, Y_in, X_in, Y2_k, X2_k, Y2_in, X2_in, iAll, MW, rho_real, rho2_real, rho_old, rho2_old, D, D_k, D2_k, P_k, f_old, f2_old, sink, n] = species_init(NPI1);

% make a vector with initial values for all non-specie dependent parameters 
[u, u2, rho, d_u, b, SP, Su, relax_rho, relax_f, Dt, u_old] = param_init(NPI1, u_in);


%% generating disturbance profile
Tstart  = [30 60 90  120 160];
Tend    = [55 85 115 155 185];
down1   = disturbance(Tend(1), Tstart(1), Dt, Y_in(1), 0.05);
up1     = disturbance(Tend(2), Tstart(2), Dt, Y_in(1), 0.18);
down2   = disturbance(Tend(3), Tstart(3), Dt, Y_in(1), 0.08);
up2     = disturbance(Tend(4), Tstart(4), Dt, Y_in(1), 0.18);
up3     = disturbance(Tend(5), Tstart(5), Dt, Y_in(1), 0.30);

profile = [down1 up1 down2 up2 up3];
tt      = 1; 
%% storage 
store_times = [0:0.5:Total_time];      % define sample points to save data
ii          = 1;
M2 =    (rho2_real(2)*u_in) ;        % initial mass at permeate side
%--------------------------------------------------------------------------
%% The main calculation part
for time = 0:Dt:Total_time

%     if time> Tstart(tt) && time < Tend(tt)
%         z = z+1;
        Y_in(1) = 0.5 +0.5*sin(0.3*time).*cos(0.15*time).*cos(0.5*time).*cos(0.02*time).*sin(0.25*time);
        Y_in(2) = 1 - Y_in(1);
%     end
%     if time > Tend(tt) && tt < length(Tend)
%         tt =tt+1;
%     end
%-------------------------------------------------------------------------%
%% 1st membrane
%-------------------------------------------------------------------------%
     Y_in = species_bound(NPI1, n, Y_in);
    [Y_k, Y2_k, rho_real, rho2_real] = bound(NPI1, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);                          % Apply boundary condtions
    
    %% Species retenate
    for i = 1:n
        [aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Istart_f, Y_sink(i,:)] = Fcoeff(NPI1, rho_real, A, x1, x_u1, u, Y_k(i,:), Y2_k(i,:),T, rho_real, rho2_real, P_k(i,:), D_k(i,:), relax_f, Dt, f_old(i,:), Dx, rho_old, sink(i), Pp, Pr, MW1, MW2, w);
        Y_k(i,:) = solve_eq(NPI1,aE_f(i,:), aW_f(i,:), aP_f(i,:), b_f(i,:), Y_k(i,:), 2);
    end
    for j = 1:n
        for i = 1:NPI1+2
            if Y_k(j,i) <0
                Y_k(j,i) = 0;
                Y_sink(j,i) = 0;
            end               
        end
    end
    Y_k = species_bound(NPI1, n, Y_k);
    
    [Y_k, Y2_k,rho_real, rho2_real ] = bound(NPI1, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);  
    
    % determine new value for X_k MW rho_real and D_k 
    [X_k, X2_k, D_k, rho_real, rho2_real, MW1, MW2]          = mole(NPI1, n,Y_k, Y2_k, MW, D, rho_s);
    
    u2 = -sum(Y_sink)./rho_real(2:end);
    u2(1)= 0;
%     u2(2)= 0;
    u2(end+1) = u2(end);
    
    
    %% Species Permeate
    for i = 1:n
        [aE_f2(i,:), aW_f2(i,:), aP_f2(i,:), b_f2(i,:), Istart_f, Y2_sink(i,:)] = F2coeff(NPI1, rho,rho_s, A, x1, x_u1, u2, Y_k(i,:), Y2_k(i,:),T, rho_real, rho2_real, P_k(i,:), D_k(i,:), relax_f, Dt, f2_old(i,:), Dx, rho2_old, sink(i), Y_sink(i,:), sum(Y_sink), n);
        Y2_k(i,:) = solve_eq(NPI1,aE_f2(i,:), aW_f2(i,:), aP_f2(i,:), b_f2(i,:), Y2_k(i,:), 2);
    end    
    
    Y2_k = species_bound(NPI1, n, Y2_k);
    [Y_k, Y2_k,rho_real, rho2_real ]        = bound(NPI1, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);  
    
    %% determine stage cut of first membrane
    for i = 2:NPI1+2
        [x_dummy1(i,:), theta1(i)] = stagecut((rho_real(1)*u_in), Y_sink(:,1:i-1));
    end

     Y_in = species_bound(NPI1, n, Y_in);
    [Y_k, Y2_k, rho_real, rho2_real] = bound(NPI1, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);                          % Apply boundary condtions

%-------------------------------------------------------------------------%
%% 2nd membrane   
%-------------------------------------------------------------------------%

Y3_in = Y2_k(:,end);

[Y3_k, Y4_k,rho3_real, rho4_real ] = bound(NPI2, Y3_in, Y2_in, Y3_k, Y4_k, rho3_real, rho4_real);  
Y3_k = species_bound(NPI2, n, Y3_k);

    %% Species retenate
    for i = 1:n
        [aE_f3(i,:), aW_f3(i,:), aP_f3(i,:), b_f3(i,:), Istart_f3, Y_sink2(i,:)] = Fcoeff(NPI2, rho3_real, A, x2, x_u2, u, Y3_k(i,:), Y4_k(i,:),T, rho3_real, rho4_real, P_k(i,:), D3_k(i,:), relax_f, Dt, f3_old(i,:), Dx, rho_old, sink(i), Pp, Pr, MW3, MW4, w);
        Y3_k(i,:) = solve_eq(NPI2,aE_f3(i,:), aW_f3(i,:), aP_f3(i,:), b_f3(i,:), Y3_k(i,:), 2);
    end
    for j = 1:n
        for i = 1:NPI2+2
            if Y3_k(j,i) <0
                Y3_k(j,i) = 0;
                Y_sink2(j,i) = 0;
            end               
        end
    end
    Y3_k = species_bound(NPI2, n, Y3_k);
    
    [Y3_k, Y4_k,rho3_real, rho4_real ] = bound(NPI2, Y3_in, Y2_in, Y3_k, Y4_k, rho3_real, rho4_real);  
    
    % determine new value for X_k MW rho_real and D_k 
    [X3_k, X4_k, D3_k, rho3_real, rho4_real, MW3, MW4]          = mole(NPI2, n,Y3_k, Y4_k, MW, D, rho_s);
    
    u2 = -sum(Y_sink2)./rho3_real(2:end);
    u2(1)= 0;
%     u2(2)= 0;
    u2(end+1) = u2(end);
        
    %% Species Permeate
    for i = 1:n
        [aE_f4(i,:), aW_f4(i,:), aP_f4(i,:), b_f4(i,:), Istart_f, Y2_sink2(i,:)] = F2coeff(NPI2, rho,rho_s, A, x2, x_u2, u2, Y3_k(i,:), Y4_k(i,:),T, rho3_real, rho4_real, P_k(i,:), D4_k(i,:), relax_f, Dt, f4_old(i,:), Dx, rho4_old, sink(i), Y_sink2(i,:), sum(Y_sink2), n);
        Y4_k(i,:) = solve_eq(NPI2,aE_f4(i,:), aW_f4(i,:), aP_f4(i,:), b_f4(i,:), Y4_k(i,:), 2);
    end    
    
    Y4_k = species_bound(NPI2, n, Y4_k);
    [Y3_k, Y4_k,rho3_real, rho4_real ]        = bound(NPI2, Y3_in, Y4_in, Y3_k, Y4_k, rho3_real, rho4_real);  
    
    %% determine stage cut of 2nd membrane
    for i = 2:NPI2+2
        [x_dummy2(i,:), theta2(i)] = stagecut((rho3_real(1)*u_in), Y_sink2(:,1:i-1));
    end

%-------------------------------------------------------------------------%
%% 3th membrane   
%-------------------------------------------------------------------------%

Y5_in = Y4_k(:,end);

[Y5_k, Y6_k,rho5_real, rho6_real ] = bound(NPI3, Y5_in, Y6_in, Y5_k, Y6_k, rho5_real, rho6_real);  
Y5_k = species_bound(NPI3, n, Y5_k);

    %% Species retenate
    for i = 1:n
        [aE_f5(i,:), aW_f5(i,:), aP_f5(i,:), b_f5(i,:), Istart_f5, Y_sink3(i,:)] = Fcoeff(NPI3, rho5_real, A, x3, x_u3, u, Y5_k(i,:), Y6_k(i,:),T, rho5_real, rho6_real, P_k(i,:), D5_k(i,:), relax_f, Dt, f5_old(i,:), Dx, rho_old, sink(i), Pp, Pr, MW5, MW6, w);
        Y5_k(i,:) = solve_eq(NPI3,aE_f5(i,:), aW_f5(i,:), aP_f5(i,:), b_f5(i,:), Y5_k(i,:), 2);
    end
    for j = 1:n
        for i = 1:NPI3+2
            if Y5_k(j,i) <0
                Y5_k(j,i) = 0;
                Y_sink3(j,i) = 0;
            end               
        end
    end
    Y5_k = species_bound(NPI3, n, Y5_k);
    
    [Y5_k, Y6_k,rho5_real, rho6_real ] = bound(NPI3, Y5_in, Y6_in, Y5_k, Y6_k, rho5_real, rho6_real);  
    
    % determine new value for X_k MW rho_real and D_k 
    [X5_k, X6_k, D5_k, rho5_real, rho6_real, MW5, MW6]          = mole(NPI3, n,Y5_k, Y6_k, MW, D, rho_s);
    
    u2 = -sum(Y_sink3)./rho5_real(2:end);
    u2(1)= 0;
%     u2(2)= 0;
    u2(end+1) = u2(end);
        
    %% Species Permeate
    for i = 1:n
        [aE_f6(i,:), aW_f6(i,:), aP_f6(i,:), b_f6(i,:), Istart_f, Y2_sink3(i,:)] = F2coeff(NPI3, rho,rho_s, A, x3, x_u3, u2, Y5_k(i,:), Y6_k(i,:),T, rho5_real, rho6_real, P_k(i,:), D6_k(i,:), relax_f, Dt, f6_old(i,:), Dx, rho6_old, sink(i), Y_sink3(i,:), sum(Y_sink3), n);
        Y6_k(i,:) = solve_eq(NPI3,aE_f6(i,:), aW_f6(i,:), aP_f6(i,:), b_f6(i,:), Y6_k(i,:), 2);
    end    
    
    Y6_k = species_bound(NPI3, n, Y6_k);
    [Y5_k, Y6_k,rho5_real, rho6_real ]        = bound(NPI3, Y5_in, Y6_in, Y5_k, Y6_k, rho5_real, rho6_real);  
    
    %% determine stage cut of 3th membrane
    for i = 2:NPI3+2
        [x_dummy3(i,:), theta3(i)] = stagecut((rho5_real(1)*u_in), Y_sink3(:,1:i-1));
    end
    
%% updating + storing data    
   [X_k, X2_k, D_k, rho_real, rho2_real, MW1, MW2]          = mole(NPI1, n,Y_k, Y2_k, MW, D, rho_s);
   [X3_k, X4_k, D3_k, rho3_real, rho4_real, MW3, MW4]          = mole(NPI2, n,Y3_k, Y4_k, MW, D, rho_s);
   [X5_k, X6_k, D5_k, rho5_real, rho6_real, MW5, MW6]          = mole(NPI3, n,Y5_k, Y6_k, MW, D, rho_s);

    % store results of this run as old results for next iteration
    [rho_old, rho2_old, f_old, f2_old] = storeresults(NPI1, rho_real, rho2_real, Y_k, Y2_k, rho_old, rho2_old, f_old, f2_old, n);
    [rho3_old, rho4_old, f3_old, f4_old] = storeresults(NPI2, rho3_real, rho4_real, Y3_k, Y4_k, rho3_old, rho4_old, f3_old, f4_old, n);
    [rho5_old, rho6_old, f5_old, f6_old] = storeresults(NPI3, rho5_real, rho6_real, Y5_k, Y6_k, rho5_old, rho6_old, f5_old, f6_old, n);
    
    % determine new value for X_k rho_real and D_k 
%     [X_k, X2_k,D_k, rho_real, rho2_real]          = mole(NPI1, n, Y_k, Y2_k, MW, D, rho_s);
%     [Y_k, Y2_k,rho_real, rho2_real ]        = bound(NPI1, Y_in, Y2_in, Y_k, Y2_k, rho_real, rho2_real);  

%     store data at different time steps
    if time == store_times(ii)
        time_x = time*ones(1,length(u));
        test = fopen(path_Results1,'a');
        fprintf(test,'%-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f \n',[time_x; x1; x_u1; X_k(1,:); X_k(2,:); X2_k(1,:); X2_k(2,:); theta1]);    
        fprintf(test,'\n');        
        fclose(test);

        file2 = fopen(path_Results2,'a');
        fprintf(file2,'%-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f \n',[x2; x_u2; X3_k(1,:); X3_k(2,:); X4_k(1,:); X4_k(2,:); theta2]);    
        fprintf(file2,'\n');        
        fclose(file2);
        
        file3 = fopen(path_Results3,'a');
        fprintf(file2,'%-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f %-12.12f \n',[x3; x_u3; X5_k(1,:); X5_k(2,:); X6_k(1,:); X6_k(2,:); theta3]);    
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
axis([0 L1+Dx 0 2]);
p2 = plot(x_u1(2:NPI1+2),u(2:NPI1+2),'r','LineWidth',2);
legend(p2,'Velocity','Location','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(3)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction retenate [-] ','LineWidth', 2)
axis([0 L1+Dx 0 1]);
p3 = plot(x1(1:NPI1+1),Y_k(1,1:NPI1+1),'r','LineWidth',2);
p4 = plot(x1(1:NPI1+1),Y_k(2,1:NPI1+1),'b','LineWidth',2);
% p5 = plot(x1(1:NPI1+1),Y_k(3,1:NPI1+1),'k','LineWidth',2);
% p6 = plot(x1(1:NPI1+1),Y_k(4,1:NPI1+1),'c','LineWidth',2);
legend([p3, p4],'CO2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(4)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction permeate[-] ','LineWidth', 2)
axis([0 L1+Dx 0 2]);
p3 = plot(x1(2:NPI1+1),Y2_k(1,2:NPI1+1) ,'r','LineWidth',2);
p4 = plot(x1(2:NPI1+1),Y2_k(2,2:NPI1+1) ,'b','LineWidth',2);
% p5 = plot(x1(2:NPI1+1),Y2_k(3,2:NPI1+1) ,'k','LineWidth',2);
% p6 = plot(x1(2:NPI1+1),Y2_k(4,2:NPI1+1) ,'c','LineWidth',2);
legend([p3, p4],'CO2','AR','NorthEast')
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
% 
% figure(8)
% hold on
% grid on
% xlabel('Stage cut, \Theta [-] ','LineWidth', 2)
% ylabel('X_{permeate} [-] ','LineWidth', 2)
% % axis([0 0.3 0 1]);
% ylim([0 1])
% p11 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X2_k(1,2:NPI+2),'r','LineWidth',2);
% p12 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X2_k(2,2:NPI+2),'b','LineWidth',2);
% p13 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X2_k(3,2:NPI+2),'k','LineWidth',2);
% p14 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X2_k(4,2:NPI+2),'c','LineWidth',2);
% set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% % legend([p11, p12, p13, p14],'O_2','CO2','H_2O','AR','Location','Best')
% 
% figure(9)
% hold on
% grid on
% xlabel('Stage cut, \Theta [-] ','LineWidth', 2)
% ylabel('X_{retentate} [-] ','LineWidth', 2)
% ylim([0 1])
% % axis([0 0.3 0 1]);
% p11 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X_k(1,2:NPI+2),'r','LineWidth',2);
% p12 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X_k(2,2:NPI+2),'b','LineWidth',2);
% p13 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X_k(3,2:NPI+2),'k','LineWidth',2);
% p14 = plot((x_dummy(1:NPI+1,2)')/(rho_real(1)*u_in), X_k(4,2:NPI+2),'c','LineWidth',2);
% set(gca,'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% % legend([p11, p12, p13, p14],'O_2','CO2','H_2O','AR','Location','Best')

% Mp = Y_sink(2,:)+Y_sink(4,:);

figure(11)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction retenate [-] ','LineWidth', 2)
axis([0 L2+Dx 0 1]);
p3 = plot(x2(1:NPI2+1),X3_k(1,1:NPI2+1),'r','LineWidth',2);
p4 = plot(x2(1:NPI2+1),X3_k(2,1:NPI2+1),'b','LineWidth',2);
title('Membrane 2')
legend([p3, p4],'CO2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(12)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction permeate [-] ','LineWidth', 2)
axis([0 L2+Dx 0 1]);
p3 = plot(x2(1:NPI2+1),X4_k(1,1:NPI2+1),'r','LineWidth',2);
p4 = plot(x2(1:NPI2+1),X4_k(2,1:NPI2+1),'b','LineWidth',2);
legend([p3, p4],'CO2','AR','NorthEast')
title('Membrane 2')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(13)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction retenate [-] ','LineWidth', 2)
axis([0 L1+Dx 0 1]);
p3 = plot(x1(1:NPI1+1),X_k(1,1:NPI1+1),'r','LineWidth',2);
p4 = plot(x1(1:NPI1+1),X_k(2,1:NPI1+1),'b','LineWidth',2);
title('Membrane 1')
legend([p3, p4],'CO2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(14)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction permeate [-] ','LineWidth', 2)
axis([0 L1+Dx 0 1]);
p3 = plot(x1(1:NPI1+1),X2_k(1,1:NPI1+1),'r','LineWidth',2);
p4 = plot(x1(1:NPI1+1),X2_k(2,1:NPI1+1),'b','LineWidth',2);
title('Membrane 1')
legend([p3, p4],'CO2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(15)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction retenate [-] ','LineWidth', 2)
axis([0 L3+Dx 0 1]);
p3 = plot(x3(1:NPI3+1),X5_k(1,1:NPI3+1),'r','LineWidth',2);
p4 = plot(x3(1:NPI3+1),X5_k(2,1:NPI3+1),'b','LineWidth',2);
title('Membrane 3')
legend([p3, p4],'CO2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(16)
hold on
grid on
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction permeate [-] ','LineWidth', 2)
axis([0 L3+Dx 0 1]);
p3 = plot(x3(1:NPI3+1),X6_k(1,1:NPI3+1),'r','LineWidth',2);
p4 = plot(x3(1:NPI3+1),X6_k(2,1:NPI3+1),'b','LineWidth',2);
title('Membrane 3')
legend([p3, p4],'CO2','AR','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)