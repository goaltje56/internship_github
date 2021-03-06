% This program solves steady convection-diffustion problems
% using the simple algorithm described in ch. 6.4 in 
% "Computational Fluid Dynamics" by H.K. Versteeg and W. Malalasekera. 
% Symbols and equations cited are from this reference unless
% mentioned otherwise

% Reference: Computational Fluid Dynamics, H.K. Versteeg and W. 
%            Malalasekera, Longman Group Ltd, 1995

clear all;
close all;
clc;

% set timer to indicate the computation time
timerVal = tic;

% set path to place where values must be stored, clean the old file
% add name taggs to created file and close file again.
path = 'C:\Users\s137280\Documents\Master_tue\Internship\internship_github\1D_unsteady_incompressible\results\output.txt';
test = fopen(path,'w');
fprintf(test,'%-12s %-12s %-12s %-12s %-13s\n', 'Time','Position','velocity','Temperature', 'Pressure');
fclose(test);

%% initializing
NPI = 20;        % number of grid cells in x-direction [-] 
XMAX = 1;       % length of the domain [m]
P_atm = 101000; % athmosphesric pressure [Pa]
u_in = 1;      % inflow velocity [m/s]
A    = 1;       % area of one cell
m_in = 1;       % mass flow in
m_out = 1;      % mass flow out
Total_time = 2;

% make a vector with initial values for all parameters
[u, u_guess, p, pc, T, rho, mu, Cp, Gamma, d_u, b, SP, Su, relax_u, relax_pc, relax_T, relax_rho, Dt, u_old, T_old, pc_old, rho_old, p_old] = param_init(NPI, u_in);

%% grid generation
[Dx, x, x_u] = grid_gen(NPI,XMAX);   % create staggered grid


%% The main calculation part
for time = 0:Dt:Total_time
    for z =1:100
    [u, u_guess, T, m_in, m_out, p] = bound(NPI,rho,x,x_u,A,u, u_in, T, u_guess, p);
    
    % momentum
    [u_guess, d_u] = upseudo(NPI, rho, x, x_u, u, A, relax_u, d_u,mu, Dt, Dx, u_guess);

    % pressure (modified form of continuity equation)
    [aE_p, aW_p, aP_p, b_p, Istart_p, p] = pcoeff(NPI, rho, A, x, x_u, u_guess, d_u, p, p_old, Dx, Dt);
    p = solve_eq(NPI-1, aE_p, aW_p, aP_p, b_p, p, 2);

    [u, u_guess, T, m_in, m_out, p] = bound(NPI,rho,x,x_u,A,u, u_in, T, u_guess, p);
    
    % momentum
    [aP_u, aE_u, aW_u, b_u, d_u, Istart_u, u] = ucoeff(NPI, rho, x, x_u, u, p, A, relax_u, d_u, mu, u_in, T, Dt, u_old, Dx);
    u = solve_eq(NPI, aE_u, aW_u, aP_u, b_u, u, 3);

    [u, u_guess, T, m_in, m_out, p] = bound(NPI,rho,x,x_u,A,u, u_in, T, u_guess, p);

    % pressure correction (modified form of continuity equation)
    [aE_pc, aW_pc, aP_pc, b_pc, Istart_pc, pc] = pccoeff(NPI, rho, A, x, x_u, u, d_u, pc);
    pc = solve_eq(NPI-1, aE_pc, aW_pc, aP_pc, b_pc, pc, 2);

    [u, u_guess, T, m_in, m_out, p] = bound(NPI,rho,x,x_u,A,u, u_in, T, u_guess, p);
    
    % correction for pressure and velocity
    u = velcorr(NPI, pc, u, d_u);

    % Temperature
    [aE_T, aW_T, aP_T, b_T, Istart_T] = Tcoeff(NPI, rho, A, x, x_u, u, T, Gamma, relax_T, Dt, T_old, Dx);
%     [TR, r_T] = GS_solve(NPI+1,T, aW_T, aE_T, aP_T, b_T, 10^(-6));
 	T = solve_eq(NPI,aE_T, aW_T, aP_T, b_T, T, 2);

    [u, u_guess, T, m_in, m_out, p] = bound(NPI,rho,x,x_u,A,u, u_in, T, u_guess, p);    

    end
    
    % store results of this run as old restults for next iteration
    [u_old, pc_old, T_old, rho_old, p_old] = storeresults(NPI, u, pc,p, T, rho, u_old, pc_old, T_old, rho_old, p_old);
    
    % store data it different time steps
    if time < 10*Dt
        time_x = time*ones(1,length(u));
        test = fopen(path,'a');
        fprintf(test,'%-12.2f %-12.2f %-12.2f %-12.2f %-12.2f \n',[time_x; x; u; T; p]);    
        fprintf(test,'\n');        
        fclose(test);
        
    elseif mod(Total_time,time) == 0 
        time_x = time*ones(1,length(u));
        test = fopen(path,'a');
        fprintf(test,'%-12.2f %-12.2f %-12.2f %-12.2f %-12.2f \n',[time_x; x; u; T; p]);    
        fprintf(test,'\n');        
        fclose(test);
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
axis([0 XMAX+Dx 0 600]);
plot(x(2:NPI),p(2:NPI),'b','LineWidth',2)
plot(x(1:NPI+1),T(1:NPI+1),'k','LineWidth',2)
% plot(x_u(2:NPI+2),u(2:NPI+2),'sr','LineWidth',2);
% plot(x(1:NPI+1),T2(1:NPI+1),'r','LineWidth',2)

% plot(x(1:NPI+2),pc(1:NPI+2),'sb','LineWidth',2)
% plot(x(2:NPI+1),rho(2:NPI+1),':c','LineWidth',2)
% plot(x(2:NPI+1),d_u(2:NPI+1),':k','LineWidth',2)
% legend('P','u','P_c','\rho','d_u','Location','SouthWest')
legend('P','Location','NorthEast')

% for i = 2:NPI+1
% mdot(i) = rho(i)*u(i)*A;
% end

