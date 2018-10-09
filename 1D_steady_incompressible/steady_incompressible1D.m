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
timerVal = tic;
%% initializing
NPI = 10;        % number of grid cells in x-direction [-] 
XMAX = 1;       % length of the domain [m]
P_atm = 101000; % athmosphesric pressure [Pa]
u_in = 2;      % inflow velocity [m/s]
A    = 1;       % area of one cell
m_in = 1;       % mass flow in
m_out = 1;      % mass flow out

% make a vector with initial values for all parameters
[u, p, pc, T, rho, mu, Cp, Gamma, d_u, b, SP, Su, relax_u, relax_pc, relax_T, relax_rho] = param_init(NPI, u_in);

%% grid generation
[Dx, x, x_u] = grid_gen(NPI,XMAX);   % create staggered grid


%% The main calculation part
for z =1:20
[u T m_in m_out] = bound(NPI,rho,x,x_u,A,u, u_in, T);
    
% momentum
[aP_u aE_u aW_u b_u d_u Istart_u u T] = ucoeff(NPI, rho, x, x_u, u, p, A, relax_u, d_u, mu, u_in, T);
[u r_u] = GS_solve2(NPI+1, u, aW_u, aE_u, aP_u, b_u, 10^(-6));
[u T m_in m_out] = bound(NPI,rho,x,x_u,A,u, u_in, T);
% pressure correction (modified form of continuity equation)
[aE_pc aW_pc aP_pc b_pc Istart_pc pc] = pccoeff(NPI, rho, A, x, x_u, u, d_u, pc);
[pc r_pc] = GS_solve(NPI+1,pc, aW_pc, aE_pc, aP_pc, b_pc, 10^(-6));

% correction for pressure and velocity
[p u pc] = velcorr(NPI, pc, p, u, relax_pc, d_u);

% Temperature
[aE_T aW_T aP_T b_T Istart_T] = Tcoeff(NPI, rho, A, x, x_u, u, T, Gamma, relax_T);
[T r_T] = GS_solve(NPI+1,T, aW_T, aE_T, aP_T, b_T, 10^(-1));

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
axis([0 XMAX+Dx 0 400]);
% plot(x(2:NPI+1),p(2:NPI+1),'b','LineWidth',2)
plot(x(1:NPI+1),T(1:NPI+1),'b','LineWidth',2)
% plot(x_u(2:NPI+2),u(2:NPI+2),'sr','LineWidth',2);
% plot(x(1:NPI+2),pc(1:NPI+2),'sb','LineWidth',2)
% plot(x(2:NPI+1),rho(2:NPI+1),':c','LineWidth',2)
% plot(x(2:NPI+1),d_u(2:NPI+1),':k','LineWidth',2)
% legend('P','u','P_c','\rho','d_u','Location','SouthWest')
legend('T','Location','NorthEast')

% for i = 2:NPI+1
% mdot(i) = rho(i)*u(i)*A;
% end

% 
% for i = 2: length(aP_u)
%    A(i-1,i-1) = aP_u(i);
%    B(i-1,i-1) = aP_pc(i);
%  
%    if i > 2 && i< length(aP_u)
%         A(i-1,i) = aE_u(i);
%         B(i-1,i) = aE_pc(i);  
%    end
%    if i >= 3
%        B(i-1,i-2) = aW_pc(i);
%        A(i-1,i-2) = aW_u(i);
%    end
%    b_new(i-1) = b_u(i);
%    b_newpc(i-1) = b_pc(i);
% 
% end