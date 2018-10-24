% this script is made to post process the time dependent data
% received from the matlab calculation. 

clear all;
close all;
clc;

A = importdata('output.txt');
B = importdata('diff.txt');

yy = B.data;
y = A.data;

% Assign some pointers
it   = find(strcmpi('Time',A.colheaders));
ix   = find(strcmpi('Position',A.colheaders));
ix_u   = find(strcmpi('u_Position',A.colheaders));
iu   = find(strcmpi('Velocity',A.colheaders));
iT   = find(strcmpi('Temperature',A.colheaders));
iP   = find(strcmpi('Pressure',A.colheaders));
if1  = find(strcmpi('species1',A.colheaders));
if2  = find(strcmpi('species2',A.colheaders));
if3  = find(strcmpi('species3',A.colheaders));
if4  = find(strcmpi('species4',A.colheaders));
irho  = find(strcmpi('density',A.colheaders));

iD1  = find(strcmpi('Dspecies1',B.colheaders));
iD2  = find(strcmpi('Dspecies2',B.colheaders));
iD3  = find(strcmpi('Dspecies3',B.colheaders));
iD4  = find(strcmpi('Dspecies4',B.colheaders));

% store the values in the appropriate array
time = y(:,it);
x = y(:,ix);
x_u = y(:,ix_u);
u = y(:,iu);
T = y(:,iT);
P = y(:,iP);
f1 = y(:,if1);
f2 = y(:,if2);
f3 = y(:,if3);
f4 = y(:,if4);
rho = y(:,irho);

D1 = yy(:,iD1);
D2 = yy(:,iD2);
D3 = yy(:,iD3);
D4 = yy(:,iD4);

% reshape the array such that they can be plotted within
% the computational domain. 
col = length(unique(time));
row = length(u)/col;

x_new = reshape(x,[row, col]);
x_u_new = reshape(x_u,[row, col]);
u_new = reshape(u,[row, col]);
T_new = reshape(T,[row, col]);
P_new = reshape(P,[row, col]);
f1_new = reshape(f1,[row, col]);
f2_new = reshape(f2,[row, col]);
f3_new = reshape(f3,[row, col]);
f4_new = reshape(f4,[row, col]);
rho_new = reshape(rho,[row, col]);

D1_new = reshape(D1,[row, col]);
D2_new = reshape(D2,[row, col]);
D3_new = reshape(D3,[row, col]);
D4_new = reshape(D4,[row, col]);

% plot the results
% figure(1)
% hold on
% plot(x_new(2:row-1,:),T_new(2:row-1,:),'-k','LineWidth',2)
% xlabel('Geometric position [m] ','LineWidth', 2)
% ylabel('Temperature [K] ','LineWidth', 2)
% axis([0 1 0 0.002]);
% legend('T','NorthEast')
% set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(2)
hold on
grid on
plot(x_u_new(2:row-1,:),u_new(2:row-1,:),'-r','LineWidth',2)
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Velocity [m/s] ','LineWidth', 2)
axis([0 1 0 0.002]);
legend('u','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(3)
hold on
plot(x_new(2:row-1,:),P_new(2:row-1,:),'-b','LineWidth',2)
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Pressure [Pa] ','LineWidth', 2)
axis([0 1 0 1000]);
legend('P','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(4)
hold on
grid on
subplot(2,2,1)
p1 = plot(x_new(1:row-1,:),f1_new(1:row-1,:),'-b','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 0.5])
legend(p1,'O_2', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction [-] ','LineWidth', 2)

subplot(2,2,2)
p2 = plot(x_new(1:row-1,:),f2_new(1:row-1,:),'-r','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 0.001])
legend(p2,'CO_2', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction [-] ','LineWidth', 2)

subplot(2,2,3)
p3 = plot(x_new(1:row-1,:),f3_new(1:row-1,:),'-k','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 1])
legend(p3,'N_2', 'Location','NorthEast')

subplot(2,2,4)
p4 = plot(x_new(1:row-1,:),f4_new(1:row-1,:),'-c','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 1])
legend(p4,'AR', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction [-] ','LineWidth', 2)
% legend('O_2','CO_2','N_2','Ar','location','NorthEast')
 
figure(5)
hold on
plot(x_new(1:row-1,:),rho_new(1:row-1,:),'-k','LineWidth',2)
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Density [kg/m^3] ','LineWidth', 2)
axis([0 1 0 2]);
legend('\rho','Location','NorthEast')
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

%%
figure(6)
hold on
grid on
subplot(2,2,1)
p1 = plot(x_new(1:row-1,:),D1_new(1:row-1,:),'-b','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 0.5])
legend(p1,'D_{O_2}', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity ','LineWidth', 2)

subplot(2,2,2)
p2 = plot(x_new(1:row-1,:),D2_new(1:row-1,:),'-r','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 0.001])
legend(p2,'D_{CO_2}', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity','LineWidth', 2)

subplot(2,2,3)
p3 = plot(x_new(1:row-1,:),D3_new(1:row-1,:),'-k','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 1])
legend(p3,'D_{N_2}', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity','LineWidth', 2)

subplot(2,2,4)
p4 = plot(x_new(1:row-1,:),D4_new(1:row-1,:),'-c','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 1])
legend(p4,'D_{AR}', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity ','LineWidth', 2)