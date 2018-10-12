% this script is made to post process the time dependent data
% received from the matlab calculation. 

clear all;
close all;
clc;

A = importdata('output.txt');
y = A.data;

% Assign some pointers
it   = find(strcmpi('Time',A.colheaders));
ix   = find(strcmpi('Position',A.colheaders));
iu   = find(strcmpi('Velocity',A.colheaders));
iT   = find(strcmpi('Temperature',A.colheaders));
iP   = find(strcmpi('Pressure',A.colheaders));
if1  = find(strcmpi('species1',A.colheaders));
if2  = find(strcmpi('species2',A.colheaders));
irho  = find(strcmpi('density',A.colheaders));


% store the values in the appropriate array
time = y(:,it);
x = y(:,ix);
u = y(:,iu);
T = y(:,iT);
P = y(:,iP);
f1 = y(:,if1);
f2 = y(:,if2);
rho = y(:,irho);

% reshape the array such that they can be plotted within
% the computational domain. 
col = length(unique(time));
row = length(u)/col;

x_new = reshape(x,[row, col]);
u_new = reshape(u,[row, col]);
T_new = reshape(T,[row, col]);
P_new = reshape(P,[row, col]);
f1_new = reshape(f1,[row, col]);
f2_new = reshape(f2,[row, col]);
rho_new = reshape(rho,[row, col]);
% plot the results
figure(1)
hold on
plot(x_new(2:row-1,:),T_new(2:row-1,:),'-k','LineWidth',2)
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(2)
hold on
plot(x_new(2:row-1,:),u_new(2:row-1,:),'-r','LineWidth',2)
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(3)
hold on
plot(x_new(2:row-1,:),P_new(2:row-1,:),'-b','LineWidth',2)
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(4)
hold on
plot(x_new(2:row-1,:),f1_new(2:row-1,:),'-b','LineWidth',2)
plot(x_new(2:row-1,:),f2_new(2:row-1,:),'-r','LineWidth',2)
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
 
figure(5)
hold on
plot(x_new(1:row-1,:),rho_new(1:row-1,:),'-k','LineWidth',2)
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)