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

% store the values in the appropriate array
time = y(:,it);
x = y(:,ix);
u = y(:,iu);
T = y(:,iT);
P = y(:,iP);

% reshape the array such that they can be plotted within
% the computational domain. 
col = length(unique(time));
row = length(u)/col;

x_new = reshape(x,[row, col]);
u_new = reshape(u,[row, col]);
T_new = reshape(T,[row, col]);
P_new = reshape(P,[row, col]);

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
 