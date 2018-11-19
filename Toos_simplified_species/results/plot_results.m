% this script is made to post process the time dependent data
% received from the matlab calculation. 

clear all;
close all;
clc;

A = importdata('output.txt');%_NPI100_time1600.txt');
B = importdata('diff.txt');%_NPI100_time1600.txt');
C = importdata('x_dummy.txt');%_NPI100_time1600.txt');

yy = B.data;
y = A.data;
yyy = C.data;

% Assign some pointers
it   = find(strcmpi('Time',A.colheaders));
ix   = find(strcmpi('Position',A.colheaders));
% ix_u   = find(strcmpi('u_Position',A.colheaders));
% iu   = find(strcmpi('Velocity',A.colheaders));
% iT   = find(strcmpi('Temperature',A.colheaders));
% iP   = find(strcmpi('Pressure',A.colheaders));
if1  = find(strcmpi('species1',A.colheaders));
if2  = find(strcmpi('species2',A.colheaders));
if3  = find(strcmpi('species3',A.colheaders));
if4  = find(strcmpi('species4',A.colheaders));
% irho  = find(strcmpi('density',A.colheaders));

iD1  = find(strcmpi('Dspecies1',B.colheaders));
iD2  = find(strcmpi('Dspecies2',B.colheaders));
iD3  = find(strcmpi('Dspecies3',B.colheaders));
iD4  = find(strcmpi('Dspecies4',B.colheaders));
iX1  = find(strcmpi('Yp1',B.colheaders));
iX2  = find(strcmpi('Yp2',B.colheaders));
iX3  = find(strcmpi('Yp3',B.colheaders));
iX4  = find(strcmpi('Yp4',B.colheaders));

iMr  = find(strcmpi('Mr',C.colheaders));
iMp  = find(strcmpi('Mp',C.colheaders));
iYr1 = find(strcmpi('permated1',C.colheaders));
iYr2 = find(strcmpi('permated2',C.colheaders));
iYr3 = find(strcmpi('permated3',C.colheaders));
iYr4 = find(strcmpi('permated4',C.colheaders));

% store the values in the appropriate array
time = y(:,it);
x = y(:,ix);
% x_u = y(:,ix_u);
% u = y(:,iu);
% T = y(:,iT);
% P = y(:,iP);
f1 = y(:,if1);
f2 = y(:,if2);
f3 = y(:,if3);
f4 = y(:,if4);
% rho = y(:,irho);

D1 = yy(:,iD1);
D2 = yy(:,iD2);
D3 = yy(:,iD3);
D4 = yy(:,iD4);
X1 = yy(:,iX1);
X2 = yy(:,iX2);
X3 = yy(:,iX3);
X4 = yy(:,iX4);

Mr  = yyy(:,iMr);
Mp  = yyy(:,iMp);
Yr1 = yyy(:,iYr1);
Yr2 = yyy(:,iYr2);
Yr3 = yyy(:,iYr3);
Yr4 = yyy(:,iYr4);

% reshape the array such that they can be plotted within
% the computational domain. 
col = length(unique(time));
row = length(f1)/col;
row2 = length(Mr)/col;

x_new = reshape(x,[row, col]);
% x_u_new = reshape(x_u,[row, col]);
% u_new = reshape(u,[row, col]);
% T_new = reshape(T,[row, col]);
% P_new = reshape(P,[row, col]);
f1_new = reshape(f1,[row, col]);
f2_new = reshape(f2,[row, col]);
f3_new = reshape(f3,[row, col]);
f4_new = reshape(f4,[row, col]);
% rho_new = reshape(rho,[row, col]);

D1_new = reshape(D1,[row, col]);
D2_new = reshape(D2,[row, col]);
D3_new = reshape(D3,[row, col]);
D4_new = reshape(D4,[row, col]);
X1_new = reshape(X1,[row, col]);
X2_new = reshape(X2,[row, col]);
X3_new = reshape(X3,[row, col]);
X4_new = reshape(X4,[row, col]);

Mr_new  = reshape(Mr,[row2, col]);
Mp_new  = reshape(Mp,[row2, col]);
Yr1_new = reshape(Yr1,[row2, col]);
Yr2_new = reshape(Yr2,[row2, col]);
Yr3_new = reshape(Yr3,[row2, col]);
Yr4_new = reshape(Yr4,[row2, col]);
% plot the results
% figure(1)
% hold on
% plot(x_new(2:row-1,:),T_new(2:row-1,:),'-k','LineWidth',2)
% xlabel('Geometric position [m] ','LineWidth', 2)
% ylabel('Temperature [K] ','LineWidth', 2)
% axis([0 1 0 0.002]);
% legend('T','NorthEast')
% set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

% figure(2);
% hold on
% grid on
% plot(x_u_new(2:row-1,:),u_new(2:row-1,:),'-r','LineWidth',2);
% xlabel('Geometric position [m] ','LineWidth', 2)
% ylabel('Velocity [m/s] ','LineWidth', 2)
% axis([0 1 0 0.002]);
% legend('u','Location','NorthEast')
% set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

% figure(3);
% hold on
% plot(x_new(2:row-1,:),P_new(2:row-1,:),'-b','LineWidth',2);
% xlabel('Geometric position [m] ','LineWidth', 2)
% ylabel('Pressure [Pa] ','LineWidth', 2)
% axis([0 1 0 1]);
% legend('P','Location','NorthEast')
% set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

figure(4);
hold on
grid on
subplot(2,2,1);
p1 = plot(x_new(1:row-1,:),f1_new(1:row-1,:),'-b','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 1])
legend(p1,'O_2', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction [-] ','LineWidth', 2)

subplot(2,2,2)
p2 = plot(x_new(1:row-1,:),f2_new(1:row-1,:),'-r','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 1])
legend(p2,'CO_2', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction [-] ','LineWidth', 2)

subplot(2,2,3);
p3 = plot(x_new(1:row-1,:),f3_new(1:row-1,:),'-k','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 1])
legend(p3,'N_2', 'Location','NorthEast')

subplot(2,2,4);
p4 = plot(x_new(1:row-1,:),f4_new(1:row-1,:),'-c','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 0 1])
legend(p4,'AR', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mass fraction [-] ','LineWidth', 2)
% legend('O_2','CO_2','N_2','Ar','location','NorthEast')
 
% figure(5);
% hold on
% plot(x_new(1:row-1,:),rho_new(1:row-1,:),'-k','LineWidth',2);
% xlabel('Geometric position [m] ','LineWidth', 2)
% ylabel('Density [kg/m^3] ','LineWidth', 2)
% axis([0 1 0 2]);
% legend('\rho','Location','NorthEast')
% set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)

%%
figure(6);
hold on
grid on
subplot(2,2,1);
p1 = plot(x_new(1:row-1,:),D1_new(1:row-1,:),'-b','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 1.4*10^(-5) 2.1*10^(-5)])
legend(p1,'D_{O_2}', 'Location','SouthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity ','LineWidth', 2)

subplot(2,2,2);
p2 = plot(x_new(1:row-1,:),D2_new(1:row-1,:),'-r','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 1.4*10^(-5) 2.1*10^(-5)])
legend(p2,'D_{CO_2}', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity','LineWidth', 2)

subplot(2,2,3);
p3 = plot(x_new(1:row-1,:),D3_new(1:row-1,:),'-k','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 1.4*10^(-5) 2.1*10^(-5)])
legend(p3,'D_{N_2}', 'Location','SouthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity','LineWidth', 2)

subplot(2,2,4);
p4 = plot(x_new(1:row-1,:),D4_new(1:row-1,:),'-c','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
axis([0 1 1.4*10^(-5) 2.1*10^(-5)])
legend(p4,'D_{AR}', 'Location','SouthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Diffusivity ','LineWidth', 2)

%%
figure(7);
hold on
grid on
subplot(2,2,1);
p1 = plot(x_new(2:row-1,:),X1_new(2:row-1,:),'-b','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 1])
legend(p1,'X_{O_2}', 'Location','SouthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mole fraction [-]','LineWidth', 2)

subplot(2,2,2);
p2 = plot(x_new(2:row-1,:),X2_new(2:row-1,:),'-r','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 1])
legend(p2,'X_{CO_2}', 'Location','NorthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mole fraction','LineWidth', 2)

subplot(2,2,3);
p3 = plot(x_new(2:row-1,:),X3_new(2:row-1,:),'-k','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 1])
legend(p3,'X_{N_2}', 'Location','SouthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mole fraction','LineWidth', 2)

subplot(2,2,4);
p4 = plot(x_new(2:row-1,:),X4_new(2:row-1,:),'-c','LineWidth',2);
set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% axis([0 1 0 1])
legend(p4,'X_{AR}', 'Location','SouthEast')
xlabel('Geometric position [m] ','LineWidth', 2)
ylabel('Mole fraction ','LineWidth', 2)

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             here starts the video generation                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fall = struct('O2',f1_new,'CO2',f2_new,'N2',f3_new,'AR',f4_new);
Xall = struct('O2',X1_new,'CO2',X2_new,'N2',X3_new,'AR',X4_new);
Yall = struct('O2',Yr1_new,'CO2',Yr2_new,'N2',Yr3_new,'AR',Yr4_new);

fields = fieldnames(Yall);

for i = 1%:numel(fields)
    fields(i);
    path_Results = strcat('C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_simplified_species\results\Ypermeate',num2str(i),'.avi');
    filename = VideoWriter(path_Results)
    filename.FrameRate = 4;

    axis tight manual % this ensures that getframe() returns a consistent size
    for n = 1:1:col
        
    %% mass fraction
    h1 = figure(10);
    hold on
    p1 = plot(x_new(2:row-1,n),X1_new(2:row-1,n),'-r','LineWidth',2);
    p2 = plot(x_new(2:row-1,n),X2_new(2:row-1,n),'-b','LineWidth',2);
    p3 = plot(x_new(2:row-1,n),X3_new(2:row-1,n),'-k','LineWidth',2);
    p4 = plot(x_new(2:row-1,n),X4_new(2:row-1,n),'-c','LineWidth',2);

    set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15);

    axis([0 1 0 1]);
%     legend([p1 p2], strcat(fields{i}, ' M_r'), strcat(fields{i}, ' M_p'), 'Location','NorthEast');
    xlabel('Geometric position [m] ','LineWidth', 2);
    ylabel('Mass fraction [-] ','LineWidth', 2);
    title(strcat('Time:',num2str(time(row*n))));
    frame1 = getframe(h1);
    open(filename);
    writeVideo(filename,frame1); 
    close(figure(10)); 
    end
    
    close(filename);
    
end

% for i = 1:numel(fields)
%     fields(i);
%     path_Results = strcat('C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Species_sink\results\Y',num2str(i),'.avi');
%     filename = VideoWriter(path_Results)
%     filename.FrameRate = 5;
%     
%     path_Results2 = strcat('C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Species_sink\results\X',num2str(i),'.avi');
%     filename2 = VideoWriter(path_Results2)
%     filename2.FrameRate = 5;
%     
%     h1 = figure(1);
%     h2 = figure(2);
%     axis tight manual % this ensures that getframe() returns a consistent size
%     for n = 1:1:col
%         
%     %% mass fraction
%     h1 = figure(10);
%     plot(x_new(1:row-1,n),Xall.(fields{i})(1:row-1,n),'-b','LineWidth',2);
%     set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% 
%     axis([0 1 0 1]);
%     legend(fields{i}, 'Location','NorthEast')
%     xlabel('Geometric position [m] ','LineWidth', 2)
%     ylabel('Mass fraction [-] ','LineWidth', 2)
%     title(strcat('Time:',num2str(time(row*n))));
%     frame1 = getframe(h1);
%     open(filename);
%     writeVideo(filename,frame1); 
%     
%     %% mole fraction 
%     h2 = figure(20)
%     plot(x_new(1:row-1,n),Xall.(fields{i})(1:row-1,n),'-r','LineWidth',2);
%     set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 15)
% 
%     axis([0 1 0 1]);
%     legend(fields{i}, 'Location','NorthEast')
%     xlabel('Geometric position [m] ','LineWidth', 2)
%     ylabel('Mole fraction [-] ','LineWidth', 2)
%     title(strcat('Time:',num2str(time(row*n))));
%       % Capture the plot as an image 
%         
%     frame2 = getframe(h2);
%     open(filename2);
%     writeVideo(filename2,frame2); 
%       % Write to the avi File 
%           
% %       else 
%      %      end 
%     end
%     close(filename);
%     close(filename2);
% 
% end

% play1 = implay('species1.avi')
% play2 = implay('species2.avi')
% play3 = implay('species3.avi')
% play4 = implay('species4.avi')
