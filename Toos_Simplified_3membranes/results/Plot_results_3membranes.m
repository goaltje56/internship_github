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
iX_k1  = find(strcmpi('X_k1',A.colheaders));
iX_k2  = find(strcmpi('X_k2',A.colheaders));
iX2_k1  = find(strcmpi('X2_k1',A.colheaders));
iX2_k2  = find(strcmpi('X2_k2',A.colheaders));
istagecut = find(strcmpi('stagecut1',A.colheaders));

% irho  = find(strcmpi('density',A.colheaders));
% 
% iD1  = find(strcmpi('Dspecies1',B.colheaders));
% iD2  = find(strcmpi('Dspecies2',B.colheaders));
% iD3  = find(strcmpi('Dspecies3',B.colheaders));
% iD4  = find(strcmpi('Dspecies4',B.colheaders));
ix2   = find(strcmpi('Position2',B.colheaders));
iX3_k1  = find(strcmpi('X3_k1',B.colheaders));
iX3_k2  = find(strcmpi('X3_k2',B.colheaders));
iX4_k1  = find(strcmpi('X4_k1',B.colheaders));
iX4_k2  = find(strcmpi('X4_k2',B.colheaders));
istagecut2 = find(strcmpi('stagecut2',B.colheaders));

% iMr  = find(strcmpi('Mr',C.colheaders));
% iMp  = find(strcmpi('Mp',C.colheaders));
ix3   = find(strcmpi('Position3',C.colheaders));
iX5_k1 = find(strcmpi('X5_k1',C.colheaders));
iX5_k2 = find(strcmpi('X5_k2',C.colheaders));
iX6_k1 = find(strcmpi('X6_k1',C.colheaders));
iX6_k2 = find(strcmpi('X6_k2',C.colheaders));
istagecut3 = find(strcmpi('stagecut3',C.colheaders));

% store the values in the appropriate array
time = y(:,it);
x = y(:,ix);
% x_u = y(:,ix_u);
% u = y(:,iu);
% T = y(:,iT);
% P = y(:,iP);
X_k1 = y(:,iX_k1);
X_k2 = y(:,iX_k2);
X2_k1 = y(:,iX2_k1);
X2_k2 = y(:,iX2_k2);
% rho = y(:,irho);

x2 = yy(:,ix2);
X3_k1 = yy(:,iX3_k1);
X3_k2 = yy(:,iX3_k2);
X4_k1 = yy(:,iX4_k1);
X4_k2 = yy(:,iX4_k2);
% X1 = yy(:,iX1);
% X2 = yy(:,iX2);
% X3 = yy(:,iX3);
% X4 = yy(:,iX4);
x3 = yyy(:,ix3);
X5_k1 = yyy(:,iX5_k1);
X5_k2 = yyy(:,iX5_k2);
X6_k1 = yyy(:,iX6_k1);
X6_k2 = yyy(:,iX6_k2);

stagecut = y(:,istagecut);
stagecut2 = yy(:,istagecut2);
stagecut3 = yyy(:,istagecut3);
% reshape the array such that they can be plotted within
% the computational domain. 
col = length(unique(time));
row = length(X_k1)/col;
row2 = length(X3_k1)/col;
row3 = length(X5_k1)/col;

time_new = reshape(time, [row, col]);
x_new = reshape(x,[row, col]);
% x_u_new = reshape(x_u,[row, col]);
% u_new = reshape(u,[row, col]);
% T_new = reshape(T,[row, col]);
% P_new = reshape(P,[row, col]);
X_k1_new = reshape(X_k1,[row, col]);
X_k2_new = reshape(X_k2,[row, col]);
X2_k1_new = reshape(X2_k1,[row, col]);
X2_k2_new = reshape(X2_k2,[row, col]);
% rho_new = reshape(rho,[row, col]);

% D1_new = reshape(D1,[row, col]);
% D2_new = reshape(D2,[row, col]);
% D3_new = reshape(D3,[row, col]);
% D4_new = reshape(D4,[row, col]);
x2_new = reshape(x2,[row2, col]);
X3_k1_new = reshape(X3_k1,[row2, col]);
X3_k2_new = reshape(X3_k2,[row2, col]);
X4_k1_new = reshape(X4_k1,[row2, col]);
X4_k2_new = reshape(X4_k2,[row2, col]);

% Mr_new  = reshape(Mr,[row2, col]);
% Mp_new  = reshape(Mp,[row2, col]);
x3_new = reshape(x3,[row3, col]);
X5_k1_new = reshape(X5_k1,[row3, col]);
X5_k2_new = reshape(X5_k2,[row3, col]);
X6_k1_new = reshape(X6_k1,[row3, col]);
X6_k2_new = reshape(X6_k2,[row3, col]);

stagecut_new = reshape(stagecut,[row, col]);
stagecut2_new = reshape(stagecut2,[row2, col]);
stagecut3_new = reshape(stagecut3,[row3, col]);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             here starts the video generation                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fall = struct('CO2',X_k1_new,'AR',X_k2_new);
% Xall = struct('O2',X1_new,'CO2',X2_new,'N2',X3_new,'AR',X4_new);
% Yall = struct('O2',Yr1_new,'CO2',Yr2_new,'N2',Yr3_new,'AR',Yr4_new);

fields = fieldnames(Fall);

for i = 3%:numel(fields)
%     fields(i);
    path_Results = strcat('C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_Simplified_3membranes\results\3membranes_transient',num2str(i),'.avi');
    filename = VideoWriter(path_Results)
    filename.FrameRate = 3;

    axis tight manual % this ensures that getframe() returns a consistent size
    for n = 1:1:col
        
    %% mass fraction
    h1 = figure('visible','off');

    subplot(2,2,1)
    hold on
    p1 = plot(time_new(1,1:n), X_k1_new(1,1:n),'-r','LineWidth',1);
    p2 = plot(time_new(1,1:n), X_k2_new(1,1:n),'-b','LineWidth',1);
%     p3 = plot(time_new(1,1:n), f3_new(1,1:n),'-k','LineWidth',1);
%     p4 = plot(time_new(1,1:n), f4_new(1,1:n),'-c','LineWidth',1);
    title(strcat('Time:',num2str(time(row*n))));
    xlabel('time [s]','LineWidth', 1);
    ylabel('X_{feed} [-] ','LineWidth', 1);
    axis([0 200 0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0','0.2','0.4','0.6','0.8','1'})      
    xticks([0 25 50 75 100 125 150 175 200])
    xticklabels({'0','25','50','75','100','125','150','175','200'})    
%     yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
%     yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})      
%     xticks([0 20 40 60 80 100 120 140 160 180 200])
%     xticklabels({'0','20','40','60','80','100','120','140','160','180','200'})
    set(gca, 'box', 'on', 'LineWidth', 1, 'FontSize', 8);

    subplot(2,2,4)
    hold on
    p1 = plot(stagecut_new(2:row-1,n),X2_k1_new(2:row-1,n),'-r','LineWidth',1);
    p2 = plot(stagecut_new(2:row-1,n),X2_k2_new(2:row-1,n),'-b','LineWidth',1);
%     p3 = plot(stagecut_new(2:row-1,n),X3_new(2:row-1,n),'-k','LineWidth',1);
%     p4 = plot(stagecut_new(2:row-1,n),X4_new(2:row-1,n),'-c','LineWidth',1);
    xlabel('stage cut \Theta [-]','LineWidth', 1);
    ylabel('X [-] ','LineWidth', 1);
    title('Permeate');
    axis([0 0.7 0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0','0.2','0.4','0.6','0.8','1'})      
%     xticks([0 0.05 0.1 0.15 0.2])
%     xticklabels({'0','0.5','0.1','0.15','0.2'})    
    set(gca, 'box', 'on', 'LineWidth', 1, 'FontSize', 8);

    ax1 = subplot(2,2,3)
    hold on
    p1 = plot(1,1,'-r','LineWidth',1);
    p2 = plot(1,1,'-b','LineWidth',1);
%     p3 = plot(1,1,'-k','LineWidth',1);
%     p4 = plot(1,1,'-c','LineWidth',1);
    set(gca, 'box', 'on', 'LineWidth', 2, 'FontSize', 8);
    legend(ax1,[p1;p2],'CO_2','AR','Location','best') 
    legend('boxoff')
    axis([0 1 0 1]);
    xticks([0 1])
    xticklabels({'', ''});
    yticks([0 1])
    yticklabels({'', ''});
    
    subplot(2,2,2)
    hold on
    p1 = plot(stagecut_new(2:row-1,n),X_k1_new(2:row-1,n),'-r','LineWidth',1);
    p2 = plot(stagecut_new(2:row-1,n),X_k2_new(2:row-1,n),'-b','LineWidth',1);
%     p3 = plot(stagecut_new(2:row-1,n),f3_new(2:row-1,n),'-k','LineWidth',1);
%     p4 = plot(stagecut_new(2:row-1,n),f4_new(2:row-1,n),'-c','LineWidth',1);
    axis([0 0.7 0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0','0.2','0.4','0.6','0.8','1'}) 
%     xticks([0 0.05 0.1 0.15 0.7])
%     xticklabels({'0','0.5','0.1','0.15','0.2'})      
    title('Retentate')
    xlabel('stage cut \Theta [-]','LineWidth', 1);
    ylabel('X [-] ','LineWidth', 1);
    set(gca, 'box', 'on', 'LineWidth', 1, 'FontSize', 8);
%     legend([p1;p2;p3;p4],'O_2','CO_2','H_2O','AR','Location','best') 


%     legend([p1 p2], strcat(fields{i}, ' M_r'), strcat(fields{i}, ' M_p'), 'Location','NorthEast');

    frame1 = getframe(h1);
    open(filename);
    writeVideo(filename,frame1); 
    close(h1); 
    end
    
    close(filename);
    
end

for i = 4%:numel(fields)
%     fields(i);
    path_Results = strcat('C:\Users\s137280\Documents\Master_tue\Internship\internship_github\Toos_Simplified_3membranes\results\3membranes_transient',num2str(i),'.avi');
    filename = VideoWriter(path_Results)
    filename.FrameRate = 3;

    axis tight manual % this ensures that getframe() returns a consistent size
    for n = 1:1:col
        
    %% mass fraction
    h1 = figure('visible','off');
    
    subplot(2,2,1)
    hold on
    p1 = plot(stagecut2_new(2:row2-1,n), X3_k1_new(2:row2-1,n),'-r','LineWidth',1);
    p2 = plot(stagecut2_new(2:row2-1,n), X3_k2_new(2:row2-1,n),'-b','LineWidth',1);
%     p3 = plot(time_new(1,1:n), f3_new(1,1:n),'-k','LineWidth',1);
%     p4 = plot(time_new(1,1:n), f4_new(1,1:n),'-c','LineWidth',1);
    title(strcat('Retenate Membrane 2, Time:',num2str(time(row*n))));
    xlabel('\Theta [-]','LineWidth', 1);
    ylabel('X [-] ','LineWidth', 1);
    axis([0 0.21 0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0','0.2','0.4','0.6','0.8','1'})      
    xticks([0 0.05 0.1 0.15 0.2])
    xticklabels({'0','0.5','0.1','0.15','0.2'})      
%     yticks([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
%     yticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})      
%     xticks([0 20 40 60 80 100 120 140 160 180 200])
%     xticklabels({'0','20','40','60','80','100','120','140','160','180','200'})
    set(gca, 'box', 'on', 'LineWidth', 1, 'FontSize', 8);

    subplot(2,2,4)
    hold on
    p1 = plot(stagecut3_new(2:row3-1,n),X6_k1_new(2:row3-1,n),'-r','LineWidth',1);
    p2 = plot(stagecut3_new(2:row3-1,n),X6_k2_new(2:row3-1,n),'-b','LineWidth',1);
%     p3 = plot(stagecut_new(2:row-1,n),X3_new(2:row-1,n),'-k','LineWidth',1);
%     p4 = plot(stagecut_new(2:row-1,n),X4_new(2:row-1,n),'-c','LineWidth',1);
    xlabel('stage cut \Theta [-]','LineWidth', 1);
    ylabel('X [-] ','LineWidth', 1);
    title('Permeate Membrane 3');
    axis([0 0.8 0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0','0.2','0.4','0.6','0.8','1'})      
%     xticks([0 0.05 0.1 0.15 0.2])
%     xticklabels({'0','0.5','0.1','0.15','0.2'})    
    set(gca, 'box', 'on', 'LineWidth', 1, 'FontSize', 8);

    ax1 = subplot(2,2,3)
    hold on
   p1 = plot(stagecut2_new(2:row2-1,n),X4_k1_new(2:row2-1,n),'-r','LineWidth',1);
   p2 = plot(stagecut2_new(2:row2-1,n),X4_k2_new(2:row2-1,n),'-b','LineWidth',1);
%     p3 = plot(1,1,'-k','LineWidth',1);
%     p4 = plot(1,1,'-c','LineWidth',1);
    xlabel('stage cut \Theta [-]','LineWidth', 1);
    ylabel('X [-] ','LineWidth', 1);
    title('Permeate Membrane 2');
    axis([0 0.21 0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0','0.2','0.4','0.6','0.8','1'})      
    xticks([0 0.05 0.1 0.15 0.2])
    xticklabels({'0','0.5','0.1','0.15','0.2'})    
    set(gca, 'box', 'on', 'LineWidth', 1, 'FontSize', 8);
    
    subplot(2,2,2)
    hold on
    p1 = plot(stagecut3_new(2:row3-1,n),X5_k1_new(2:row3-1,n),'-r','LineWidth',1);
    p2 = plot(stagecut3_new(2:row3-1,n),X5_k2_new(2:row3-1,n),'-b','LineWidth',1);
%     p3 = plot(stagecut_new(2:row-1,n),f3_new(2:row-1,n),'-k','LineWidth',1);
%     p4 = plot(stagecut_new(2:row-1,n),f4_new(2:row-1,n),'-c','LineWidth',1);
    axis([0 0.8 0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    yticklabels({'0','0.2','0.4','0.6','0.8','1'}) 
%     xticks([0 0.05 0.1 0.15 0.7])
%     xticklabels({'0','0.5','0.1','0.15','0.2'})      
    title('Retentate Membrane 3')
    xlabel('stage cut \Theta [-]','LineWidth', 1);
    ylabel('X [-] ','LineWidth', 1);
    set(gca, 'box', 'on', 'LineWidth', 1, 'FontSize', 8);
%     legend([p1;p2;p3;p4],'O_2','CO_2','H_2O','AR','Location','best') 


%     legend([p1 p2], strcat(fields{i}, ' M_r'), strcat(fields{i}, ' M_p'), 'Location','NorthEast');

    frame1 = getframe(h1);
    open(filename);
    writeVideo(filename,frame1); 
    close(h1); 
    end
    
    close(filename);
    
end