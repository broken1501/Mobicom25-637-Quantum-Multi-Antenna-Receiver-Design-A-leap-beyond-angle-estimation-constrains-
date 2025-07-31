%% basic setting
clear all;clc;

theta=linspace(-pi/2,pi/2,100000);
fs = 5e4;
ts = 1/fs;

freq=[3.0603*10^9,6.099*10^9,9.2243*10^9,11.6223*10^9,14.9313*10^9];
lamda=3*10^8./freq;
d=0.05;
d_lamda=d./lamda;
%% Data Reading

path = '.\data\';
namelist = dir([path,'*.txt']);
l = length(namelist);
P = cell(1,l);
filename=cell(1,l);
data=cell(2,l);
h1_data=cell(1,l);
h2_data=cell(1,l);
antenna_num=l+1;
%% Data alignment

for i = 1:l
    filename{1,i} = [path,namelist(i).name];
    P{1,i} = load(filename{i});
    h1_data{1,i}=P{1,i}(1:length(P{1,i}),1);
    h2_data{1,i}=P{1,i}(1:length(P{1,i}),2);
end
min_length=length(P{1,1});
for i = 2:l
    file_length = length(P{1,i});
    if file_length < min_length
        min_length = file_length;
    end
end
for i=1:l
    h1_data{1,i}=h1_data{1,i}(1:min_length);
    h2_data{1,i}=h2_data{1,i}(1:min_length);   
end


%% Data Preprocessing using IQ down-convertion
h1_data_uni=cell(1,l);
h2_data_uni=cell(1,l);
final_data_h1=cell(1,l);
final_data_h2=cell(1,l);
f_dif = 999.99866; 
t = linspace(0,min_length/fs,min_length);
y = exp(1j*2*pi*f_dif*t);
fil = fir1(500,0.02);
for i=1:l
    h1_data_uni{1,i} = h1_data{1,i} - mean(h1_data{1,i});
    h2_data_uni{1,i} = h2_data{1,i} - mean(h2_data{1,i});
    
    final_data_h1{1,i} = h1_data_uni{1,i}'./y;
    final_data_h2{1,i} = h2_data_uni{1,i}'./y;

    final_data_h1{1,i} = filtfilt(fil,1,final_data_h1{1,i});
    final_data_h2{1,i} = filtfilt(fil,1,final_data_h2{1,i});
end
total_time=t(length(t)-1);

%% Phase Offset Cancellation
phase1=cell(1:l);
phase2=cell(1:l);
phase_diff=cell(1:l);
for i=1:l
    phase1{1,i}=angle(final_data_h1{1,i});
    phase2{1,i}=angle(final_data_h2{1,i});
    phase_diff{1,i}=mod(phase2{1,i}-phase1{1,i}+pi,2*pi)-pi;
end

m_phase1=zeros(1,length(phase_diff{1,1}(1,:)));%Using the first laser as a reference

start_data=zeros(l);
for i=1:l
    start_data(i)=mean(phase_diff{1,i}(2*length(t)/total_time:5*length(t)/total_time));
    phase_diff{1,i}=phase_diff{1,i}-start_data(i);
end

phase_diff{1,1}=(phase_diff{1,1})*(-1);
phase_diff{1,2}=(phase_diff{1,2})*(-1);
phase_diff{1,3}=(phase_diff{1,3})*(-1);
phase_diff{1,4}=phase_diff{1,4}*(-1);
phase_diff{1,5}=(phase_diff{1,5})*(-1)-2*pi;


%% Beamforming with ratio beam
freq=[3.0603*10^9,6.099*10^9,9.2243*10^9,11.6223*10^9,14.9313*10^9];
lamda=3*10^8./freq;
d=0.05;
d_lamda=d./lamda;
for i = 2:antenna_num
    eval(['m_phase',num2str(i),'=','phase_diff{1,i-1}',';']);
end


time=1.1;
theta_truth=0;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');                     
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;                   %sum beam
        plus_sub=w_sub'*a;                   %difference beam
        p(1,j)=plus_add/plus_sub;            %ratio beam
    end
    % 在所有子图代码之前，创建2×4布局，设置间距
    figure('Position', [200 200 2000 600]);  % 适当加大图窗宽度，避免子图拥挤
   
    tl = tiledlayout(2, 4, 'TileSpacing', 'tight','Padding', 'compact');  % 子图间距：紧凑（比'compact'更密）
    
    [max_data,num]=max(p);
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth); 
    r = 1;
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');  
    
    error_0=abs(max_angle);
    title(sprintf('AoA groudtruth：0°,AoA error:%.2f°',abs(max_angle)),'FontSize',14);
 
end
%
time=27.63;
theta_truth=10;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;
        plus_sub=w_sub'*a;
        p(1,j)=plus_add/plus_sub;
    end
    [max_data,num]=max(p);
  
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth); 
    r = 1; 
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    error_1=abs(max_angle-10);
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');  
    title(sprintf('AoA groudtruth：10°,AoA error:%.2f°',abs(max_angle-10)),'FontSize',14);
 
end
%
time=43.5;%
theta_truth=20;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;
        plus_sub=w_sub'*a;
        p(1,j)=plus_add/plus_sub;
    end
    
    [max_data,num]=max(p);
   
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth); 
    r = 1; %
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    error_2=abs(max_angle-20);
    hold on;
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');
    title(sprintf('AoA groudtruth：20°,AoA error:%.2f°',abs(max_angle-20)),'FontSize',14);
 
end
%
time=61.73;
theta_truth=30;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;
        plus_sub=w_sub'*a;
        p(1,j)=plus_add/plus_sub;
    end
   
    [max_data,num]=max(p);
 
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth); 
    r = 1;
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    error_3=abs(max_angle-30);
    hold on;
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');
    title(sprintf('AoA groudtruth：30°,AoA error:%.2f°',abs(max_angle-30)),'FontSize',14);
 
end
%
time=77.92;
theta_truth=40;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;
        plus_sub=w_sub'*a;
        p(1,j)=plus_add/plus_sub;
    end
  
    [max_data,num]=max(p);
   
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth); 
    r = 1; 
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    error_4=abs(max_angle-40);
    hold on;
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');
    title(sprintf('AoA groudtruth：40°,AoA error:%.2f°',abs(max_angle-40)),'FontSize',14);
end
%
time=93.91;
theta_truth=50;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;
        plus_sub=w_sub'*a;
        p(1,j)=plus_add/plus_sub;
    end
  
    [max_data,num]=max(p);
   
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth); 
    r = 1;
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    error_5=abs(max_angle-50);
    hold on;
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');
    title(sprintf('AoA groudtruth：50°,AoA error:%.2f°',abs(max_angle-50)),'FontSize',14);
 
end
%
time=113.511;
theta_truth=60;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;
        plus_sub=w_sub'*a;
        p(1,j)=plus_add/plus_sub;
    end
    
    [max_data,num]=max(p);
    
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth);
    r = 1; 
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    error_6=abs(max_angle-60);
    hold on;
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');
    title(sprintf('AoA groudtruth：60°,AoA error:%.2f°',abs(max_angle-60)),'FontSize',14);
 
end
%
time=125.38;
theta_truth=70;
final_angle=[];
final_angle_p=[];
for k=fix(time*length(t)/total_time)
    
    w_0=[];
    w_1=[];
    w_2=[];
    for i=1:antenna_num
        eval(['w_0','=[','w_0',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=1:round(antenna_num/2)
        eval(['w_1','=[','w_1',',','m_phase',num2str(i),'(1,k)',']',';']);
    end

    for i=round(antenna_num/2)+1:antenna_num
        eval(['w_2','=[','w_2',',','m_phase',num2str(i),'(1,k)',']',';']);
    end
    w_add=exp(1i.*w_0');
    w_sub=[exp(1i.*w_1')',-exp(1i.*w_2')']';
  
    p_add =zeros(1,length(theta));
    p_sub=zeros(1,length(theta),1);
    p=zeros(1,length(theta),1);

    for  j=1:length(theta)
        a_0=exp(1i*2*pi*sin(theta(j)).*d_lamda');
        a=[1,a_0']';
        plus_add=w_add'*a;
        plus_sub=w_sub'*a;
        p(1,j)=plus_add/plus_sub;
    end
    
    [max_data,num]=max(p);
   
    max_angle=rad2deg(theta(1,num));
    final_angle=[final_angle,max_angle];
    normalized_data = (abs(p) - min(abs(p)) )/ ((max(abs(p)) - min(abs(p))));
    nexttile
    polarplot(theta,normalized_data,'LineWidth', 2);
    pax = gca;
    pax.ThetaZeroLocation = 'top';  
    thetalim([-90 90])
    hold on;
    theta_add = deg2rad(theta_truth); 
    r = 1; 
    polarplot(theta_add, r, 'ro', 'MarkerSize', 5, 'LineWidth', 2);
    ax = gca;
    ax.RTick = [];    
    ax.RAxis.Visible = 'off';  
    error_7=abs(max_angle-70);
    hold on;
    legend('Measured AoA','Ground truth','FontSize',14,'NumColumns', 2,'Location', 'southoutside');
    title(sprintf('AoA groudtruth：70°,AoA error:%.2f°',abs(max_angle-70)),'FontSize',14);
 
end

%
vars=[error_0,error_1,error_2,error_3,error_4,error_5,error_6,error_7]
MAE=mean(vars)