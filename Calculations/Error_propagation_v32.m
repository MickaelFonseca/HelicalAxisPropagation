clear all
close all
clc
% cd 'C:\Users\mcdf\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\Gitlab\Helical Axis\Code'

%% 1. Definition of the symbolic function to calculate uncertainty of theta1, theta2, theta3 (syms)
syms theta kx ky e1x e1y e3x e3z real

syms ssu_theta ssu_kx ssu_ky real
syms ssu_e1x ssu_e1y real
syms ssu_e3x ssu_e3z real

k = [kx, ky, sqrt(1-kx^2-ky^2)]; 
e1 = [e1x, e1y, sqrt(1-e1x^2-e1y^2)];
assume(e1, 'real')
e3 = [e3x, sqrt(1-e3x^2-e3z), e3z];
assume(e3, 'real')
e2 = cross(e3,e1)/sqrt(sum(cross(e3,e1).^2,2));

% Theta1 -----------------------------------------------------------------
ftheta1 = dot(cross(e2,e3),k)*theta;
gtheta1 = dot(cross(e1,e2),e3);
%
dftheta1dtheta = gradient(ftheta1,theta);
dftheta1dkx = gradient(ftheta1,kx);
dftheta1dky = gradient(ftheta1, ky);
dftheta1de1x = gradient(ftheta1, e1x);
dftheta1de1y = gradient(ftheta1, e1y);
dftheta1de3x = gradient(ftheta1,e3x); 
dftheta1de3z = gradient(ftheta1,e3z);

%
dgtheta1dtheta = gradient(gtheta1,theta);
dgtheta1dkx = gradient(gtheta1,kx);
dgtheta1dky = gradient(gtheta1, ky);
dgtheta1de1x = gradient(gtheta1, e1x);
dgtheta1de1y = gradient(gtheta1, e1y);
dgtheta1de3x = gradient(gtheta1, e3x);
dgtheta1de3z = gradient(gtheta1, e3z);

ssu_theta1 = (1/gtheta1)^4*...
        ((gtheta1*dftheta1dtheta - ftheta1*dgtheta1dtheta)^2*ssu_theta + ...
        (gtheta1*dftheta1dkx - ftheta1*dgtheta1dkx)^2*ssu_kx + ...
        (gtheta1*dftheta1dky - ftheta1*dgtheta1dky)^2*ssu_ky + ...
        (gtheta1*dftheta1de1x - ftheta1*dgtheta1de1x)^2*ssu_e1x + ...
        (gtheta1*dftheta1de1y - ftheta1*dgtheta1de1y)^2*ssu_e1y + ...
        (gtheta1*dftheta1de3x - ftheta1*dgtheta1de3x)^2*ssu_e3x + ...
        (gtheta1*dftheta1de3z - ftheta1*dgtheta1de3z)^2*ssu_e3z);
% Theta2 ------------------------------------------------------------------
ftheta2 = dot(cross(e1,e3),k)*theta;
gtheta2 = dot(cross(e1,e2),e3);
%
dftheta2dtheta = gradient(ftheta2,theta);
dftheta2dkx = gradient(ftheta2,kx);
dftheta2dky = gradient(ftheta2, ky);
dftheta2de1x = gradient(ftheta2, e1x);
dftheta2de1y = gradient(ftheta2, e1y);
dftheta2de3x = gradient(ftheta2,e3x); 
dftheta2de3z = gradient(ftheta2,e3z);
%
dgtheta2dtheta = gradient(gtheta2,theta);
dgtheta2dkx = gradient(gtheta2,kx);
dgtheta2dky = gradient(gtheta2, ky);
dgtheta2de1x = gradient(gtheta2, e1x);
dgtheta2de1y = gradient(gtheta2, e1y);
dgtheta2de3x = gradient(gtheta2, e3x);
dgtheta2de3z = gradient(gtheta2, e3z);

ssu_theta2 = (1/gtheta2)^4*...
        ((gtheta2*dftheta2dtheta - ftheta2*dgtheta2dtheta)^2*ssu_theta + ...
        (gtheta2*dftheta2dkx - ftheta2*dgtheta2dkx)^2*ssu_kx + ...
        (gtheta2*dftheta2dky - ftheta2*dgtheta2dky)^2*ssu_ky + ...
        (gtheta2*dftheta2de1x - ftheta2*dgtheta2de1x)^2*ssu_e1x + ...
        (gtheta2*dftheta2de1y - ftheta2*dgtheta2de1y)^2*ssu_e1y + ...
        (gtheta2*dftheta2de3x - ftheta2*dgtheta2de3x)^2*ssu_e3x + ...
        (gtheta2*dftheta2de3z - ftheta2*dgtheta2de3z)^2*ssu_e3z);

% Theta3 ------------------------------------------------------------------
ftheta3 = dot(cross(e1,e2),k)*theta;
gtheta3 = dot(cross(e1,e2),e3);
%
dftheta3dtheta = gradient(ftheta3,theta);
dftheta3dkx = gradient(ftheta3,kx);
dftheta3dky = gradient(ftheta3, ky);
dftheta3de1x = gradient(ftheta3, e1x);
dftheta3de1y = gradient(ftheta3, e1y);
dftheta3de3x = gradient(ftheta3,e3x); 
dftheta3de3z = gradient(ftheta3,e3z);
%
dgtheta3dtheta = gradient(gtheta3,theta);
dgtheta3dkx = gradient(gtheta3,kx);
dgtheta3dky = gradient(gtheta3, ky);
dgtheta3de1x = gradient(gtheta3, e1x);
dgtheta3de1y = gradient(gtheta3, e1y);
dgtheta3de3x = gradient(gtheta3, e3x);
dgtheta3de3z = gradient(gtheta3, e3z);

ssu_theta3 = (1/gtheta3)^4*...
        ((gtheta3*dftheta3dtheta - ftheta3*dgtheta3dtheta)^2*ssu_theta + ...
        (gtheta3*dftheta3dkx - ftheta3*dgtheta3dkx)^2*ssu_kx + ...
        (gtheta3*dftheta3dky - ftheta3*dgtheta3dky)^2*ssu_ky + ...
        (gtheta3*dftheta3de1x - ftheta3*dgtheta3de1x)^2*ssu_e1x + ...
        (gtheta3*dftheta3de1y - ftheta3*dgtheta3de1y)^2*ssu_e1y + ...
        (gtheta3*dftheta3de3x - ftheta3*dgtheta3de3x)^2*ssu_e3x + ...
        (gtheta3*dftheta3de3z - ftheta3*dgtheta3de3z)^2*ssu_e3z);

% %% 2. Load and calculate input variables theta, k, e1, e3
% % cd 'D:\Helical Axis\Test\'
% % addpath('D:\Helical Axis\Code')
% % Select .c3d files 
% % [C3D_filenames, C3D_path, FilterIndex] = uigetfile('.c3d'); %, 'Select .c3d files',['D:\Helical Axis\Data' '/'],'MultiSelect','on');
% [C3D_filenames, C3D_path, FilterIndex] = uigetfile('.c3d', 'Select .c3d files',['C:\Users\mcdf\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\Gitlab\Helical Axis\Data' '/'],'MultiSelect','on');
% 
% C3D_filenames = cellstr(C3D_filenames)
% for file = 1:length(C3D_filenames) 
% %     Create_DataStruct(C3D_filenames{file}, C3D_path); % data Montreal
%     Create_DataStruct_2(C3D_filenames{file}, C3D_path); % data Geneve
% 
%     % Data example for gait
%     load Segment.mat    
%     Segment(2).T = Q2Tuv_array3(Segment(2).Q); % Foot (i = 2)
%     Segment(3).T = Q2Tuv_array3(Segment(3).Q); % Shank (i = 3)
%     Segment(4).T = Q2Tuv_array3(Segment(4).Q); % Thigh (i = 4)
%     Segment(5).T = Q2Tuv_array3(Segment(5).Q); % Pelvis (i = 5)
%         
%     % Knee joint
%     Joint(3).T = Mprod_array3(Tinv_array3(Segment(4).T),Segment(3).T); % Knee joint (i = 3)
%     R = Joint(3).T(1:3,1:3,:);
% 
%     % Helical angle and axis
%     theta_v = (acos(((R(1,1,:) + R(2,2,:) + R(3,3,:)) ... % Trace of R
%         -1)/2)); % Possible singularity when the angle is 0
%     k_v = [(R(3,2,:)- R(2,3,:))./(2*sin(theta_v));...
%         (R(1,3,:)-R(3,1,:))./(2*sin(theta_v));...
%         (R(2,1,:)-R(1,2,:))./(2*sin(theta_v))]; % In thigh SCS
%     kx_v=k_v(1,1,:); 
%     ky_v=k_v(2,1,:);
%     e1x_v=zeros(1,1,size(R,3)); 
%     e1y_v=zeros(1,1,size(R,3));   
%     e3x_v=R(1,2,:); 
%     e3z_v=R(3,2,:);

%% Load reference data

load('Parameters_KEV.mat');
theta_v = REF.theta;
kx_v = REF.kx;
ky_v = REF.ky;
e1x_v = REF.e1x;
e1y_v = REF.e1y;
e3x_v = REF.e3x;
e3z_v = REF.e3z;

% Display
disp(['Mean experimental SD on theta: ',num2str(mean(rad2deg(std(theta_v,[],2))))]); % Experimental standard deviation
disp(['Mean experimental SD on kx: ',num2str(mean(rad2deg(std(atan(kx_v),[],2))))]); % Experimental standard deviation
disp(['Mean experimental SD on ky: ',num2str(mean(rad2deg(std(atan(ky_v),[],2))))]); % Experimental standard deviation
disp(['Mean experimental SD on e3x: ',num2str(mean(rad2deg(std(atan(e3x_v),[],2))))]); % Experimental standard deviation
disp(['Mean experimental SD on e3z: ',num2str(mean(rad2deg(std(atan(e3z_v),[],2))))]); % Experimental standard deviation


for t = 1:size(theta_v,2)
    
    theta1_res(:,t) = double(subs(ftheta1/gtheta1,...
        {theta,kx,ky,e1x,e1y,e3x,e3z,}, ...
        {theta_v(:,t),kx_v(:,t),ky_v(:,t),e1x_v(:,t),e1y_v(:,t),e3x_v(:,t),e3z_v(:,t)}));
    
    theta2_res(:,t) = double(subs(ftheta2/gtheta2,...
        {theta,kx,ky,e1x,e1y,e3x,e3z,}, ...
        {theta_v(:,t),kx_v(:,t),ky_v(:,t),e1x_v(:,t),e1y_v(:,t),e3x_v(:,t),e3z_v(:,t)}));
    
    theta3_res(:,t) = double(subs(ftheta3/gtheta3,...
        {theta,kx,ky,e1x,e1y,e3x,e3z,}, ...
        {theta_v(:,t),kx_v(:,t),ky_v(:,t),e1x_v(:,t),e1y_v(:,t),e3x_v(:,t),e3z_v(:,t)}));
    
end
   
%% 3. Replace input variables in the function of squared standard uncertainty (subs)

% Uncertainty
% Squared standard uncertainty ((standard deviation)^2)
ui = deg2rad(0)^2; % Intrinsic 
ue = deg2rad(0)^2; % Extrinsic

utheta = deg2rad(0)^2;
uk = deg2rad(0)^2;
ue1 = deg2rad(0)^2;
ue3 = deg2rad(5)^2;

% % Define ssu for each imput variable
% % Intrinsic
% ssu_theta_v = ui;
% ssu_kx_v = tan(ui); ssu_ky_v = tan(ui);
% % Extrinsic
% ssu_e1x_v = tan(ue); ssu_e1y_v = tan(ue);
% ssu_e3x_v =tan(ue); ssu_e3z_v = tan(ue);

% Define ssu for each imput variable
% Intrinsic
ssu_theta_v = utheta;
ssu_kx_v = tan(uk); ssu_ky_v = tan(uk);
% Extrinsic
ssu_e1x_v = tan(ue1); ssu_e1y_v = tan(ue1);
ssu_e3x_v =tan(ue3); ssu_e3z_v = tan(ue3);

% Define ssu for theta 1
ssu_theta1_res = double(subs(ssu_theta1, ...
    {theta,kx,ky,e1x,e1y,e3x,e3z,...
    ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
    {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
    ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));

% Define ssu for theta 2
ssu_theta2_res = double(subs(ssu_theta2, ...
    {theta,kx,ky,e1x,e1y,e3x,e3z,...
    ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
    {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
    ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));

% Define ssu for theta 3
ssu_theta3_res = double(subs(ssu_theta3, ...
    {theta,kx,ky,e1x,e1y,e3x,e3z,...
    ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
    {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
    ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));


%% 4. Plot theta1 results togheter with a theoretical corridor (mean_theta, std)

h = figure(1)
%
sp1 = subplot(3,1,1)
A = rad2deg(mean(theta1_res,2));
B = rad2deg(std(theta1_res,[],2));
Max = A+B; Min = A-B;
x = [1:101 101:-1:1];
y = [Max;Min(end:-1:1)];
hold on
plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
xlim([1 101]);
B = rad2deg(sqrt(ssu_theta1_res)); % Standard uncertainty
Max = A+B; Min = A-B;
x = [1:101 101:-1:1];
y = [Max;Min(end:-1:1)];
fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
ylabel('Theta 1')
title('SSU e3 = 5')
sp1.XAxis.TickValues = [0 20 40 60 80 100];
sp1.FontSize = 20

sp2 = subplot(3,1,2)
A = rad2deg(mean(theta2_res,2));
B = rad2deg(std(theta2_res,[],2));
Max = A+B; Min = A-B;
x = [1:101 101:-1:1];
y = [Max;Min(end:-1:1)];
hold on
plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
xlim([1 101]);
B = rad2deg(sqrt(ssu_theta2_res)); % Standard uncertainty
Max = A+B; Min = A-B;
x = [1:101 101:-1:1];
y = [Max;Min(end:-1:1)];
fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
ylabel('Theta 2')
sp2.XAxis.TickValues = [0 20 40 60 80 100];
sp2.FontSize = 20

sp3 = subplot(3,1,3)
A = rad2deg(mean(theta3_res,2));
B = rad2deg(std(theta3_res,[],2));
Max = A+B; Min = A-B;
x = [1:101 101:-1:1];
y = [Max;Min(end:-1:1)];
hold on
plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
xlim([1 101]);
B = rad2deg(sqrt(ssu_theta3_res)); % Standard uncertainty
Max = A+B; Min = A-B;
x = [1:101 101:-1:1];
y = [Max;Min(end:-1:1)];
fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
ylabel('Theta 3')
xlabel('% of Gait Cycle')
legend({'Experimental (mean)','Experimental (SD)','Theoritical (SU)'})
sp3.XAxis.TickValues = [0 20 40 60 80 100];
sp3.FontSize = 20

h = gcf;
set(h, 'PaperPositionMode', 'auto');
set(h, 'PaperOrientation', 'landscape');
set(h, 'Position', [50 50 1200 800]);
print(gcf, '-dpdf', 'SSUe35.pdf')

%% 5. Sensitivity to parameters
SSUtheta = [0,1,3,5];
SSUk     = [0,1,3,5];
SSUe1    = [0,1,3,5];
SSUe3    = [0,1,3,5];

% calculate experimental SD
theta1sd = rad2deg(std(theta1_res,[],2))';
theta2sd = rad2deg(std(theta2_res,[],2))';
theta3sd = rad2deg(std(theta3_res,[],2))';

RMSD_t1_cc=zeros(length(SSUtheta)*length(SSUk)*length(SSUe1)*length(SSUe3),1);
RMSD_t1_st=zeros(length(RMSD_t1_cc),1);
RMSD_t1_sw=zeros(length(RMSD_t1_cc),1);
RMSD_t2_cc=zeros(length(RMSD_t1_cc),1);
RMSD_t2_st=zeros(length(RMSD_t1_cc),1);
RMSD_t2_sw=zeros(length(RMSD_t1_cc),1);
RMSD_t3_cc=zeros(length(RMSD_t1_cc),1);
RMSD_t3_st=zeros(length(RMSD_t1_cc),1);
RMSD_t3_sw=zeros(length(RMSD_t1_cc),1);
cssut = zeros(length(RMSD_t1_cc),1);
cssuk = zeros(length(RMSD_t1_cc),1);
cssue1 = zeros(length(RMSD_t1_cc),1);
cssue3 = zeros(length(RMSD_t1_cc),1);

count = 1;
for ssut = 1:length(SSUtheta)
    for ssuk = 1:length(SSUk)
        for ssue1 = 1:length(SSUe1)
            for ssue3 = 1:length(SSUe3)
                tic;
                utheta = deg2rad(SSUtheta(ssut))^2;
                uk = deg2rad(SSUk(ssuk))^2;
                ue1 = deg2rad(SSUe1(ssue1))^2;
                ue3 = deg2rad(SSUe3(ssue3))^2;
                
                % Define ssu for each imput variable
                % Intrinsic
                ssu_theta_v = utheta;
                ssu_kx_v = tan(uk); ssu_ky_v = tan(uk);
                % Extrinsic
                ssu_e1x_v = tan(ue1); ssu_e1y_v = tan(ue1);
                ssu_e3x_v =tan(ue3); ssu_e3z_v = tan(ue3);

                % Define ssu for theta 1
                ssu_theta1_res = double(subs(ssu_theta1, ...
                    {theta,kx,ky,e1x,e1y,e3x,e3z,...
                    ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
                    {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
                    ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));

                % Define ssu for theta 2
                ssu_theta2_res = double(subs(ssu_theta2, ...
                    {theta,kx,ky,e1x,e1y,e3x,e3z,...
                    ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
                    {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
                    ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));

                % Define ssu for theta 3
                ssu_theta3_res = double(subs(ssu_theta3, ...
                    {theta,kx,ky,e1x,e1y,e3x,e3z,...
                    ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
                    {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
                    ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));
                
                % calculate SD theoretical corridor
                ssutheta1 = rad2deg(sqrt(ssu_theta1_res))';
                ssutheta2 = rad2deg(sqrt(ssu_theta2_res))';
                ssutheta3 = rad2deg(sqrt(ssu_theta3_res))';
                
                RMSD_t1_cc(count) = sqrt(mean(theta1sd - ssutheta1).^2);
                RMSD_t1_st(count) = sqrt(mean(theta1sd(1:60) - ssutheta1(1:60)).^2);
                RMSD_t1_sw(count) = sqrt(mean(theta1sd(61:end) - ssutheta1(61:end)).^2);
                
                RMSD_t2_cc(count) = sqrt(mean(theta2sd - ssutheta2).^2);
                RMSD_t2_st(count) = sqrt(mean(theta2sd(1:60) - ssutheta2(1:60)).^2);
                RMSD_t2_sw(count) = sqrt(mean(theta2sd(61:end) - ssutheta2(61:end)).^2);
                
                RMSD_t3_cc(count) = sqrt(mean(theta3sd - ssutheta3).^2);
                RMSD_t3_st(count) = sqrt(mean(theta3sd(1:60) - ssutheta3(1:60)).^2);
                RMSD_t3_sw(count) = sqrt(mean(theta3sd(61:end) - ssutheta3(61:end)).^2);
                
                cssut(count)  = SSUtheta(ssut);
                cssuk(count)  = SSUtheta(ssuk);
                cssue1(count) = SSUtheta(ssue1);
                cssue3(count) = SSUtheta(ssue3);
                
                disp(count)
                count = count+1;
                toc
            end
        end
    end
end

filename = 'SSU_sensitivity.mat';
save(filename)
load(filename)

% create table
T = table(cssuk, cssut, cssue1, cssue3, abs(RMSD_t1_cc), RMSD_t1_st, RMSD_t1_sw, RMSD_t2_cc, RMSD_t2_st, RMSD_t2_sw, RMSD_t3_cc, RMSD_t3_st, RMSD_t3_sw)
T.Properties.VariableNames = {'cssuk', 'cssut', 'cssue1', 'cssue3', 'RMSD_t1_cc', 'RMSD_t1_st', 'RMSD_t1_sw', 'RMSD_t2_cc', 'RMSD_t2_st', 'RMSD_t2_sw', 'RMSD_t3_cc', 'RMSD_t3_st', 'RMSD_t3_sw'}