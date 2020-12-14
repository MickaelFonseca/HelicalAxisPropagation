clear all
close all
clc
cd 'C:\Users\mcdf\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\Gitlab\Helical Axis\Code'

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

%% 2. Load and calculate input variables theta, k, e1, e3
% cd 'D:\Helical Axis\Test\'
% addpath('D:\Helical Axis\Code')
% Select .c3d files 
% [C3D_filenames, C3D_path, FilterIndex] = uigetfile('.c3d'); %, 'Select .c3d files',['D:\Helical Axis\Data' '/'],'MultiSelect','on');
[C3D_filenames, C3D_path, FilterIndex] = uigetfile('.c3d', 'Select .c3d files',['C:\Users\mcdf\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\Gitlab\Helical Axis\Data' '/'],'MultiSelect','on');

C3D_filenames = cellstr(C3D_filenames)
for file = 1:length(C3D_filenames) 
%     Create_DataStruct(C3D_filenames{file}, C3D_path); % data Montreal
    Create_DataStruct_2(C3D_filenames{file}, C3D_path); % data Geneve

    % Data example for gait
    load Segment.mat    
    Segment(2).T = Q2Tuv_array3(Segment(2).Q); % Foot (i = 2)
    Segment(3).T = Q2Tuv_array3(Segment(3).Q); % Shank (i = 3)
    Segment(4).T = Q2Tuv_array3(Segment(4).Q); % Thigh (i = 4)
    Segment(5).T = Q2Tuv_array3(Segment(5).Q); % Pelvis (i = 5)
        
    % Knee joint
    Joint(3).T = Mprod_array3(Tinv_array3(Segment(4).T),Segment(3).T); % Knee joint (i = 3)
    R = Joint(3).T(1:3,1:3,:);

    % Helical angle and axis
    theta_v = (acos(((R(1,1,:) + R(2,2,:) + R(3,3,:)) ... % Trace of R
        -1)/2)); % Possible singularity when the angle is 0
    k_v = [(R(3,2,:)- R(2,3,:))./(2*sin(theta_v));...
        (R(1,3,:)-R(3,1,:))./(2*sin(theta_v));...
        (R(2,1,:)-R(1,2,:))./(2*sin(theta_v))]; % In thigh SCS
    kx_v=k_v(1,1,:); 
    ky_v=k_v(2,1,:);
    e1x_v=zeros(1,1,size(R,3)); 
    e1y_v=zeros(1,1,size(R,3));   
    e3x_v=R(1,2,:); 
    e3z_v=R(3,2,:);
    
%% 3. Replace input variables in the function of squared standard uncertainty (subs)

    % Define ssu for each imput variable 
    ssu_theta_v=(1); 
    ssu_kx_v =tan(deg2rad(5)); ssu_ky_v =tan(deg2rad(5)); 
    ssu_e1x_v=tan(deg2rad(5)); ssu_e1y_v=tan(deg2rad(5)); 
    ssu_e3x_v=tan(deg2rad(5)); ssu_e3z_v=tan(deg2rad(5));
               
%     Define ssu for theta 1 (return complex numbers, used real)
    ssu_theta1_res(:,file) = interpol(squeeze(double(subs(ssu_theta1, {theta,kx,ky,e1x,e1y,e3x,e3z,...
        ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, {theta_v,kx_v,ky_v,e1x_v,e1y_v,...
        e3x_v,e3z_v,ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}))),101);

%     Define ssu for theta 2
    ssu_theta2_res(:,file) = interpol(squeeze(double(subs(ssu_theta2, {theta,kx,ky,e1x,e1y,e3x,e3z,...
        ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, {theta_v,kx_v,ky_v,e1x_v,e1y_v,....
        e3x_v,e3z_v,ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}))),101);
%     Define ssu for theta 3

    ssu_theta3_res(:,file) = interpol(squeeze(double(subs(ssu_theta3, {theta,kx,ky,e1x,e1y,e3x,e3z,...
        ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, {theta_v,kx_v,ky_v,e1x_v,e1y_v,...
        e3x_v,e3z_v,ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}))),101);
    % 4. Plot theta1 results togheter with a theoretical corridor (mean_theta, std)
    figure(1)
    subplot(3,1,1)
    plot(rad2deg(ssu_theta1_res(:,file)), 'b')%plot ssu_theta1
    hold on

    subplot(3,1,2)
    plot(rad2deg(ssu_theta2_res(:,file)), 'b')%plot ssu_theta2
    hold on

    subplot(3,1,3)
    plot(rad2deg(ssu_theta3_res(:,file)), 'b')%plot ssu_theta3
    hold on
    
%     % Store reference parameters
%     REF.theta(:,file)=interpol(squeeze(theta_v),101);
%     REF.kx(:,file)=interpol(squeeze(kx_v),101);
%     REF.ky(:,file)=interpol(squeeze(ky_v),101);
%     REF.e1x(:,file)=interpol(squeeze(e1x_v),101);
%     REF.e1y(:,file)=interpol(squeeze(e1y_v),101);
%     REF.e3x(:,file)=interpol(squeeze(e3x_v),101);
%     REF.e3z(:,file)=interpol(squeeze(e3z_v),101);
end
subplot(3,1,1)
mean_shaded_std(rad2deg(mean(ssu_theta1_res')), rad2deg(std(ssu_theta1_res')))
xlim([0 100])
title('ssu theta 1')
subplot(3,1,2)
mean_shaded_std(rad2deg(mean(ssu_theta2_res')), rad2deg(std(ssu_theta2_res')))
xlim([0 100])
title('ssu theta 2')
subplot(3,1,3)
mean_shaded_std(rad2deg(mean(ssu_theta3_res')), rad2deg(std(ssu_theta3_res')))
xlim([0 100])
title('ssu theta 3')
suptitle('Montreal PN10')

% save('Parameters_KEV.mat', 'REF') 
