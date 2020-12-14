clear all
close all
clc


%% 1. Definition of the symbolic function to calculate uncertainty of theta1, theta2, theta3 (syms)
syms theta kx ky kz e1x e1y e1z e2x e2y e2z e3x e3y e3z real

k = [kx, ky, kz];
e1 = [e1x, e1y, e1z];
% e2 = [e2x, e2y, e2z];
e3 = [e3x, e3y, e3z];
e2 = cross(e3,e1)/sqrt(sum(cross(e3,e1).^2,2));

syms ssu_theta ssu_kx ssu_ky ssu_kz ssu_e1x ssu_e1y ssu_e1z real
syms ssu_e2x ssu_e2y ssu_e2z ssu_e3x ssu_e3y ssu_e3z real

ssu_k = [ssu_kx, ssu_ky, ssu_kz];
ssu_e1 = [ssu_e1x, ssu_e1y, ssu_e1z];
ssu_e2 = [ssu_e2x, ssu_e2y, ssu_e2z];
ssu_e3 = [ssu_e3x, ssu_e3y, ssu_e3z];

% Theta1 -----------------------------------------------------------------
ftheta1 = dot(cross(e2,e3),k)*theta;
gtheta1 = dot(cross(e1,e2),e3);
%
dftheta1dtheta = gradient(ftheta1,theta)
dftheta1dk = gradient(ftheta1,k); % 0
dftheta1de1 = gradient(ftheta1,e1); 
dftheta1de2 = gradient(ftheta1,e2);
dftheta1de3 = gradient(ftheta1,e3);
%
dgtheta1dtheta = gradient(gtheta1,theta); % 0
dgtheta1dk = gradient(gtheta1,k); % 0
dgtheta1de1 = gradient(gtheta1,e1); 
dgtheta1de2 = gradient(gtheta1,e2);
dgtheta1de3 = gradient(gtheta1,e3);

ssu_theta1 = (1/gtheta1)^4*...
        ((gtheta1*dftheta1dtheta - ftheta1*dgtheta1dtheta)^2*ssu_theta + ...
        (gtheta1*dftheta1dk(1) - ftheta1*dgtheta1dk(1))^2*ssu_k(1) + ...
        (gtheta1*dftheta1dk(2) - ftheta1*dgtheta1dk(2))^2*ssu_k(2) + ...
        (gtheta1*dftheta1dk(3) - ftheta1*dgtheta1dk(3))^2*ssu_k(3) + ...
        (gtheta1*dftheta1de1(1) - ftheta1*dgtheta1de1(1))^2*ssu_e1(1) + ...
        (gtheta1*dftheta1de1(2) - ftheta1*dgtheta1de1(1))^2*ssu_e1(2) + ...
        (gtheta1*dftheta1de1(3) - ftheta1*dgtheta1de1(1))^2*ssu_e1(3) + ...
        (gtheta1*dftheta1de2(1) - ftheta1*dgtheta1de2(1))^2*ssu_e2(1) + ...
        (gtheta1*dftheta1de2(2) - ftheta1*dgtheta1de2(2))^2*ssu_e2(2) + ...
        (gtheta1*dftheta1de2(3) - ftheta1*dgtheta1de2(3))^2*ssu_e2(3) + ...
        (gtheta1*dftheta1de3(1) - ftheta1*dgtheta1de3(1))^2*ssu_e3(1) + ...
        (gtheta1*dftheta1de3(2) - ftheta1*dgtheta1de3(2))^2*ssu_e3(2) + ...
        (gtheta1*dftheta1de3(3) - ftheta1*dgtheta1de3(1))^2*ssu_e3(3));

% Theta2 ------------------------------------------------------------------
ftheta2 = dot(cross(e1,e3),k)*theta;
gtheta2 = dot(cross(e1,e2),e3);
%
dftheta2dtheta = gradient(ftheta2,theta)
dftheta2dk = gradient(ftheta2,k); % 0
dftheta2de1 = gradient(ftheta2,e1); 
dftheta2de2 = gradient(ftheta2,e2);
dftheta2de3 = gradient(ftheta2,e3);
%
dgtheta2dtheta = gradient(gtheta2,theta); % 0
dgtheta2dk = gradient(gtheta2,k); % 0
dgtheta2de1 = gradient(gtheta2,e1); 
dgtheta2de2 = gradient(gtheta2,e2);
dgtheta2de3 = gradient(gtheta2,e3);

ssu_theta2 = (1/gtheta2)^4*...
        ((gtheta2*dftheta2dtheta - ftheta2*dgtheta2dtheta)^2*ssu_theta + ...
        (gtheta2*dftheta2dk(1) - ftheta2*dgtheta2dk(1))^2*ssu_k(1) + ...
        (gtheta2*dftheta2dk(2) - ftheta2*dgtheta2dk(2))^2*ssu_k(2) + ...
        (gtheta2*dftheta2dk(3) - ftheta2*dgtheta2dk(3))^2*ssu_k(3) + ...
        (gtheta2*dftheta2de1(1) - ftheta2*dgtheta2de1(1))^2*ssu_e1(1) + ...
        (gtheta2*dftheta2de1(2) - ftheta2*dgtheta2de1(1))^2*ssu_e1(2) + ...
        (gtheta2*dftheta2de1(3) - ftheta2*dgtheta2de1(1))^2*ssu_e1(3) + ...
        (gtheta2*dftheta2de2(1) - ftheta2*dgtheta2de2(1))^2*ssu_e2(1) + ...
        (gtheta2*dftheta2de2(2) - ftheta2*dgtheta2de2(2))^2*ssu_e2(2) + ...
        (gtheta2*dftheta2de2(3) - ftheta2*dgtheta2de2(3))^2*ssu_e2(3) + ...
        (gtheta2*dftheta2de3(1) - ftheta2*dgtheta2de3(1))^2*ssu_e3(1) + ...
        (gtheta2*dftheta2de3(2) - ftheta2*dgtheta2de3(2))^2*ssu_e3(2) + ...
        (gtheta2*dftheta2de3(3) - ftheta2*dgtheta2de3(1))^2*ssu_e3(3));

% Theta3 ------------------------------------------------------------------
ftheta3 = dot(cross(e1,e2),k)*theta;
gtheta3 = dot(cross(e1,e2),e3);
%
dftheta3dtheta = gradient(ftheta3,theta)
dftheta3dk = gradient(ftheta3,k); % 0
dftheta3de1 = gradient(ftheta3,e1); 
dftheta3de2 = gradient(ftheta3,e2);
dftheta3de3 = gradient(ftheta3,e3);
%
dgtheta3dtheta = gradient(gtheta3,theta); % 0
dgtheta3dk = gradient(gtheta3,k); % 0
dgtheta3de1 = gradient(gtheta3,e1); 
dgtheta3de2 = gradient(gtheta3,e2);
dgtheta3de3 = gradient(gtheta3,e3);

ssu_theta3 = (1/gtheta3)^4*...
        ((gtheta3*dftheta3dtheta - ftheta3*dgtheta3dtheta)^2*ssu_theta + ...
        (gtheta3*dftheta3dk(1) - ftheta3*dgtheta3dk(1))^2*ssu_k(1) + ...
        (gtheta3*dftheta3dk(2) - ftheta3*dgtheta3dk(2))^2*ssu_k(2) + ...
        (gtheta3*dftheta3dk(3) - ftheta3*dgtheta3dk(3))^2*ssu_k(3) + ...
        (gtheta3*dftheta3de1(1) - ftheta3*dgtheta3de1(1))^2*ssu_e1(1) + ...
        (gtheta3*dftheta3de1(2) - ftheta3*dgtheta3de1(1))^2*ssu_e1(2) + ...
        (gtheta3*dftheta3de1(3) - ftheta3*dgtheta3de1(1))^2*ssu_e1(3) + ...
        (gtheta3*dftheta3de2(1) - ftheta3*dgtheta3de2(1))^2*ssu_e2(1) + ...
        (gtheta3*dftheta3de2(2) - ftheta3*dgtheta3de2(2))^2*ssu_e2(2) + ...
        (gtheta3*dftheta3de2(3) - ftheta3*dgtheta3de2(3))^2*ssu_e2(3) + ...
        (gtheta3*dftheta3de3(1) - ftheta3*dgtheta3de3(1))^2*ssu_e3(1) + ...
        (gtheta3*dftheta3de3(2) - ftheta3*dgtheta3de3(2))^2*ssu_e3(2) + ...
        (gtheta3*dftheta3de3(3) - ftheta3*dgtheta3de3(1))^2*ssu_e3(3));
    
    
%% 2. Load and calculate input variables theta, k, e1, e3
cd 'D:\Helical Axis\Test\'
addpath('D:\Helical Axis\Code')
% Select .c3d files 
[C3D_filenames, C3D_path, FilterIndex] = uigetfile('.c3d', 'Select .c3d files',['D:\Helical Axis\Data' '/'],'MultiSelect','on');
for file = 1:length(C3D_filenames) 
    Create_DataStruct(C3D_filenames{file}, C3D_path); % data Montreal
    %Create_DataStruct_2(C3D_filenames{file}, C3D_path); % data Geneve

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
    theta_v = acos(((R(1,1,:) + R(2,2,:) + R(3,3,:)) ... % Trace of R
        -1)/2); % Possible singularity when the angle is 0
    k_v = [(R(3,2,:)- R(2,3,:))./(2*sin(theta_v));...
        (R(1,3,:)-R(3,1,:))./(2*sin(theta_v));...
        (R(2,1,:)-R(1,2,:))./(2*sin(theta_v))]; % In thigh SCS
    kx_v=k_v(1,1,:); ky_v=k_v(2,1,:); kz_v=k_v(3,1,:);
    e1x_v=R(1,3,:); e1y_v=R(2,3,:); e1z_v=R(3,3,:);
    e3x_v=R(1,1,:); e3y_v=R(2,1,:); e3z_v=R(3,1,:);
    
    k_v = [kx_v, ky_v, kz_v];
    e1_v = [e1x_v, e1y_v, e1z_v];
    e3_v = [e3x_v, e3y_v, e3z_v];
    e2_v = zeros(1,3,length(e1x_v));
    for i=1:length(e1x_v)
        e2_v(:,:,i) = (cross(e3_v(:,:,i),e1_v(:,:,i)))/(norm(cross(e3_v(:,:,i),e1_v(:,:,i))));
    end
    e2x_v=e2_v(:,1,:); e2y_v=e2_v(:,2,:); e2z_v=e2_v(:,3,:);
    
%% 3. Replace input variables in the function of squared standard uncertainty (subs)

%     for i=1:size(k,3)

%         e1 = [e1x(:,:,i), e1y(:,:,i), e1z(:,:,i)];
%         e3 = [e3x(:,:,i), e3y(:,:,i), e3z(:,:,i)];
%         e2 = cross(e3,e1)/(norm(cross(e3,e1)))
% 
%         % Theta1
%         ftheta1 = dot(cross(e2,e3),k(:,:,i))*theta(i);
%         gtheta1 = dot(cross(e1,e2),e3);
%         %
%         dftheta1dtheta = gradient(ftheta1,theta(i))
%         dftheta1dk = [gradient(ftheta1,k(1,1,i)), gradient(ftheta1,k(2,1,i)), gradient(ftheta1,k(3,1,i))]; % 0
%         dftheta1de1 = [gradient(ftheta1,e1x), gradient(ftheta1,e1y), gradient(ftheta1,e1z)]; 
%         dftheta1de2 = [gradient(ftheta1,e2(1)), gradient(ftheta1,e2(2)), gradient(ftheta1,e2(3))];
%         dftheta1de3 = [gradient(ftheta1,e3x), gradient(ftheta1,e3y), gradient(ftheta1,e3z)];
%         %
%         dgtheta1dtheta = gradient(gtheta1,theta(i)); % 0
%         dgtheta1dk = [gradient(gtheta1,k(1,1,i)), gradient(gtheta1,k(2,1,i)), gradient(gtheta1,k(3,1,i))]; % 0
%         dgtheta1de1 = [gradient(gtheta1,e1x), gradient(gtheta1,e1y), gradient(gtheta1,e1z)]; 
%         dgtheta1de2 = [gradient(gtheta1,e2(1)), gradient(gtheta1,e2(2)), gradient(gtheta1,e2(3))];
%         dgtheta1de3 = [gradient(gtheta1,e3x), gradient(gtheta1,e3y), gradient(gtheta1,e3z)];

%     

    % Define ssu for each imput variable theta 1
    ssu_theta_v=1;
    ssu_kx_v =1; ssu_ky_v =1; ssu_kz_v =1;
    ssu_e1x_v=1; ssu_e1y_v=1; ssu_e1z_v=1;
    ssu_e2x_v=1; ssu_e2y_v=1; ssu_e2z_v=1;
    ssu_e3x_v=1; ssu_e3y_v=1; ssu_e3z_v=1;
    var_x = [0, 1, 0, 0, 0, 1, 1, 1];
    var_y = [0, 0, 1, 0, 1, 0, 1, 1];
    var_z = [0, 0, 0, 1, 1, 1, 0, 1];
    for t = 1:length(var_x)
        
        % altering e1
%         ssu_kx_v = var_x(t);
%         ssu_ky_v = var_y(t);
%         ssu_kz_v = var_z(t);
        
        
        % Define ssu for theta 1
        ssu_theta1_res = interpol(squeeze(double(subs(ssu_theta1, {theta,kx,ky,kz,e1x,e1y,e1z,e2x,e2y,e2z,e3x,e3y,e3z,...
            ssu_e1x,ssu_e1y,ssu_e1z,ssu_e2x,ssu_e2y,ssu_e2z,ssu_e3x,ssu_e3y,...
            ssu_e3z, ssu_theta, ssu_kx, ssu_ky, ssu_kz}, {theta_v,kx_v,ky_v,...
            kz_v,e1x_v,e1y_v,e1z_v,e2x_v,e2y_v,e2z_v,e3x_v,e3y_v,e3z_v,ssu_e1x_v,...
            ssu_e1y_v,ssu_e1z_v,ssu_e2x_v,ssu_e2y_v,ssu_e2z_v,ssu_e3x_v,ssu_e3y_v,...
            ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v, ssu_kz_v}))),101);

        % Define ssu for theta 2
        ssu_theta2_res = interpol(squeeze(double(subs(ssu_theta2, {theta,kx,ky,kz,e1x,e1y,e1z,e2x,e2y,e2z,e3x,e3y,e3z,...
            ssu_e1x,ssu_e1y,ssu_e1z,ssu_e2x,ssu_e2y,ssu_e2z,ssu_e3x,ssu_e3y,...
            ssu_e3z, ssu_theta, ssu_kx, ssu_ky, ssu_kz}, {theta_v,kx_v,ky_v,...
            kz_v,e1x_v,e1y_v,e1z_v,e2x_v,e2y_v,e2z_v,e3x_v,e3y_v,e3z_v,ssu_e1x_v,...
            ssu_e1y_v,ssu_e1z_v,ssu_e2x_v,ssu_e2y_v,ssu_e2z_v,ssu_e3x_v,ssu_e3y_v,...
            ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v, ssu_kz_v}))),101);
        % Define ssu for theta 3

        ssu_theta3_res = interpol(squeeze(double(subs(ssu_theta3, {theta,kx,ky,kz,e1x,e1y,e1z,e2x,e2y,e2z,e3x,e3y,e3z,...
            ssu_e1x,ssu_e1y,ssu_e1z,ssu_e2x,ssu_e2y,ssu_e2z,ssu_e3x,ssu_e3y,...
            ssu_e3z, ssu_theta, ssu_kx, ssu_ky, ssu_kz}, {theta_v,kx_v,ky_v,...
            kz_v,e1x_v,e1y_v,e1z_v,e2x_v,e2y_v,e2z_v,e3x_v,e3y_v,e3z_v,ssu_e1x_v,...
            ssu_e1y_v,ssu_e1z_v,ssu_e2x_v,ssu_e2y_v,ssu_e2z_v,ssu_e3x_v,ssu_e3y_v,...
            ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v, ssu_kz_v}))),101);
            %% 4. Plot theta1 results togheter with a theoretical corridor (mean_theta, std)

        figure(1)
        subplot(3,1,1)
        plot(ssu_theta1_res, 'b')%plot ssu_theta1
        hold on
        
        subplot(3,1,2)
        plot(ssu_theta2_res, 'b')%plot ssu_theta2
        hold on 
        subplot(3,1,3)
        plot(ssu_theta3_res, 'b')%plot ssu_theta3
        hold on
    end
    %plot theta1(mean)
    subplot(3,1,1)
    plot(interpol(squeeze(theta_v),101), 'r')
    title('ssu theta 1')
    xlim([0 100])
    subplot(3,1,2)
    plot(interpol(squeeze(theta_v),101), 'r')
    title('ssu theta 2')
    xlim([0 100])
    subplot(3,1,3)
    plot(interpol(squeeze(theta_v),101), 'r')
    title('ssu theta 3')
    xlim([0 100])
    suptitle('Altering ssu theta')

    %plot theta1(+/-std)
    
end