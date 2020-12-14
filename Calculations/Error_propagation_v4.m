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
% intrinsic vs extrinsic ssu
ssui = [0,1,5];
ssue = [0,1,5];

% Tssui = []; Tssue = [];
% Tssutheta1_m = []; Tssutheta2_m = []; Tssutheta3_m=[];
% Tssutheta1_std = []; Tssutheta2_std = []; Tssutheta3_std=[];
% CU = [];
% EU = [];
% exp_factor = 2;
% for i=1:length(ssui)
%     for e=1:length(ssue)
%         
%         % Uncertainty
%         % Squared standard uncertainty ((standard deviation)^2)
%         ui = deg2rad(ssui(i))^2; % Intrinsic 
%         ue = deg2rad(ssue(e))^2; % Extrinsic
% 
%         % Define ssu for each imput variable
%         % Intrinsic
%         ssu_theta_v = ui;
%         ssu_kx_v = tan(ui); ssu_ky_v = tan(ui);
%         % Extrinsic
%         ssu_e1x_v = tan(ue); ssu_e1y_v = tan(ue);
%         ssu_e3x_v =tan(ue); ssu_e3z_v = tan(ue);
% 
%         % Define ssu for theta 1
%         ssu_theta1_res = double(subs(ssu_theta1, ...
%             {theta,kx,ky,e1x,e1y,e3x,e3z,...
%             ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
%             {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
%             ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));
% 
%         % Define ssu for theta 2
%         ssu_theta2_res = double(subs(ssu_theta2, ...
%             {theta,kx,ky,e1x,e1y,e3x,e3z,...
%             ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
%             {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
%             ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));
% 
%         % Define ssu for theta 3
%         ssu_theta3_res = double(subs(ssu_theta3, ...
%             {theta,kx,ky,e1x,e1y,e3x,e3z,...
%             ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
%             {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
%             ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));
% 
% 
%         %% 4. Plot theta1 results togheter with a theoretical corridor (mean_theta, std)
% 
%         gcf = figure(1)
%         %
%         subplot(3,1,1)
%         A = rad2deg(mean(theta1_res,2));
%         B = rad2deg(std(theta1_res,[],2));
%         Max = A+B; Min = A-B;
%         x = [1:101 101:-1:1];
%         y = [Max;Min(end:-1:1)];
%         hold on
%         plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
%         fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%         xlim([0 101]);
%         B = rad2deg(sqrt(ssu_theta1_res)); % Standard uncertainty
%         Max = A+B; Min = A-B;
%         x = [1:101 101:-1:1];
%         y = [Max;Min(end:-1:1)];
%         fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%         ylabel('Theta 1')
%         str = sprintf('Theoritical and Experimental Variabilities ssui:%d ssue: %d', [ssui(i), ssue(e)])
%         legend({'Experimental (mean)','Experimental (SD)','Theoritical (SU)'}, 'Location', 'east')
%         title(str)
%         
% 
%         subplot(3,1,2)
%         A = rad2deg(mean(theta2_res,2));
%         B = rad2deg(std(theta2_res,[],2));
%         Max = A+B; Min = A-B;
%         x = [1:101 101:-1:1];
%         y = [Max;Min(end:-1:1)];
%         hold on
%         plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
%         fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%         xlim([0 101]);
%         B = rad2deg(sqrt(ssu_theta2_res)); % Standard uncertainty
%         Max = A+B; Min = A-B;
%         x = [1:101 101:-1:1];
%         y = [Max;Min(end:-1:1)];
%         fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%         ylabel('Theta 2')
% 
%         subplot(3,1,3)
%         A = rad2deg(mean(theta3_res,2));
%         B = rad2deg(std(theta3_res,[],2));
%         Max = A+B; Min = A-B;
%         x = [1:101 101:-1:1];
%         y = [Max;Min(end:-1:1)];
%         hold on
%         plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
%         fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%         xlim([0 101]);
%         B = rad2deg(sqrt(ssu_theta3_res)); % Standard uncertainty
%         Max = A+B; Min = A-B;
%         x = [1:101 101:-1:1];
%         y = [Max;Min(end:-1:1)];
%         fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%         ylabel('Theta 3')
%         xlabel('% of Gait Cycle')
%         strsvg = sprintf('ssui_%d_ssue_%d', [ssui(i), ssue(e)]);
%         saveas(gcf, [strsvg '.png'], 'png')
%         close(gcf)
%         
%         Tssui = [Tssui; ssui(i)];
%         Tssue = [Tssue; ssue(e)];
%         Tssutheta1_m = [Tssutheta1_m; rad2deg(mean(sqrt(ssu_theta1_res)))];
%         Tssutheta1_std = [Tssutheta1_std; rad2deg(std(sqrt(ssu_theta1_res)))];
%         Tssutheta2_m = [Tssutheta2_m; rad2deg(mean(sqrt(ssu_theta2_res)))];
%         Tssutheta2_std = [Tssutheta2_std; rad2deg(std(sqrt(ssu_theta2_res)))];
%         Tssutheta3_m = [Tssutheta3_m; rad2deg(mean(sqrt(ssu_theta3_res)))];
%         Tssutheta3_std = [Tssutheta3_std; rad2deg(std(sqrt(ssu_theta3_res)))];
%         %Combined standard uncertainty
%         CU = [CU; sqrt(rad2deg(mean(sqrt(ssu_theta1_res))).^2 + rad2deg(mean(sqrt(ssu_theta2_res))).^2 + rad2deg(mean(sqrt(ssu_theta3_res))).^2)];          
%         % Expanded Uncertainty
%         EU = [EU; exp_factor*sqrt(rad2deg(mean(sqrt(ssu_theta1_res))).^2 + rad2deg(mean(sqrt(ssu_theta2_res))).^2 + rad2deg(mean(sqrt(ssu_theta3_res))).^2)];
%     end
% end
% 
% R = table(Tssui, Tssue, Tssutheta1_m, Tssutheta1_std, Tssutheta2_m, Tssutheta2_std, Tssutheta3_m, Tssutheta3_std, CU, EU);
% writetable(R, 'standard uncertainties.xlsx', 'Sheet',1)

%% 4. Replace input variables in the function of squared standard uncertainty ssui ssue1 ssue3 (subs)
ssui = [0];
ssue1 = [0,1,5];
ssue3 = [0,1,5];

Tssui = []; Tssue1 = []; Tssue3 = [];
Tssutheta1_m = []; Tssutheta2_m = []; Tssutheta3_m=[];
Tssutheta1_std = []; Tssutheta2_std = []; Tssutheta3_std=[];
CU = [];
EU = [];
exp_factor = 2;
for i=1:length(ssui)
    for e1=1:length(ssue1)
        for e3 = 1:length(ssue3)
        
            % Uncertainty
            % Squared standard uncertainty ((standard deviation)^2)
            ui = deg2rad(ssui(i))^2; % Intrinsic 
            ue1 = deg2rad(ssue(e1))^2; % Extrinsic
            ue3 = deg2rad(ssue(e3))^2;
            
            % Define ssu for each imput variable
            % Intrinsic
            ssu_theta_v = ui;
            ssu_kx_v = tan(ui); ssu_ky_v = tan(ui);
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

            figure(1)
            %
            subplot(3,1,1)
            A = rad2deg(mean(theta1_res,2));
            B = rad2deg(std(theta1_res,[],2));
            Max = A+B; Min = A-B;
            x = [1:101 101:-1:1];
            y = [Max;Min(end:-1:1)];
            hold on
            plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
            fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
            xlim([0 101]);
            B = rad2deg(sqrt(ssu_theta1_res)); % Standard uncertainty
            Max = A+B; Min = A-B;
            x = [1:101 101:-1:1];
            y = [Max;Min(end:-1:1)];
            fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
            ylabel('Theta 1')
            str = sprintf('Theoritical and Experimental Variabilities ssui:%d ssue1: %d ssue3: %d', [ssui, ssue1(e1), ssue3(e3)])
            legend({'Experimental (mean)','Experimental (SD)','Theoritical (SU)'}, 'Location', 'east')
            title(str)
            
            subplot(3,1,2)
            A = rad2deg(mean(theta2_res,2));
            B = rad2deg(std(theta2_res,[],2));
            Max = A+B; Min = A-B;
            x = [1:101 101:-1:1];
            y = [Max;Min(end:-1:1)];
            hold on
            plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
            fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
            xlim([0 101]);
            B = rad2deg(sqrt(ssu_theta2_res)); % Standard uncertainty
            Max = A+B; Min = A-B;
            x = [1:101 101:-1:1];
            y = [Max;Min(end:-1:1)];
            fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
            ylabel('Theta 2')

            subplot(3,1,3)
            A = rad2deg(mean(theta3_res,2));
            B = rad2deg(std(theta3_res,[],2));
            Max = A+B; Min = A-B;
            x = [1:101 101:-1:1];
            y = [Max;Min(end:-1:1)];
            hold on
            plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
            fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
            xlim([0 101]);
            B = rad2deg(sqrt(ssu_theta3_res)); % Standard uncertainty
            Max = A+B; Min = A-B;
            x = [1:101 101:-1:1];
            y = [Max;Min(end:-1:1)];
            fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
            ylabel('Theta 3')
            xlabel('% of Gait Cycle')
            strsvg = sprintf('ssui_%d_ssue1_%d_ssue3_%d', [ssui, ssue1(e1), ssue3(e3)]);
            saveas(gcf, [strsvg '.png'], 'png')
            close(gcf)
            
            Tssui = [Tssui; ssui(i)];
            Tssue1 = [Tssue1; ssue1(e1)];
            Tssue3 = [Tssue3; ssue3(e3)];
            Tssutheta1_m = [Tssutheta1_m; rad2deg(mean(sqrt(ssu_theta1_res)))];
            Tssutheta1_std = [Tssutheta1_std; rad2deg(std(sqrt(ssu_theta1_res)))];
            Tssutheta2_m = [Tssutheta2_m; rad2deg(mean(sqrt(ssu_theta2_res)))];
            Tssutheta2_std = [Tssutheta2_std; rad2deg(std(sqrt(ssu_theta2_res)))];
            Tssutheta3_m = [Tssutheta3_m; rad2deg(mean(sqrt(ssu_theta3_res)))];
            Tssutheta3_std = [Tssutheta3_std; rad2deg(std(sqrt(ssu_theta3_res)))];
            %Combined standard uncertainty
            CU = [CU; sqrt(rad2deg(mean(sqrt(ssu_theta1_res))).^2 + rad2deg(mean(sqrt(ssu_theta2_res))).^2 + rad2deg(mean(sqrt(ssu_theta3_res))).^2)];          
            % Expanded Uncertainty
            EU = [EU; exp_factor*sqrt(rad2deg(mean(sqrt(ssu_theta1_res))).^2 + rad2deg(mean(sqrt(ssu_theta2_res))).^2 + rad2deg(mean(sqrt(ssu_theta3_res))).^2)];
        end
    end
end

R2 = table(Tssui, Tssue1, Tssue3, Tssutheta1_m, Tssutheta1_std, Tssutheta2_m, Tssutheta2_std, Tssutheta3_m, Tssutheta3_std, CU, EU);
writetable(R2, 'standard uncertainties.xlsx', 'Sheet',2)     

% %% 5. Replace input variables in the function of squared standard uncertainty ssui ssue1x ssue1y ssue3x ssue3z (subs)
% ssui = [0,1,3,5];
% ssue1x = [0,1,3,5];
% ssue1y = [0,1,3,5];
% ssue3x = [0,1,3,5];
% ssue3z = [0,0,1,3,5];
% 
% Tssui = []; Tssue1x = []; Tssue1y = []; Tssue3x = []; Tssue3z = [];
% Tssutheta1_m = []; Tssutheta2_m = []; Tssutheta3_m=[];
% Tssutheta1_std = []; Tssutheta2_std = []; Tssutheta3_std=[];
% CU = [];
% EU = [];
% exp_factor = 2;
% for i=1:length(ssui)
%     for e1x=1:length(ssue1x)
%         for e1y=1:length(ssue1y)
%             for e3x=1:length(ssue3x)
%                 for e3z = 1:length(ssue3z)
%         
%                     % Uncertainty
%                     % Squared standard uncertainty ((standard deviation)^2)
%                     ui = deg2rad(ssui(i))^2; % Intrinsic 
%                     ue1x = deg2rad(ssue(e1x))^2; % Extrinsic
%                     ue1y = deg2rad(ssue(e1y))^2;
%                     ue3x = deg2rad(ssue(e3x))^2;
%                     ue3z = deg2rad(ssue(e3z))^2;
%                     disp(ue3z)
%                     
%                     % Define ssu for each imput variable
%                     % Intrinsic
%                     ssu_theta_v = ui;
%                     ssu_kx_v = tan(ui); ssu_ky_v = tan(ui);
%                     % Extrinsic
%                     ssu_e1x_v = tan(ue1x); ssu_e1y_v = tan(ue1y);
%                     ssu_e3x_v =tan(ue3x); ssu_e3z_v = tan(ue3z);
% 
%                     % Define ssu for theta 1
%                     ssu_theta1_res = double(subs(ssu_theta1, ...
%                         {theta,kx,ky,e1x,e1y,e3x,e3z,...
%                         ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
%                         {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
%                         ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));
% 
%                     % Define ssu for theta 2
%                     ssu_theta2_res = double(subs(ssu_theta2, ...
%                         {theta,kx,ky,e1x,e1y,e3x,e3z,...
%                         ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
%                         {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
%                         ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));
% 
%                     % Define ssu for theta 3
%                     ssu_theta3_res = double(subs(ssu_theta3, ...
%                         {theta,kx,ky,e1x,e1y,e3x,e3z,...
%                         ssu_e1x,ssu_e1y,ssu_e3x,ssu_e3z, ssu_theta, ssu_kx, ssu_ky}, ...
%                         {mean(theta_v,2),mean(kx_v,2),mean(ky_v,2),mean(e1x_v,2),mean(e1y_v,2),mean(e3x_v,2),mean(e3z_v,2),...
%                         ssu_e1x_v,ssu_e1y_v,ssu_e3x_v,ssu_e3z_v, ssu_theta_v, ssu_kx_v, ssu_ky_v}));
% 
% 
%                     %% 4. Plot theta1 results togheter with a theoretical corridor (mean_theta, std)
% 
%                     figure(1)
%                     %
%                     subplot(3,1,1)
%                     A = rad2deg(mean(theta1_res,2));
%                     B = rad2deg(std(theta1_res,[],2));
%                     Max = A+B; Min = A-B;
%                     x = [1:101 101:-1:1];
%                     y = [Max;Min(end:-1:1)];
%                     hold on
%                     plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
%                     fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%                     xlim([0 101]);
%                     B = rad2deg(sqrt(ssu_theta1_res)); % Standard uncertainty
%                     Max = A+B; Min = A-B;
%                     x = [1:101 101:-1:1];
%                     y = [Max;Min(end:-1:1)];
%                     fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%                     ylabel('Theta 1')
%                     title('Theoritical and Experimental Variabilities')
%         
%                     subplot(3,1,2)
%                     A = rad2deg(mean(theta2_res,2));
%                     B = rad2deg(std(theta2_res,[],2));
%                     Max = A+B; Min = A-B;
%                     x = [1:101 101:-1:1];
%                     y = [Max;Min(end:-1:1)];
%                     hold on
%                     plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
%                     fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%                     xlim([0 101]);
%                     B = rad2deg(sqrt(ssu_theta2_res)); % Standard uncertainty
%                     Max = A+B; Min = A-B;
%                     x = [1:101 101:-1:1];
%                     y = [Max;Min(end:-1:1)];
%                     fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%                     ylabel('Theta 2')
%         
%                     subplot(3,1,3)
%                     A = rad2deg(mean(theta3_res,2));
%                     B = rad2deg(std(theta3_res,[],2));
%                     Max = A+B; Min = A-B;
%                     x = [1:101 101:-1:1];
%                     y = [Max;Min(end:-1:1)];
%                     hold on
%                     plot(A,'Color',[0 0.4470 0.7410],'linestyle','-','Marker','none');
%                     fill(x,y,[0 0.4470 0.7410],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%                     xlim([0 101]);
%                     B = rad2deg(sqrt(ssu_theta3_res)); % Standard uncertainty
%                     Max = A+B; Min = A-B;
%                     x = [1:101 101:-1:1];
%                     y = [Max;Min(end:-1:1)];
%                     fill(x,y,[0.8500 0.3250 0.0980],'linestyle','-','edgecolor','none','faceAlpha',0.5)
%                     ylabel('Theta 3')
%                     xlabel('% of Gait Cycle')
%                     legend({'Experimental (mean)','Experimental (SD)','Theoritical (SU)'})
% 
%                     Tssui = [Tssui; ssui(i)];
%                     Tssue1x = [Tssue1x; ssue1x(e1x)];
%                     Tssue1y = [Tssue1y; ssue1y(e1y)];
%                     Tssue3x = [Tssue3x; ssue3x(e3x)];
%                     Tssue3z = [Tssue3z; ssue3z(e3z)];
%                     Tssutheta1_m = [Tssutheta1_m; rad2deg(mean(sqrt(ssu_theta1_res)))];
%                     Tssutheta1_std = [Tssutheta1_std; rad2deg(std(sqrt(ssu_theta1_res)))];
%                     Tssutheta2_m = [Tssutheta2_m; rad2deg(mean(sqrt(ssu_theta2_res)))];
%                     Tssutheta2_std = [Tssutheta2_std; rad2deg(std(sqrt(ssu_theta2_res)))];
%                     Tssutheta3_m = [Tssutheta3_m; rad2deg(mean(sqrt(ssu_theta3_res)))];
%                     Tssutheta3_std = [Tssutheta3_std; rad2deg(std(sqrt(ssu_theta3_res)))];
%                     %Combined standard uncertainty
%                     CU = [CU; sqrt(rad2deg(mean(sqrt(ssu_theta1_res))).^2 + rad2deg(mean(sqrt(ssu_theta2_res))).^2 + rad2deg(mean(sqrt(ssu_theta3_res))).^2)];          
%                     % Expanded Uncertainty
%                     EU = [EU; exp_factor*sqrt(rad2deg(mean(sqrt(ssu_theta1_res))).^2 + rad2deg(mean(sqrt(ssu_theta2_res))).^2 + rad2deg(mean(sqrt(ssu_theta3_res))).^2)];
%                 end
%             end
%         end    
%     end
% end
% % 
% R3 = table(Tssui, Tssue1x, Tssue1y, Tssue3x, Tssue3z, Tssutheta1_m, Tssutheta1_std, Tssutheta2_m, Tssutheta2_std, Tssutheta3_m, Tssutheta3_std, CU, EU);
% writetable(R3, 'standard uncertainties.xlsx', 'Sheet',3)