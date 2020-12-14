% Gait data
% Processing to obtain the rotation matrix from thigh to shank
% Based on:
% https://fr.mathworks.com/matlabcentral/fileexchange/58021-3d-kinematics-and-inverse-dynamics
clear all
close all
clc

cd 'D:\Helical Axis\Test\'
addpath('D:\Helical Axis\Code')
% Select .c3d files 
[C3D_filenames, C3D_path, FilterIndex] = uigetfile('.c3d', 'Select .c3d files',['D:\Helical Axis\Data' '/'],'MultiSelect','on');

title_joint = {'ankle', 'knee', 'hip'};
for file = 1:length(C3D_filenames) 
    Create_DataStruct(C3D_filenames{file}, C3D_path); % data Montreal
    %Create_DataStruct_2(C3D_filenames{file}, C3D_path); % data Geneve

    % Data example for gait
    load Segment.mat    
    Segment(2).T = Q2Tuv_array3(Segment(2).Q); % Foot (i = 2)
    Segment(3).T = Q2Tuv_array3(Segment(3).Q); % Shank (i = 3)
    Segment(4).T = Q2Tuv_array3(Segment(4).Q); % Thigh (i = 4)
    Segment(5).T = Q2Tuv_array3(Segment(5).Q); % Pelvis (i = 5)
        
    for i = 2:4
        if i == 2
            % Ankle joint
            Joint(2).T = Mprod_array3(Tinv_array3(Segment(3).T),Segment(2).T); % Ankle joint (i = 2)
            R = Joint(2).T(1:3,1:3,:);     
        elseif i == 3
            % Knee joint
            Joint(3).T = Mprod_array3(Tinv_array3(Segment(4).T),Segment(3).T); % Knee joint (i = 3)
            R = Joint(3).T(1:3,1:3,:);
        else 
            % Hip joint
            Joint(4).T = Mprod_array3(Tinv_array3(Segment(5).T),Segment(4).T); % Hip joint (i = 4)
            R = Joint(4).T(1:3,1:3,:);
        end

        %% Helical angle and axis
        theta = acos(((R(1,1,:) + R(2,2,:) + R(3,3,:)) ... % Trace of R
            -1)/2); % Possible singularity when the angle is 0
        k = [(R(3,2,:)- R(2,3,:))./(2*sin(theta));...
            (R(1,3,:)-R(3,1,:))./(2*sin(theta));...
            (R(2,1,:)-R(1,2,:))./(2*sin(theta))]; % In thigh SCS

        % Mean helical axis
        mk = mean(k,3);
        mk = mk/norm(mk); % Normed
        % Number of frames
        n = size(Joint(i).T,3);
        alpha = acos(dot(k,repmat(mk,[1,1,n]))); % Variation of axis orientation

        %% Figure
%         figure(i)
%         hold on
%         subplot(2,3,1)
%         plot(squeeze(theta*180/pi))
%         title ('Helical Angle (theta)')
%         xlabel('Sampled instant of time')
%         ylabel('Angle (in degree)')
%         subplot(2,3,2)
%         plot(squeeze(alpha*180/pi))
%         title ('Variation of Helical Axis w.r.t. Mean Axis');
%         xlabel('Sampled instant of time')
%         ylabel('Angle (in degree)')
%         suptitle(strcat(title_joint(i-1), '  joint'))

        %% Attitude vector about the Euler angle axes
        ktheraj = Vnop_array3(Mprod_array3(theta,k),... % to be projected
            repmat([0;0;1],[1,1,n]),... % Z axis of proximal segment
            Vnorm_array3(cross(R(1:3,2,:),repmat([0;0;1],[1,1,n]))),... X = YxZ
            R(1:3,2,:)); % Y axis of distal segment

        %% Euler angles
        Euler = R2mobileZXY_array3(R);

        %% Store data
        x = 1:length(ktheraj);
        THETA = 0;
        ALPHA = 0;
        for len = 1:length(ktheraj)
            THETA(len) = theta(1,1,len);
            ALPHA(len) = alpha(1,1,len);
        end
        if i == 2
            Joint_name = {'Ankle'};
            Data.theta(file).Theta_ankle = interp1(x', THETA', linspace(1, length(x), 101));
            Data.alpha(file).Alpha_ankle = interp1(x', ALPHA, linspace(1, length(x), 101));
        elseif i == 3
            Joint_name = {'Knee'};
            Data.theta(file).Theta_knee(:) = (interp1(x', THETA', linspace(1, length(x), 101)))';
            Data.alpha(file).Alpha_knee = interp1(x', ALPHA, linspace(1, length(x), 101));
        else
            Joint_name = {'Hip'};
            Data.theta(file).Theta_hip = interp1(x', THETA', linspace(1, length(x), 101));
            Data.alpha(file).Alpha_hip = interp1(x', ALPHA, linspace(1, length(x), 101));
        end

        for lin = 1:length(ktheraj)
            Ktheraj(lin,1) = ktheraj(1,1,lin)*180/pi; Ktheraj(lin,2) = ktheraj(2,1,lin)*180/pi; Ktheraj(lin,3) = ktheraj(3,1,lin)*180/pi;
            Eulert(lin,1)   = Euler(1,1,lin)*180/pi; Eulert(lin,2)   = Euler(1,2,lin)*180/pi; Eulert(lin,3)   = Euler(1,3,lin)*180/pi;
        end
        x = 1:length(Eulert);
        Data.(Joint_name{1})(file).HelicalAngles(:,1) = interp1(x',Ktheraj(:,1), linspace(1, length(x),101));
        Data.(Joint_name{1})(file).HelicalAngles(:,2) = interp1(x',Ktheraj(:,2), linspace(1, length(x),101));
        Data.(Joint_name{1})(file).HelicalAngles(:,3) = interp1(x',Ktheraj(:,3), linspace(1, length(x),101));
        Data.(Joint_name{1})(file).Euler(:,1) = interp1(x',Eulert(:,1), linspace(1, length(x), 101));
        Data.(Joint_name{1})(file).Euler(:,2) = interp1(x',Eulert(:,2), linspace(1, length(x), 101));
        Data.(Joint_name{1})(file).Euler(:,3) = interp1(x',Eulert(:,3), linspace(1, length(x), 101));
            
        %% Figure (cont.)
%         figure (i)
%         subplot(2,3,4)
%         hold on
%         plot(squeeze(ktheraj(1,1,:)*180/pi))
%         plot(squeeze(Euler(1,1,:)*180/pi))
%         title ('Flexion-Extension')
%         xlabel('Sampled instant of time')
%         ylabel('Angle (in degree)')
%         legend({'Projected ktheta', 'Euler angle'})
% 
%         subplot(2,3,5)
%         hold on
%         plot(squeeze(ktheraj(2,1,:)*180/pi))
%         plot(squeeze(Euler(1,2,:)*180/pi))
%         title ('Abduction-Adduction')
%         xlabel('Sampled instant of time')
%         ylabel('Angle (in degree)')
%         legend({'Projected ktheta', 'Euler angle'})
% 
%         subplot(2,3,6)
%         hold on
%         plot(squeeze(ktheraj(3,1,:)*180/pi))
%         plot(squeeze(Euler(1,3,:)*180/pi))
%         title ('Internal-External Rotation')
%         xlabel('Sampled instant of time')
%         ylabel('Angle (in degree)')
%         legend({'Projected ktheta', 'Euler angle'})
%         hold on

    end
end

%% Calculate Variability for KTHETA angles
%Mean +/- SD 
for i = 1:length(C3D_filenames)
   % Projected ktheta in 3 axis
   Res.HF_KT(:,i) = Data.Hip(i).HelicalAngles(:,1);
   Res.HA_KT(:,i) = Data.Hip(i).HelicalAngles(:,2);
   Res.HR_KT(:,i) = Data.Hip(i).HelicalAngles(:,3);
   Res.KF_KT(:,i) = Data.Knee(i).HelicalAngles(:,1);
   Res.KA_KT(:,i) = Data.Knee(i).HelicalAngles(:,2);
   Res.KR_KT(:,i) = Data.Knee(i).HelicalAngles(:,3);  
   Res.AF_KT(:,i) = Data.Ankle(i).HelicalAngles(:,1);
   Res.AA_KT(:,i) = Data.Ankle(i).HelicalAngles(:,2);
   Res.AR_KT(:,i) = Data.Ankle(i).HelicalAngles(:,3);  
   
   % Euler angles
   Res.HF_Eul(:,i) = Data.Hip(i).Euler(:,1);
   Res.HA_Eul(:,i) = Data.Hip(i).Euler(:,2);
   Res.HR_Eul(:,i) = Data.Hip(i).Euler(:,3);
   Res.KF_Eul(:,i) = Data.Knee(i).Euler(:,1);
   Res.KA_Eul(:,i) = Data.Knee(i).Euler(:,2);
   Res.KR_Eul(:,i) = Data.Knee(i).Euler(:,3);  
   Res.AF_Eul(:,i) = Data.Ankle(i).Euler(:,1);
   Res.AA_Eul(:,i) = Data.Ankle(i).Euler(:,2);
   Res.AR_Eul(:,i) = Data.Ankle(i).Euler(:,3);
   
   % Helical angle
   Res.THETA_ank(:,i) = (Data.theta(i).Theta_ankle)*180/pi;
   Res.THETA_knee(:,i) = (Data.theta(i).Theta_knee)*180/pi;
   Res.THETA_hip(:,i) = (Data.theta(i).Theta_hip)*180/pi;
   
   % Helical axis variation w.r.t. mean axis
   Res.ALPHA_ank(:,i) = (Data.alpha(i).Alpha_ankle)*180/pi;
   Res.ALPHA_knee(:,i) = (Data.alpha(i).Alpha_knee)*180/pi;
   Res.ALPHA_hip(:,i) = (Data.alpha(i).Alpha_hip)*180/pi;
end

Res.m_HF_KT  = mean(Res.HF_KT,2); Res.m_HA_KT  = mean(Res.HA_KT,2); Res.m_HR_KT  = mean(Res.HR_KT,2);
Res.sd_HF_KT = std(Res.HF_KT,1,2); Res.sd_HA_KT = std(Res.HA_KT,1,2);Res.sd_HR_KT = std(Res.HR_KT,1,2);
Res.m_KF_KT  = mean(Res.KF_KT,2); Res.m_KA_KT  = mean(Res.KA_KT,2); Res.m_KR_KT  = mean(Res.KR_KT,2);
Res.sd_KF_KT = std(Res.KF_KT,1,2); Res.sd_KA_KT = std(Res.KA_KT,1,2);Res.sd_KR_KT = std(Res.KR_KT,1,2);
Res.m_AF_KT  = mean(Res.AF_KT,2); Res.m_AA_KT  = mean(Res.AA_KT,2); Res.m_AR_KT  = mean(Res.AR_KT,2);
Res.sd_AF_KT = std(Res.AF_KT,1,2); Res.sd_AA_KT = std(Res.AA_KT,1,2);Res.sd_AR_KT = std(Res.AR_KT,1,2);

Res.m_HF_Eul  = mean(Res.HF_Eul,2); Res.m_HA_Eul  = mean(Res.HA_Eul,2); Res.m_HR_Eul  = mean(Res.HR_Eul,2);
Res.sd_HF_Eul = std(Res.HF_Eul,1,2); Res.sd_HA_Eul = std(Res.HA_Eul,1,2);Res.sd_HR_Eul = std(Res.HR_Eul,1,2);
Res.m_KF_Eul  = mean(Res.KF_Eul,2); Res.m_KA_Eul  = mean(Res.KA_Eul,2); Res.m_KR_Eul  = mean(Res.KR_Eul,2);
Res.sd_KF_Eul = std(Res.KF_Eul,1,2); Res.sd_KA_Eul = std(Res.KA_Eul,1,2);Res.sd_KR_Eul = std(Res.KR_Eul,1,2);
Res.m_AF_Eul  = mean(Res.AF_Eul,2); Res.m_AA_Eul  = mean(Res.AA_Eul,2); Res.m_AR_Eul  = mean(Res.AR_Eul,2);
Res.sd_AF_Eul = std(Res.AF_Eul,1,2); Res.sd_AA_Eul = std(Res.AA_Eul,1,2);Res.sd_AR_Eul = std(Res.AR_Eul,1,2);

Res.m_A_Ank = mean(Res.ALPHA_ank,2); Res.m_A_Knee = mean(Res.ALPHA_knee,2); Res.m_A_Hip = mean(Res.ALPHA_hip,2);
Res.sd_A_Ank = std(Res.ALPHA_ank,1,2);Res.sd_A_Knee = std(Res.ALPHA_knee,1,2);Res.sd_A_Hip = std(Res.ALPHA_hip,1,2);
Res.m_T_Ank = mean(Res.THETA_ank,2); Res.m_T_Knee = mean(Res.THETA_knee,2); Res.m_T_Hip = mean(Res.THETA_hip,2);
Res.sd_T_Ank = std(Res.THETA_ank,1,2);Res.sd_T_Knee = std(Res.THETA_knee,1,2);Res.sd_T_Hip = std(Res.THETA_hip,1,2);

figure(5)
plot_variability(Res, string('Ankle'))

figure(6)
plot_variability(Res, string('Knee'))

figure(7)
plot_variability(Res, string('Hip'))

%% Calculate RMSD 
[RMSD_Hip, m_RMSD_Hip] = RMSD_all(Res,C3D_filenames, 'Hip');
[RMSD_Knee, m_RMSD_Knee] = RMSD_all(Res,C3D_filenames, 'Knee');
[RMSD_Ankle, m_RMSD_Ankle] = RMSD_all(Res, C3D_filenames, 'Ankle');

T_RMSD = [m_RMSD_Hip; m_RMSD_Knee; m_RMSD_Ankle];
T_RMSD.Properties.RowNames = {'Hip', 'Knee', 'Ankle'};

writetable(T_RMSD, 'RMSD Variability IntraSession Geneve.xlsx', 'WriteRowNames', true)
%% Calculate LFT for alpha
[LFT.a1.A_Hip, LFT.a0.A_Hip, LFT.R2.A_Hip, LFT.Ya.A_Hip]= Linear_Fit_Method_Alpha(Res, 'hip', C3D_filenames);
[LFT.a1.A_Knee, LFT.a0.A_Knee, LFT.R2.A_Knee, LFT.Ya.A_Knee]= Linear_Fit_Method_Alpha(Res, 'knee', C3D_filenames);
[LFT.a1.A_Ankle, LFT.a0.A_Ankle, LFT.R2.A_Ankle, LFT.Ya.A_Ankle]= Linear_Fit_Method_Alpha(Res, 'ank', C3D_filenames);

%% Calculate LFT for Ktheta 
[LFT.a1.KT_HipRot, LFT.a0.KT_HipRot, LFT.R2.KT_HipRot, LFT.Ya.KT_HipRot]= Linear_Fit_Method_KTEul(Res, 'HR_KT', C3D_filenames);
[LFT.a1.KT_HipAdd, LFT.a0.KT_HipAdd, LFT.R2.KT_HipAdd, LFT.Ya.KT_HipAdd]= Linear_Fit_Method_KTEul(Res, 'HA_KT', C3D_filenames);
[LFT.a1.KT_HipFle, LFT.a0.KT_HipFle, LFT.R2.KT_HipFle, LFT.Ya.KT_HipFle]= Linear_Fit_Method_KTEul(Res, 'HF_KT', C3D_filenames);

[LFT.a1.KT_KneeRot, LFT.a0.KT_KneeRot, LFT.R2.KT_KneeRot, LFT.Ya.KT_KneeRot]= Linear_Fit_Method_KTEul(Res, 'KR_KT', C3D_filenames);
[LFT.a1.KT_KneeAdd, LFT.a0.KT_KneeAdd, LFT.R2.KT_KneeAdd, LFT.Ya.KT_KneeAdd]= Linear_Fit_Method_KTEul(Res, 'KA_KT', C3D_filenames);
[LFT.a1.KT_KneeFle, LFT.a0.KT_KneeFle, LFT.R2.KT_KneeFle, LFT.Ya.KT_KneeFle]= Linear_Fit_Method_KTEul(Res, 'KF_KT', C3D_filenames);

[LFT.a1.KT_AnkleRot, LFT.a0.KT_AnkleRot, LFT.R2.KT_AnkleRot, LFT.Ya.KT_AnkleRot]= Linear_Fit_Method_KTEul(Res, 'AR_KT', C3D_filenames);
[LFT.a1.KT_AnkleAdd, LFT.a0.KT_AnkleAdd, LFT.R2.KT_AnkleAdd, LFT.Ya.KT_AnkleAdd]= Linear_Fit_Method_KTEul(Res, 'AA_KT', C3D_filenames);
[LFT.a1.KT_AnkleFle, LFT.a0.KT_AnkleFle, LFT.R2.KT_AnkleFle, LFT.Ya.KT_AnkleFle]= Linear_Fit_Method_KTEul(Res, 'AF_KT', C3D_filenames);

%% Calculate LFT for Euler angles
[LFT.a1.Eul_HipRot, LFT.a0.Eul_HipRot, LFT.R2.Eul_HipRot, LFT.Ya.KT_HipRot]= Linear_Fit_Method_KTEul(Res, 'HR_Eul', C3D_filenames);
[LFT.a1.Eul_HipAdd, LFT.a0.Eul_HipAdd, LFT.R2.Eul_HipAdd, LFT.Ya.KT_HipAdd]= Linear_Fit_Method_KTEul(Res, 'HA_Eul', C3D_filenames);
[LFT.a1.Eul_HipFle, LFT.a0.Eul_HipFle, LFT.R2.Eul_HipFle, LFT.Ya.KT_HipFle]= Linear_Fit_Method_KTEul(Res, 'HF_Eul', C3D_filenames);

[LFT.a1.Eul_KneeRot, LFT.a0.Eul_KneeRot, LFT.R2.Eul_KneeRot, LFT.Ya.Eul_KneeRot]= Linear_Fit_Method_KTEul(Res, 'KR_Eul', C3D_filenames);
[LFT.a1.Eul_KneeAdd, LFT.a0.Eul_KneeAdd, LFT.R2.Eul_KneeAdd, LFT.Ya.Eul_KneeAdd]= Linear_Fit_Method_KTEul(Res, 'KA_Eul', C3D_filenames);
[LFT.a1.Eul_KneeFle, LFT.a0.Eul_KneeFle, LFT.R2.Eul_KneeFle, LFT.Ya.Eul_KneeFle]= Linear_Fit_Method_KTEul(Res, 'KF_Eul', C3D_filenames);

[LFT.a1.Eul_AnkleRot, LFT.a0.Eul_AnkleRot, LFT.R2.Eul_AnkleRot, LFT.Ya.Eul_AnkleRot]= Linear_Fit_Method_KTEul(Res, 'AR_Eul', C3D_filenames);
[LFT.a1.Eul_AnkleAdd, LFT.a0.Eul_AnkleAdd, LFT.R2.Eul_AnkleAdd, LFT.Ya.Eul_AnkleAdd]= Linear_Fit_Method_KTEul(Res, 'AA_Eul', C3D_filenames);
[LFT.a1.Eul_AnkleFle, LFT.a0.Eul_AnkleFle, LFT.R2.Eul_AnkleFle, LFT.Ya.Eul_AnkleFle]= Linear_Fit_Method_KTEul(Res, 'AF_Eul', C3D_filenames);

%% Compare LFT parameters




% plot mean + sd euler and ktheta angles
% figure(5)
% subplot(3,3,1)
% plot(m_HF_KT,'r')
% hold on
% corridor(m_HF_KT,sd_HF_KT, 'r', [0:100])
% plot(m_HF_Eul,'b')
% corridor(m_HF_Eul,sd_HF_Eul, 'b', [0:100])
% xlim([0 100])
% title('Hip Flex-Ext')
% 
% subplot(3,3,2)
% plot(m_HA_KT,'r')
% hold on
% corridor(m_HA_KT,sd_HA_KT, 'r', [0:100])
% plot(m_HA_Eul,'b')
% corridor(m_HA_Eul,sd_HA_Eul, 'b', [0:100])
% xlim([0 100])
% title('Hip Ad/Ab')
% 
% subplot(3,3,3)
% plot(m_HR_KT,'r')
% hold on
% corridor(m_HR_KT,sd_HR_KT, 'r', [0:100])
% plot(m_HR_Eul,'b')
% corridor(m_HR_Eul,sd_HR_Eul, 'b', [0:100])
% xlim([0 100])
% title('Hip Int/Ext Rot')
% legend('Ktheta', 'Euler')
% 
% subplot(3,3,4)
% plot(m_KF_KT,'r')
% hold on
% corridor(m_KF_KT,sd_KF_KT, 'r', [0:100])
% plot(m_KF_Eul,'b')
% corridor(m_KF_Eul,sd_KF_Eul, 'b', [0:100])
% xlim([0 100])
% title('Knee Flex-Ext')
% 
% subplot(3,3,5)
% plot(m_KA_KT,'r')
% hold on
% corridor(m_KA_KT,sd_KA_KT, 'r', [0:100])
% plot(m_KA_Eul,'b')
% corridor(m_KA_Eul,sd_KA_Eul, 'b', [0:100])
% xlim([0 100])
% title('Knee Ad/Ab')
% 
% subplot(3,3,6)
% plot(m_KR_KT,'r')
% hold on
% corridor(m_KR_KT,sd_KR_KT, 'r', [0:100])
% plot(m_KR_Eul,'b')
% corridor(m_KR_Eul,sd_KR_Eul, 'b', [0:100])
% xlim([0 100])
% title('Knee Int/Ext Rot')
% 
% subplot(3,3,7)
% plot(m_AF_KT,'r')
% hold on
% corridor(m_AF_KT,sd_AF_KT, 'r', [0:100])
% plot(m_AF_Eul,'b')
% corridor(m_AF_Eul,sd_AF_Eul, 'b', [0:100])
% xlim([0 100])
% title('Ankle Flex-Ext')
% 
% subplot(3,3,8)
% plot(m_AA_KT,'r')
% hold on
% corridor(m_AA_KT,sd_AA_KT, 'r', [0:100])
% plot(m_AA_Eul,'b')
% corridor(m_AA_Eul,sd_AA_Eul, 'b', [0:100])
% xlim([0 100])
% title('Ankle Ad/Ab')
% 
% subplot(3,3,9)
% plot(m_AR_KT,'r')
% hold on
% corridor(m_AR_KT,sd_AR_KT, 'r', [0:100])
% plot(m_AR_Eul,'b')
% corridor(m_AR_Eul,sd_AR_Eul, 'b', [0:100])
% xlim([0 100])
% title('Ankle Int/Ext Rot')
