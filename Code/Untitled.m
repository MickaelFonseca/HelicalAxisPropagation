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
for patient = 1:length(C3D_filenames) 
    Create_DataStruct(C3D_filenames{patient}, C3D_path);

    % Data example for gait
    load Segment.mat
    
    Segment(2).T = Q2Tuv_array3(Segment(2).Q); % Foot (i = 2)
    Segment(3).T = Q2Tuv_array3(Segment(3).Q); % Shank (i = 3)
    Segment(4).T = Q2Tuv_array3(Segment(4).Q); % Thigh (i = 4)
    Segment(5).T = Q2Tuv_array3(Segment(5).Q); % Pelvis (i = 5)
     
    % Ankle joint
    Joint(2).T = Mprod_array3(Tinv_array3(Segment(3).T),Segment(2).T); % Ankle joint (i = 2)
    RAnkle = Joint(2).T(1:3,1:3,:);
    
    % Knee joint
    Joint(3).T = Mprod_array3(Tinv_array3(Segment(4).T),Segment(3).T); % Knee joint (i = 3)
    RKnee = Joint(3).T(1:3,1:3,:);
    
    % Hip joint
    Joint(4).T = Mprod_array3(Tinv_array3(Segment(5).T),Segment(4).T); % Hip joint (i = 4)
    RHip = Joint(4).T(1:3,1:3,:);

        %% Helical angle and axis
        % Ankle
        thetaAnkle = acos(((RAnkle(1,1,:) + RAnkle(2,2,:) + RAnkle(3,3,:)) ... % Trace of R
            -1)/2); % Possible singularity when the angle is 0
        kAnkle = [(RAnkle(3,2,:)- RAnkle(2,3,:))./(2*sin(thetaAnkle));...
            (RAnkle(1,3,:)-RAnkle(3,1,:))./(2*sin(thetaAnkle));...
            (RAnkle(2,1,:)-RAnkle(1,2,:))./(2*sin(thetaAnkle))]; % In thigh SCS
        
        % Knee
        thetaKnee = acos(((RKnee(1,1,:) + RKnee(2,2,:) + RKnee(3,3,:)) ... % Trace of R
            -1)/2); % Possible singularity when the angle is 0
        kAnkle = [(RKnee(3,2,:)- RKnee(2,3,:))./(2*sin(thetaKnee));...
            (RKnee(1,3,:)-RKnee(3,1,:))./(2*sin(thetaKnee));...
            (RKnee(2,1,:)-RKnee(1,2,:))./(2*sin(thetaKnee))]; % In thigh SCS
        
        % Hip
        thetaHip = acos(((RHip(1,1,:) + RHip(2,2,:) + RHip(3,3,:)) ... % Trace of R
            -1)/2); % Possible singularity when the angle is 0
        kHip = [(RHip(3,2,:)- RHip(2,3,:))./(2*sin(thetaHip));...
            (RHip(1,3,:)-RHip(3,1,:))./(2*sin(thetaHip));...
            (RHip(2,1,:)-RHip(1,2,:))./(2*sin(thetaHip))]; % In thigh SCS

        
        % Mean helical axis
        mkAnkle = mean(kAnkle,3);
        mkAnkle = mkAnkle/norm(mkAnkle); % Normed
        mkKnee = mean(kKnee,3);
        mkKnee = mkKnee/norm(mkKnee); % Normed        
        mkHip = mean(kHip,3);
        mkHip = mkHip/norm(mkHip); % Normed
        
        % Number of frames
        n = size(Joint(2).T,3);
        alphaAnkle = acos(dot(k,repmat(mkAnkle,[1,1,n]))); % Variation of axis orientation
        alphaKnee = acos(dot(k,repmat(mkKnee,[1,1,n]))); % Variation of axis orientation
        alphaHip = acos(dot(k,repmat(mkHip,[1,1,n]))); % Variation of axis orientation
        
        %% Figure
        figure(i)
        subplot(2,3,1)
        plot(squeeze(theta*180/pi))
        title ('Helical Angle (theta)')
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')
        subplot(2,3,2)
        plot(squeeze(alpha*180/pi))
        title ('Variation of Helical Axis w.r.t. Mean Axis');
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')
        suptitle(strcat(title_joint(i-1), '  joint'))


        %% Attitude vector about the Euler angle axes
        ktheraj = Vnop_array3(...
            Mprod_array3(theta,k),... % to be projected
            repmat([0;0;1],[1,1,n]),... % Z axis of proximal segment
            Vnorm_array3(cross(R(1:3,2,:),repmat([0;0;1],[1,1,n]))),... X = YxZ
            R(1:3,2,:)); % Y axis of distal segment


        %% Euler angles
        Euler = R2mobileZXY_array3(R);

        %% Store data
        if C3D_filenames{patient}(9:11) == 'S01'
            if C3D_filenames{patient}(5:8) == 'OP01'
                Data(patient).Exam(1).Session(1).Hip.HelicalAngles = ktheraj
                Data(patient).Exam(1).Session(1).Knee.HelicalAngles = ktheraj
                Data(patient).Exam(1).Session(1).Ankle.HelicalAngles = ktheraj

                Data(patient).Exam(1).Session(1).Hip.Euler = Euler
                Data(patient).Exam(1).Session(1).Knee.Euler = Euler
                Data(patient).Exam(1).Session(1).Ankle.Euler = Euler
                
            elseif C3D_filenames{patient}(5:8) == 'OP02'
                Data(patient).Exam(2).Session(1).Joint.HelicalAngles = ktheraj
                Data(patient).Exam(2).Session(1).Joint.Euler = Euler
            else C3D_filenames{patient}.Exam(3).Session == 'OP03'
                Data(patient).Exam(3).Session(1).Joint.HelicalAngles = ktheraj
                Data(patient).Exam(3).Session(1).Joint.Euler = Euler
            end
        elseif C3D_filenames{patient}(9:11) == 'S02'
            if C3D_filenames{patient}(5:8) == 'OP01'
                Data(patient).Exam(1).Session(2).Joint.HelicalAngles = ktheraj
                Data(patient).Exam(1).Session(2).Joint.Euler = Euler
            elseif C3D_filenames{patient}(5:8) == 'OP02'
                Data(patient).Exam(2).Session(2).Joint.HelicalAngles = ktheraj
                Data(patient).Exam(2).Session(2).Joint.Euler = Euler
            else C3D_filenames{patient}.Exam(3).Session == 'OP03'
                Data(patient).Exam(3).Session(2).Joint.HelicalAngles = ktheraj
                Data(patient).Exam(3).Session(2).Joint.Euler = Euler
            end
        end
            
                
        % Data(patient).Exam(examiner).Session(1 or 2).Joint.Euler = Euler
        %% Figure (cont.)
        subplot(2,3,4)
        hold on
        plot(squeeze(ktheraj(1,1,:)*180/pi))
        plot(squeeze(Euler(1,1,:)*180/pi))
        title ('Flexion-Extension')
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')
        % legend({'Projected ktheta', 'Euler angle'})

        subplot(2,3,5)
        hold on
        plot(squeeze(ktheraj(2,1,:)*180/pi))
        plot(squeeze(Euler(1,2,:)*180/pi))
        title ('Abduction-Adduction')
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')
        % legend({'Projected ktheta', 'Euler angle'})

        subplot(2,3,6)
        hold on
        plot(squeeze(ktheraj(3,1,:)*180/pi))
        plot(squeeze(Euler(1,3,:)*180/pi))
        title ('Internal-External Rotation')
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')
        legend({'Projected ktheta', 'Euler angle'})

    end
end