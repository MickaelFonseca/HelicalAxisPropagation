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
        figure(i)
        hold on
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
        if i == 2
            Joint_name = {'Ankle'};
        elseif i == 3
            Joint_name = {'Knee'};
        else
            Joint_name = {'Hip'};
        end
        Patient_name = string(C3D_filenames{patient}(1:4)); 
        if C3D_filenames{patient}(9:11) == 'S01'
            if C3D_filenames{patient}(5:8) == 'OP01'
                Data.(Patient_name{1}).Exam(1).Session(1).(Joint_name{1}).HelicalAngles = ktheraj
                Data.(Patient_name{1}).Exam(1).Session(1).(Joint_name{1}).Euler = Euler
            elseif C3D_filenames{patient}(5:8) == 'OP02'
                Data.(Patient_name{1}).Exam(2).Session(1).(Joint_name{1}).HelicalAngles = ktheraj
                Data.(Patient_name{1}).Exam(2).Session(1).(Joint_name{1}).Euler = Euler
            elseif C3D_filenames{patient}(5:8) == 'OP03'
                Data.(Patient_name{1}).Exam(3).Session(1).(Joint_name{1}).HelicalAngles = ktheraj
                Data.(Patient_name{1}).Exam(3).Session(1).(Joint_name{1}).Euler = Euler
            end
        elseif C3D_filenames{patient}(9:11) == 'S02'
            if C3D_filenames{patient}(5:8) == 'OP01'
                Data.(Patient_name{1}).Exam(1).Session(2).(Joint_name{1}).HelicalAngles = ktheraj
                Data.(Patient_name{1}).Exam(1).Session(2).(Joint_name{1}).Euler = Euler
            elseif C3D_filenames{patient}(5:8) == 'OP02'
                Data.(Patient_name{1}).Exam(2).Session(2).(Joint_name{1}).HelicalAngles = ktheraj
                Data.(Patient_name{1}).Exam(2).Session(2).(Joint_name{1}).Euler = Euler
            elseif C3D_filenames{patient}(5:8) == 'OP03'
                Data.(Patient_name{1}).Exam(2).Session(2).(Joint_name{1}).HelicalAngles = ktheraj
                Data.(Patient_name{1}).Exam(2).Session(2).(Joint_name{1}).Euler = Euler
            end
        end
            
        %% Figure (cont.)
        figure (i)
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
        hold on

    end
end
