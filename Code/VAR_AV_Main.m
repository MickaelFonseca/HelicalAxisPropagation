% Study VARAV - Prove that VAR in gait comes from more from the axis
% definition than from the movement itself
% Author: MFonseca # July 2019
close all
clc

data_path  = 'D:\Helical Axis\Test\';
file   = 'PN01OP01S01SS01_CGM1_1.c3d';

acq = btkReadAcquisition(strcat(data_path, file));
markers = btkGetMarkers(acq);
angles  = btkGetAngles(acq);


%% 1. Extract rotation matrix from Euler angles

%% 2. Extract k and theta(O) from rotation matrix

%% 3. Calculate mean of k among all frames

%% 4. Compute angle between mean k and k for all frames 

%% 5. Normalize it and calculate .product (= cos(angle) between both vectors)


%% 6. Inverse cos(angle) give the helical axis


%% 7. Project KO to the Euler angles


%% 8. Calculate VAR between Euler angles and Projected KO to the Euler angles


%% 9. Compare VAR and export data 