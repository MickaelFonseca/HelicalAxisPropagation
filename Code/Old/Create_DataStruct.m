% Create segment.mat containing the parameters Q
% Force plate (0), foot (1), shank (2), thigh (3) and pelvis (4)
% 
% Q: vector containing distal and proximal points, vectors u and w

% July 2019
% ________________________________________________________________________
function Create_DataStruct(C3D_filenames, C3D_path)
% clear all
% close all
% clc
% 
% % Select .c3d files 
% [C3D_filenames, C3D_path, FilterIndex] = uigetfile('.c3d', 'Select .c3d files',['D:\Helical Axis\Data' '/'],'MultiSelect','on');

C3D  = Get_C3D_BTK(C3D_path,C3D_filenames,1);
nf   = C3D.NumberOfFrames;
HEE  = C3D.data.LHEE; % heel
MET1 = C3D.data.LFMH; % 1st metatarsal head
MET2 = C3D.data.LSMH; % 2nd metatarsal head
MET5 = C3D.data.LVMH; % 5th metatarsal head
AJC  = C3D.data.LAJC; % ankle JC
KJC  = C3D.data.LKJC; % knee JC
LHJC = C3D.data.LHJC; % left hip JC
RHJC = C3D.data.RHJC; % right hip JC
KNE  = C3D.data.LKNE; % knee lateral epicondyle
KNI  = C3D.data.LKNI; % knee medial epicondyle
ANK = C3D.data.LANK; % ankle lateral malleolus
MED = C3D.data.LMED; % ankle medial malleolus
SACR = C3D.data.SACR; % medial posterior iliac spine
midHJC = C3D.data.midHJC; % medial point between hip JCs
midASIS = C3D.data.midASIS; % medial point between anterior iliac spines
TIB = C3D.data.LTIB;
% fp1_origin = C3D.data.fp1origin;

% % % 
% C3D=convert_analog_force_plate(C3D);
% resamp=1;
% resamp_fact=2;
% filtre=1;
% % Platform center of pressure 
% FP1=struct();
% FP2=struct();
% FPNet=struct();
% FP1.center=C3D.data.fp0origin;
% FP2.center=C3D.data.fp1origin;
% a=(C3D.analog.Force_Fz1==0);
% C3D.analog.Force_Fz1(a)=-1;
% b=(C3D.analog.Force_Fz2==0);
% C3D.analog.Force_Fz2(b)=-1;
% FP1.COP(:,1)=(C3D.analog.Moment_My1+C3D.analog.Force_Fx1*FP1.center(3))./C3D.analog.Force_Fz1+FP1.center(1);
% FP1.COP(:,2)=(C3D.analog.Moment_Mx1-C3D.analog.Force_Fy1*FP1.center(3))./C3D.analog.Force_Fz1+FP1.center(2);
% FP2.COP(:,1)=(C3D.analog.Moment_My2+C3D.analog.Force_Fx2*FP1.center(3))./C3D.analog.Force_Fz2+FP2.center(1);
% FP2.COP(:,2)=(C3D.analog.Moment_Mx2-C3D.analog.Force_Fy2*FP1.center(3))./C3D.analog.Force_Fz2+FP1.center(2);
%% 1. Q
ufoot  = (MET2-HEE)./norm(MET2-HEE); % calcaneus to 2nd metatarsal head
wfoot  = (MET5-MET1)./norm(MET5-MET1); % 1st to 5th metatarsal head
ushank = cross((AJC-KJC),(TIB-KJC)); % vector normal to plan (KJC, TIB, AJC)
ushank = ushank./norm(ushank);
wshank = (ANK-MED)./norm(ANK-MED); % medial to lateral malleolus
uthigh = cross((KNE-LHJC),(KNI-LHJC)); % vector normal to plan (KNE, LHJC, KNI)
uthigh = uthigh./norm(uthigh);
wthigh = (KNE-KNI)./norm(KNE-KNI);
upelvis = (midASIS-SACR)./norm(midASIS-SACR);
wpelvis = (RHJC-LHJC)./norm(RHJC-LHJC);

for i = 1:nf
%     Segment(1).Q(:,1,i) = [0,0,0,FP1.COP(i,:),0, fp1_origin(i,:), 0,0,0];

    Segment(2).Q(:,1,i) = [ufoot(i,:), AJC(i,:), MET2(i,:), wfoot(i,:)];

    Segment(3).Q(:,1,i) = [ushank(i,:), KJC(i,:), AJC(i,:), wshank(i,:)];
    
    Segment(4).Q(:,1,i) = [uthigh(i,:), LHJC(i,:), KJC(i,:), wthigh(i,:)];

    Segment(5).Q(:,1,i) = [upelvis(i,:), SACR(i,:), midHJC(i,:), wpelvis(i,:)];
end
cd 'D:\Helical Axis\Test\'
save('Segment.mat', 'Segment')
end