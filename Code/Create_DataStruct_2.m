% Create segment.mat containing the parameters Q
% Force plate (0), foot (1), shank (2), thigh (3) and pelvis (4)
% 
% Q: vector containing distal and proximal points, vectors u and w

% July 2019
% ________________________________________________________________________
function Create_DataStruct_2(C3D_filenames, C3D_path) % For data from lab Geneve
C3D  = Get_C3D_BTK(C3D_path,C3D_filenames,1);
% C3D_filenames = [C3D_path,C3D_filenames];
% acq = btkReadAcquisition(C3D_filenames);
% C3D.data = btkGetMarkers(acq);
nf = length(C3D.data.LHEE);

PN = C3D.filename(1:4);
% Select first cycle (left)
if isfield(C3D.EventFrame, 'Left_Foot_Strike')==1
    ev_1 = C3D.EventFrame.Left_Foot_Strike(1);
    ev_2 = C3D.EventFrame.Left_Foot_Strike(2);
elseif isfield(C3D.EventFrame, strcat(PN, '_Left_Foot_Strike')) == 1
    ev_1 = C3D.EventFrame.(strcat(PN, '_Left_Foot_Strike'))(1);
    ev_2 = C3D.EventFrame.(strcat(PN, '_Left_Foot_Strike'))(2);
end
% XYZ (Qualisys) => XZ(-Y) (ISB) in m
% Left leg
HEE  = permute([C3D.data.LHEE(ev_1:ev_2,1),C3D.data.LHEE(ev_1:ev_2,3),-C3D.data.LHEE(ev_1:ev_2,2)]', [1,3,2])/1000; % heel
MET1 = permute([C3D.data.LFMH(ev_1:ev_2,1),C3D.data.LFMH(ev_1:ev_2,3),-C3D.data.LFMH(ev_1:ev_2,2)]', [1,3,2])/1000; % 1st metatarsal head
MET2 = permute([C3D.data.LTOE(ev_1:ev_2,1),C3D.data.LTOE(ev_1:ev_2,3),-C3D.data.LTOE(ev_1:ev_2,2)]', [1,3,2])/1000; % 2nd metatarsal head
if isfield(C3D.data, 'LMET')== 1
    MET5 = permute([C3D.data.LMET(ev_1:ev_2,1),C3D.data.LMET(ev_1:ev_2,3),- C3D.data.LMET(ev_1:ev_2,2)]', [1,3,2])/1000; % 5th metatarsal head
elseif isfield(C3D.data, 'LVMH') == 1
    MET5 = permute([C3D.data.LVMH(ev_1:ev_2,1),C3D.data.LVMH(ev_1:ev_2,3),- C3D.data.LVMH(ev_1:ev_2,2)]', [1,3,2])/1000; 
end
AJC   = permute([C3D.data.LAJC(ev_1:ev_2,1),C3D.data.LAJC(ev_1:ev_2,3),-C3D.data.LAJC(ev_1:ev_2,2)]', [1,3,2])/1000; % ankle JC
TIB   = permute([C3D.data.LTIB(ev_1:ev_2,1),C3D.data.LTIB(ev_1:ev_2,3),-C3D.data.LTIB(ev_1:ev_2,2)]', [1,3,2])/1000; % wond thigh
KJC   = permute([C3D.data.LKJC(ev_1:ev_2,1),C3D.data.LKJC(ev_1:ev_2,3),-C3D.data.LKJC(ev_1:ev_2,2)]', [1,3,2])/1000; % knee JC
LHJC  = permute([C3D.data.LHJC(ev_1:ev_2,1),C3D.data.LHJC(ev_1:ev_2,3),-C3D.data.LHJC(ev_1:ev_2,2)]', [1,3,2])/1000; % left hip JC
RHJC  = permute([C3D.data.RHJC(ev_1:ev_2,1),C3D.data.RHJC(ev_1:ev_2,3),-C3D.data.RHJC(ev_1:ev_2,2)]', [1,3,2])/1000; % right hip JC
KNE   = permute([C3D.data.LKNE(ev_1:ev_2,1),C3D.data.LKNE(ev_1:ev_2,3),-C3D.data.LKNE(ev_1:ev_2,2)]', [1,3,2])/1000; % knee lateral epicondyle
KNI   = permute([C3D.data.LKNI(ev_1:ev_2,1),C3D.data.LKNI(ev_1:ev_2,3),-C3D.data.LKNI(ev_1:ev_2,2)]', [1,3,2])/1000; % knee medial epicondyle
ANK   = permute([C3D.data.LANK(ev_1:ev_2,1),C3D.data.LANK(ev_1:ev_2,3),-C3D.data.LANK(ev_1:ev_2,2)]', [1,3,2])/1000; % ankle lateral malleolus
MED   = permute([C3D.data.LMED(ev_1:ev_2,1),C3D.data.LMED(ev_1:ev_2,3),-C3D.data.LMED(ev_1:ev_2,2)]', [1,3,2])/1000; % ankle medial malleolus
SACR  = permute([C3D.data.SACR(ev_1:ev_2,1),C3D.data.SACR(ev_1:ev_2,3),-C3D.data.SACR(ev_1:ev_2,2)]', [1,3,2])/1000; % medial posterior iliac spine
LTHI  = permute([C3D.data.LTHI(ev_1:ev_2,1),C3D.data.LTHI(ev_1:ev_2,3),-C3D.data.LTHI(ev_1:ev_2,2)]', [1,3,2])/1000; % wond shank
LASIS = permute([C3D.data.LASI(ev_1:ev_2,1),C3D.data.LASI(ev_1:ev_2,3),-C3D.data.LASI(ev_1:ev_2,2)]', [1,3,2])/1000; % left anterior iliac spine
RASIS = permute([C3D.data.RASI(ev_1:ev_2,1),C3D.data.RASI(ev_1:ev_2,3),-C3D.data.RASI(ev_1:ev_2,2)]', [1,3,2])/1000; % right anterior iliac spine
LPSIS = permute([C3D.data.LPSI(ev_1:ev_2,1),C3D.data.LPSI(ev_1:ev_2,3),-C3D.data.LPSI(ev_1:ev_2,2)]', [1,3,2])/1000; % left posterior iliac spine
RPSIS = permute([C3D.data.RPSI(ev_1:ev_2,1),C3D.data.RPSI(ev_1:ev_2,3),-C3D.data.RPSI(ev_1:ev_2,2)]', [1,3,2])/1000; % rigth posterior iliac spine
midHJC = permute([C3D.data.midHJC(ev_1:ev_2,1),C3D.data.midHJC(ev_1:ev_2,3),-C3D.data.midHJC(ev_1:ev_2,2)]', [1,3,2])/1000; % medial point between hip JCs
midASIS = permute([C3D.data.midASIS(ev_1:ev_2,1),C3D.data.midASIS(ev_1:ev_2,3),-C3D.data.midASIS(ev_1:ev_2,2)]', [1,3,2])/1000; % medial point between anterior iliac spines

% Segment parameters (left)
ufoot  = Vnorm_array3(MET2-HEE); % calcaneus to 2nd metatarsal head
wfoot  = Vnorm_array3(MET1-MET5); % 1st to 5th metatarsal head
ushank = Vnorm_array3(cross((TIB-AJC),(KJC-AJC))); % vector normal to plan (KJC, TIB, AJC)
wshank = Vnorm_array3(MED-ANK); % medial to lateral malleolus
uthigh = Vnorm_array3(cross((LHJC-KNE),(KNI-LHJC))); % vector normal to plan (KNE, LHJC, KNI)
wthigh = Vnorm_array3(KNI-KNE);
upelvis = Vnorm_array3(midASIS-SACR);
wpelvis = Vnorm_array3(RHJC-LHJC);

Segment(2).Q = [ufoot; AJC; MET2; wfoot];
Segment(2).rM = [HEE, MET1, MET2, MET5];
Segment(3).Q = [ushank; KJC; AJC; wshank];
Segment(3).rM = [TIB, ANK, MED];
Segment(4).Q = [uthigh; LHJC; KJC; wthigh];
Segment(4).rM = [KNE, KNI, LTHI]; 
Segment(5).Q = [upelvis; SACR; midHJC; wpelvis];
Segment(5).rM = [LASIS, RASIS, LPSIS, RPSIS]; 

%  visualisation
% Main_Segment_Visualisation
% Joint = Joint_Kinematics(Segment)
% Main_Joint_Kinematics_Curves



% cd 'D:\Helical Axis\Test\'
save('Segment.mat', 'Segment')
end