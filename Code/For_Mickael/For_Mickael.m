[a,b] = uigetfile('.c3d');
C3D_filenames = [b,a];
acq = btkReadAcquisition(C3D_filenames);
C3D.data = btkGetMarkers(acq);
nf = length(C3D.data.LHEE);

% XYZ (Qualisys) => XZ(-Y) (ISB) in m
% Left leg
HEE  = permute([C3D.data.LHEE(:,1),C3D.data.LHEE(:,3),-C3D.data.LHEE(:,2)]', [1,3,2])/1000; % heel
MET1 = permute([C3D.data.LFMH(:,1),C3D.data.LFMH(:,3),-C3D.data.LFMH(:,2)]', [1,3,2])/1000; % 1st metatarsal head
MET2 = permute([C3D.data.LSMH(:,1),C3D.data.LSMH(:,3),-C3D.data.LSMH(:,2)]', [1,3,2])/1000; % 2nd metatarsal head
MET5 = permute([C3D.data.LVMH(:,1),C3D.data.LVMH(:,3),-C3D.data.LVMH(:,2)]', [1,3,2])/1000; % 5th metatarsal head
AJC  = permute([C3D.data.LAJC(:,1),C3D.data.LAJC(:,3),-C3D.data.LAJC(:,2)]', [1,3,2])/1000; % ankle JC
TIB  = permute([C3D.data.LTIB(:,1),C3D.data.LTIB(:,3),-C3D.data.LTIB(:,2)]', [1,3,2])/1000; % wond thigh
KJC  = permute([C3D.data.LKJC(:,1),C3D.data.LKJC(:,3),-C3D.data.LKJC(:,2)]', [1,3,2])/1000; % knee JC
LHJC = permute([C3D.data.LHJC(:,1),C3D.data.LHJC(:,3),-C3D.data.LHJC(:,2)]', [1,3,2])/1000; % left hip JC
RHJC = permute([C3D.data.RHJC(:,1),C3D.data.RHJC(:,3),-C3D.data.RHJC(:,2)]', [1,3,2])/1000; % right hip JC
KNE  = permute([C3D.data.LKNE(:,1),C3D.data.LKNE(:,3),-C3D.data.LKNE(:,2)]', [1,3,2])/1000; % knee lateral epicondyle
KNI  = permute([C3D.data.LKNI(:,1),C3D.data.LKNI(:,3),-C3D.data.LKNI(:,2)]', [1,3,2])/1000; % knee medial epicondyle
ANK  = permute([C3D.data.LANK(:,1),C3D.data.LANK(:,3),-C3D.data.LANK(:,2)]', [1,3,2])/1000; % ankle lateral malleolus
MED  = permute([C3D.data.LMED(:,1),C3D.data.LMED(:,3),-C3D.data.LMED(:,2)]', [1,3,2])/1000; % ankle medial malleolus
SACR = permute([C3D.data.SACR(:,1),C3D.data.SACR(:,3),-C3D.data.SACR(:,2)]', [1,3,2])/1000; % medial posterior iliac spine
LTHI = permute([C3D.data.LTHI(:,1),C3D.data.LTHI(:,3),-C3D.data.LTHI(:,2)]', [1,3,2])/1000; % wond shank
LASIS = permute([C3D.data.LASI(:,1),C3D.data.LASI(:,3),-C3D.data.LASI(:,2)]', [1,3,2])/1000; % left anterior iliac spine
RASIS = permute([C3D.data.RASI(:,1),C3D.data.RASI(:,3),-C3D.data.RASI(:,2)]', [1,3,2])/1000; % right anterior iliac spine
LPSIS = permute([C3D.data.LPSI(:,1),C3D.data.LPSI(:,3),-C3D.data.LPSI(:,2)]', [1,3,2])/1000; % left posterior iliac spine
RPSIS = permute([C3D.data.RPSI(:,1),C3D.data.RPSI(:,3),-C3D.data.RPSI(:,2)]', [1,3,2])/1000; % rigth posterior iliac spine
midHJC = permute([C3D.data.midHJC(:,1),C3D.data.midHJC(:,3),-C3D.data.midHJC(:,2)]', [1,3,2])/1000; % medial point between hip JCs
midASIS = permute([C3D.data.midASIS(:,1),C3D.data.midASIS(:,3),-C3D.data.midASIS(:,2)]', [1,3,2])/1000; % medial point between anterior iliac spines

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

% Some visualisation
Main_Segment_Visualisation
Joint = Joint_Kinematics(Segment)
Main_Joint_Kinematics_Curves

cd 'D:\Helical Axis\Test\'
save('Segment.mat', 'Segment')

