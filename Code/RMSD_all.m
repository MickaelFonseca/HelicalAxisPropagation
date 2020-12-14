function [RMSD, m_RMSD] = RMSD_all(Res,C3D_filenames, Joint_name)
% author: M.Fonseca, November 2019

if Joint_name == string('Ankle')
    J = 'ank';
    Jang = 'A';
elseif Joint_name == string('Knee')
    J = 'knee';
    Jang = 'K';
elseif Joint_name == string('Hip')
    J = 'hip';
    Jang = 'H';
end

mean_Theta = mean(Res.(strcat('THETA_',J)),2);
mean_ALPHA = mean(Res.(strcat('ALPHA_',J)),2);
mean_Euler_X = mean(Res.(strcat(Jang,'F_Eul')),2);
mean_Euler_Y = mean(Res.(strcat(Jang,'A_Eul')),2);
mean_Euler_Z = mean(Res.(strcat(Jang,'R_Eul')),2);
mean_KT_X = mean(Res.(strcat(Jang,'F_KT')),2);
mean_KT_Y = mean(Res.(strcat(Jang,'A_KT')),2);
mean_KT_Z = mean(Res.(strcat(Jang,'R_KT')),2);

for i = 1:length(C3D_filenames)
    %% RMSD
    RMSD.HelicalAxis(i) = sqrt(mean((Res.(strcat('THETA_',J))(:,i)- mean_Theta).^2));
    RMSD.HelicalAxisVAR(i) = sqrt(mean((Res.(strcat('ALPHA_',J))(:,i)- mean_ALPHA).^2));
    RMSD.Euler_x(i) = sqrt(mean((Res.(strcat(Jang,'F_Eul'))(:,i)- mean_Euler_X).^2));
    RMSD.Euler_y(i) = sqrt(mean((Res.(strcat(Jang,'A_Eul'))(:,i)- mean_Euler_Y).^2));    
    RMSD.Euler_z(i) = sqrt(mean((Res.(strcat(Jang,'R_Eul'))(:,i)- mean_Euler_Z).^2));    
    RMSD.KT_x(i) = sqrt(mean((Res.(strcat(Jang,'F_KT'))(:,i)- mean_KT_X).^2));
    RMSD.KT_y(i) = sqrt(mean((Res.(strcat(Jang,'A_KT'))(:,i)- mean_KT_Y).^2));    
    RMSD.KT_z(i) = sqrt(mean((Res.(strcat(Jang,'R_KT'))(:,i)- mean_KT_Z).^2));        
end

m_RMSD = table(mean(RMSD.HelicalAxis), mean(RMSD.HelicalAxisVAR), mean(RMSD.Euler_x), mean(RMSD.KT_x),mean(RMSD.Euler_y), mean(RMSD.KT_y), mean(RMSD.Euler_z), mean(RMSD.KT_z));
m_RMSD.Properties.VariableNames = {'Helical_Angle', 'Variation_Helical_Axis', 'Euler_x', 'Projected_KTheta_x','Euler_y', 'Projected_KTheta_y',  'Euler_z', 'Projected_KTheta_z'};

end