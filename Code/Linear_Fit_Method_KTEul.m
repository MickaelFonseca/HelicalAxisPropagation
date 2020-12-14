function [a1, a0, R2, Ya]= Linear_Fit_Method_KTEul(Res, angle_abrev, C3D_filenames)
% author: M.Fonseca, November 2019
% based in Iosa et al. 2014

% Function applies the linear fit method
% Pa: dataset under investigation
% Pref: dataset reference (mean of overall results)
% a1: angular coefficient (mean variation of Pa w.r.t. Pref)
% a0: intercept of fitting line (predicts the value of Pa when Pref = 0 (shift))
% R2: Measure of the strength of the linear relationship between Pa & Pref
% Ya: Linear function which approximates Pa by mean of Linear
%     transformation of values of Pref

Pref = mean(Res.(angle_abrev),2);
for i = 1:length(C3D_filenames)
    Pa = Res.(angle_abrev)(:,i);
    for n = 1:length(Pref)
        a1_1(n) = ((Pref(n)-mean(Pref))*(Pa(n)-mean(Pa)));
        a1_2(n) = (Pref(n)-mean(Pref))^2;
    end
    a1(i) = (sum(a1_1)/sum(a1_2));
    a0(i) = mean(Pa)-a1(i)*(mean(Pref));
    
    for n = 1:length(Pref)
        R2_1(n) = (a0(i) + a1(i)*Pref(n)-mean(Pa));
        R2_2(n) = (Pa(n) - mean(Pa))^2;
    end
    R2(i) = sum((R2_1)/(R2_2));
    Ya(:,i) = a1(i)*Pref +a0(i);
end
% a1 = mean(a1);
% a0 = mean(a0);
% R2 = mean(R2);
end