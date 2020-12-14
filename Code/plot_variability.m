function plot_variability (Res, Joint_name)
        
        if Joint_name == 'Ankle'
            J = 'Ank';
            Jang = 'A';            
        elseif Joint_name == 'Knee'
            J = 'Knee';
            Jang = 'K';
        elseif Joint_name == 'Hip'
            J = 'Hip';
            Jang = 'H';
        end

        subplot(2,3,1)
%         suptitle(Joint_name)
        plot(Res.(strcat('m_T_',J)), 'k')
        hold on
        corridor(Res.(strcat('m_T_',J)),Res.(strcat('sd_T_',J)), 'k', [0:100])
        xlim([0 100])
        title('Helical Angle (theta)')
        
        subplot(2,3,2)
        plot(Res.(strcat('m_A_',J)), 'k')
        hold on
        corridor(Res.(strcat('m_A_',J)),Res.(strcat('sd_A_',J)), 'k', [0:100])
        xlim([0 100])
        title('Variation of Helical Axis w.r.t. Mean Axis')
        
        subplot(2,3,4)
        plot(Res.(strcat('m_',Jang,'F_KT')),'r')
        hold on
        corridor(Res.(strcat('m_',Jang,'F_KT')), Res.(strcat('sd_',Jang,'F_KT')),'r', [0:100])
        plot(Res.(strcat('m_',Jang,'F_Eul')),'b')
        corridor(Res.(strcat('m_',Jang,'F_Eul')), Res.(strcat('sd_',Jang,'F_Eul')),'b', [0:100])
        xlim([0 100])
        title ('Flexion-Extension')
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')
        
        subplot(2,3,5)
        plot(Res.(strcat('m_',Jang,'A_KT')),'r')
        hold on
        corridor(Res.(strcat('m_',Jang,'A_KT')), Res.(strcat('sd_',Jang,'A_KT')),'r', [0:100])
        plot(Res.(strcat('m_',Jang,'A_Eul')),'b')
        corridor(Res.(strcat('m_',Jang,'A_Eul')), Res.(strcat('sd_',Jang,'A_Eul')),'b', [0:100])
        xlim([0 100])
        title ('Abduction-Adduction')
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')

        subplot(2,3,6)
        plot(Res.(strcat('m_',Jang,'R_KT')),'r')
        hold on
        corridor(Res.(strcat('m_',Jang,'R_KT')), Res.(strcat('sd_',Jang,'R_KT')),'r', [0:100])
        plot(Res.(strcat('m_',Jang,'R_Eul')),'b')
        corridor(Res.(strcat('m_',Jang,'R_Eul')), Res.(strcat('sd_',Jang,'R_Eul')),'b', [0:100])
        xlim([0 100])
        title ('Internal-External Rotation')
        xlabel('Sampled instant of time')
        ylabel('Angle (in degree)')
        legend({'Projected ktheta', 'Euler angle'})
end