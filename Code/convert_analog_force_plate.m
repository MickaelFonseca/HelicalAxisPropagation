function C3D = convert_analog_force_plate(C3D)
% to complete
Name1={'Fx';'Fy';'Fz';'Mx';'My';'Mz';'Fx1';'Fy1';'Fz1';'Mx1';'My1';'Mz1';};
Name2={'Fx1';'Fy1';'Fz1';'Mx1';'My1';'Mz1';'Fx2';'Fy2';'Fz2';'Mx2';'My2';'Mz2';};
if isfield(C3D,'analog')
    if isfield(C3D.analog,'Fx2')==0
        if isfield(C3D.analog,'Fx')==1
            temp=struct();
            for i=1:length(Name1)
                temp.(Name2{i})=C3D.analog.(Name1{i});
            end
            C3D.analog = rmfield(C3D.analog, Name1);
            for i=1:length(Name1)
                C3D.analog.(Name2{i})=temp.(Name2{i});
            end
        end
        
    end
end
