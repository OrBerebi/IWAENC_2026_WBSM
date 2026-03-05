function [h, parametric, roomParams] = calc_rir(fs, roomDim, sourcePos, arrayPos, R, calc_parametric_rir_name_value_pairs, rir_from_parametric_name_value_args,Anechoic)
% for details about name-value argument see
% image_method.calc_parametric_rir.m 
% and
% rir_from_parametric.m.

% Author: Tom Shlomo, ACLab BGU, 2020

% calc paramteric info
parametric = image_method.calc_parametric_rir(roomDim, sourcePos, arrayPos,R, calc_parametric_rir_name_value_pairs{:});
if Anechoic
    %take the first reflaction and discard the rest
    parametric = parametric(1,:);
end
% convert to rir signal
[h, parametric.delay, roomParams] = image_method.rir_from_parametric(fs, parametric.delay, parametric.amp, parametric.omega, rir_from_parametric_name_value_args{:});
roomParams.rd = RoomParams.critical_distance_diffuse(roomDim, R);
end