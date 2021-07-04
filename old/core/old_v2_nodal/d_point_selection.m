
function [ind,B,interp_atx] = d_point_selection(int_method, H)
% Build the cross indices
%
switch int_method
    case {'QDEIM'}
        [~,~,e] = qr(H', 'vector');
        ind = e(1:size(H,2));
        interp_atx = H(ind,:);
        B =  H/interp_atx;
    case {'DEIM'}
        disp('Not implemented')
    case {'MaxVol'}
        [ind,B] = maxvol(H);
        interp_atx = H(ind,:);
end
if cond(interp_atx) > 1E5
    disp('warning: poor condition number in the interpolation')
end
end