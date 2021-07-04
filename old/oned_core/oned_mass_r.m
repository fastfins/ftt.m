function R = oned_mass_r(oned, interp_w)
%Evaluate the right factor of the one dimensional mass matrix
%
%Tiangang Cui, August, 2019

switch oned.type
    case{'Lagrange'}
        R = oned.mass_L'*interp_w;
    otherwise
        R = interp_w.*sqrt(oned.weights./oned.omegas);
end

end