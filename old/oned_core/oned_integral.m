function f_int = oned_integral(oned, interp_w)
%Evaluate one dimensional integrals
%
%Tiangang Cui, August, 2019

switch oned.type
    case{'Lagrange'}
        % f_int = sum( repmat(oned.weights(:), 1, size(interp_w, 2))*interp_w, 1 );
        f_int = (oned.weights(:)')*interp_w;
    otherwise
        f_int = sum( (oned.weights./oned.omegas).*interp_w, 1 );
end

end
