function rerr = local_error(f_old, f)
%Compute local error for tt block
%
%Tiangang Cui, August, 2019

diff = f_old(:) - f(:);
rerr = max(abs(diff)) / max(abs(f(:))); 
 
end