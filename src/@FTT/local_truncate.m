function [B,A,r] = local_truncate(loc_err_tol, min_rank, max_rank, oned, F)
% Truncate the svd for each TT block
%
if isa(oned, 'piecewise')
    %[U,S,V] = svd( kron(speye(rold), obj.oneds{k}.mass_L') * F, 0 );
    [U,S,V] = svd( reshape(oned.mass_R*reshape(F,oned.num_nodes,[]),size(F,1),[]), 0);
    s = diag(S);
    % truncation index r
    ind = s/s(1) > loc_err_tol;
    r = max(min_rank, min(sum(ind), max_rank));
    % interpolation basis
    %B   = full(kron(speye(rold), obj.oneds{k}.mass_L') \ U(:,1:r));
    B = reshape(oned.mass_R\reshape(U(:,1:r),oned.num_nodes,[]),size(F,1),[]);
    A = s(1:r).*V(:,1:r)';
else
    [U,S,V] = svd(F, 0);
    s = diag(S);
    % truncation index r
    ind = s/s(1) > loc_err_tol;
    r = max(min_rank, min(sum(ind), max_rank));
    B = U(:,1:r);
    A = s(1:r).*V(:,1:r)';
end
end