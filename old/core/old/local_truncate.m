function [B,A,r] = local_truncate(loc_err_tol, min_rank, max_rank, oned, F, rold)
%Truncate the svd for each TT block
%
%Tiangang Cui, August, 2019

switch oned.type
    case{'Lagrange'}
        [U,S,V] = svd( kron(speye(rold), oned.mass_L') * F, 0 );
        s   = diag(S);
        % truncation index r
        ind = s/s(1) > loc_err_tol;
        r   = max(1, min(sum(ind), max_rank));
        % interpolation basis
        B   = full(kron(speye(rold), oned.mass_L') \ U(:,1:r));
    otherwise
        [U,S,V] = svd( kron(speye(rold), oned.node2basis) * F, 0 );
        s   = diag(S);
        % truncation index r
        ind = s/s(1) > loc_err_tol;
        r   = max(min_rank, min(sum(ind), max_rank));
        % basis evaluation multiplied by coefficients
        B   = kron(speye(rold), oned.basis2node)*U(:,1:r);
end

A = s(1:r).*V(:,1:r)';

end