function ftt = build_ftt_amen(func, ftt, options, debug_x)
%Construct function tensor train using AMEN
%
%Tiangang Cui, August, 2019

% get dimensions of each TT block
[d,rs,~] = ftt_size(ftt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start

f_evals  = 0;
als_iter = 0;
errs     = ones(1,d);

while true % run ALS
    if ftt.direction > 0
        ind = 1:(d-1);
    else
        ind = d:-1:2;
    end
    
    % start
    for k = ind
        if ftt.direction > 0
            if k == 1
                Jx_left = [];
                Jr_left = [];
                Rw_left = 1;
            else
                Jx_left = ftt.interp_x{k-1};
                Jr_left = ftt.res_x{k-1};
                Rw_left = ftt.res_w{k-1};
            end

            % evaluate interpolant function at x_k nodes
            [F, nf] = local_block(options.ng_flag, ftt.oneds{k}, Jx_left, ftt.interp_x{k+1}, func);
            % update relative error and increase counter
            errs(k) = local_error(ftt.cores{k}, F);
            % evaluate res function at x_k nodes
            [Fr,nr] = local_block(options.ng_flag, ftt.oneds{k}, Jr_left, ftt.res_x{k+1}, func);
            f_evals = f_evals + nf + nr;
            
            % evaluate update function at x_k nodes
            if k > 1
                [Fu,nu] = local_block(options.ng_flag, ftt.oneds{k}, Jx_left, ftt.res_x{k+1}, func);
                f_evals = f_evals + nu;
            else
                Fu = Fr;
            end
            
            [ftt.cores{k}, ftt.interp_x{k}, ftt.res_w{k}, ftt.res_x{k}, ftt.cores{k+1}] = build_basis(options, ...
                ftt.oneds{k}, ftt.direction, Jx_left, Jr_left, Rw_left, ftt.res_w{k+1}, ftt.cores{k+1}, F, Fu, Fr);
            rs(k) = size(ftt.cores{k}, 3);
            
        else
            if k == d
                Jx_right = [];
                Jr_right = [];
                Rw_right = 1;
            else
                Jx_right = ftt.interp_x{k+1};
                Jr_right = ftt.res_x{k+1};
                Rw_right = ftt.res_w{k+1};
            end
            
            % evaluate interpolant function at x_k nodes
            [F, nf] = local_block(options.ng_flag, ftt.oneds{k}, ftt.interp_x{k-1}, Jx_right, func);
            % update relative error and increase counter
            
            errs(k) = local_error(ftt.cores{k}, F);
            % evaluate res function at x_k nodes
            [Fr,nr] = local_block(options.ng_flag, ftt.oneds{k}, ftt.res_x{k-1}, Jr_right, func);
            f_evals = f_evals + nf + nr;
            
            % different from left iteration, k + 1 is the previous index
            % evaluate update function at x_k nodes
            if k < d
                [Fu,nu] = local_block(options.ng_flag, ftt.oneds{k}, ftt.res_x{k-1}, Jx_right, func);
                f_evals = f_evals + nu;
            else
                Fu = Fr;
            end
            
            [ftt.cores{k}, ftt.interp_x{k}, ftt.res_w{k}, ftt.res_x{k}, ftt.cores{k-1}] = build_basis(options,... 
                ftt.oneds{k}, ftt.direction, Jx_right, Jr_right, ftt.res_w{k-1}, Rw_right, ftt.cores{k-1}, F, Fu, Fr);
            rs(k-1) = size(ftt.cores{k}, 1);
        end
    end
    
    als_iter = als_iter + 1;
    % evaluate the ftt and give error estimates
    debug_errs = [];
    if ~isempty(debug_x)
        exact   = func(debug_x);
        approx  = eval_ftt(ftt, debug_x);
        debug_errs(1) = max(abs(exact(:) - approx(:)))  / max(abs(exact(:)));
        debug_errs(2) = mean(abs(exact(:) - approx(:))) / max(abs(exact(:)));
    end
    print_iteration(als_iter, errs, rs, f_evals, debug_errs);
    
    if als_iter == options.max_als || max(errs) < options.err_tol
        disp('ALS completed')
        disp(rs)
        break;
    else
        % flip direction
        ftt.direction = -ftt.direction;
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [core, interp_x, res_w, res_x, core_next] = build_basis(options, oned, direction, ...
    interp_xold, res_xold, res_w_l, res_w_r, core_next, F, Fu, Fr)
% find the eigen function of the Schmidt operator
% cross case, nileft x niright kernels,

nbleft  = size(F, 1);
nbright = size(F, 3);
nrleft  = size(Fr,1);
nrright = size(Fr,3);
nnodes  = oned.num_nodes;
m       = size(F, 4);

% dimension of the next core
rn1 = size(core_next, 1);
nn  = size(core_next, 2);
rn2 = size(core_next, 3);


if direction > 0
    % from left >k is integrated
    % T*T' give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F,  nbleft, nnodes, nbright, []), [2,1,3,4]), nnodes*nbleft, []);
    Fu      = reshape(permute( reshape(Fu, nbleft, nnodes, nrright, []), [2,1,3,4]), nnodes*nbleft, []);
    rold    = nbleft;
    rnext   = rn1;
else    
    % <k is integrated
    % T'*T give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F,  nbleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []);
    Fu      = reshape(permute( reshape(Fu, nrleft, nnodes, nbright, []), [2,3,1,4]), nnodes*nbright, []);
    rold    = nbright;
    rnext   = rn2;
end

[B,A,r] = local_truncate(options.loc_err_tol, 1, max(options.max_rank-options.kick_rank,1), oned, F, rold);

% compute ress
% BSV' = FT is arranged as nodes, rold (fixed), rnew (to be updated)
% rnew dimension needs to be projected onto res basis from the update
% direction, rold dimension needs to be projected onto res basis from
% the fixed dimension
% for the right projection

if direction > 0
    % A is r x rnext x m, multiply the rnext dim with res_w_r
    tmp_r   = reshape(permute(reshape(A,r,[],m), [1,3,2]), [], rn1)*res_w_r;
    tmp_r   = reshape(permute(reshape(tmp_r,r,m,[]),[1,3,2]),r,[]);
    Fu      = Fu - B*tmp_r;
    
    % for the left projection
    tmp_l   = res_w_l*reshape(permute(reshape(B,nnodes,rold,r), [2,1,3]), rold, []);
    tmp_l   = reshape(tmp_l, nrleft*nnodes, r);
    
    % align Fr as rold (nrleft), nodes, rnew (nrright), m
    Fr      = reshape(Fr,nrleft,nnodes,nrright,[]) - reshape(tmp_l*tmp_r,nrleft,nnodes,nrright,[]);
    Fr      = reshape(permute(Fr,[2,1,3,4]), nnodes*nrleft, []);
    rrold   = nrleft;
else
    % A is r x rnext x m, multiply the rnext dim with res_w_l
    tmp_lt  = res_w_l*reshape(permute(reshape(A,r,[],m), [2,1,3]), rn2, []);
    % r x nr_l x m
    tmp_lt  = reshape(permute(reshape(tmp_lt,[],r,m),[2,1,3]),r,[]);
    % Fu is aligned with F
    Fu      = Fu - B*tmp_lt;
    
    % for the right projection
    tmp_r   = reshape(permute(reshape(B,nnodes,rold,r), [3,1,2]), [], rold)*res_w_r; % r x n x nr_r
    tmp_r   = reshape(tmp_r, r, nnodes*nrright);
    
    % align Fr as rnew (nrleft), nodes, rold (nrright), m
    Fr      = reshape(Fr,nrleft,nnodes,nrright,m) - permute(reshape(tmp_lt'*tmp_r,nrleft,m,nnodes,nrright), [1,3,4,2]);
    Fr      = reshape(permute(Fr,[2,3,1,4]), nnodes*nrright, []);
    rrold   = nrright;
end

% enrich basis

switch oned.type
    case{'Lagrange'}
        T     = kron(speye(rold), oned.mass_L') * [B, Fu];
        [Q,R] = qr(T,0);
        if size(T,1) < size(T,2)
            Q = Q(:,1:size(T,1));
        end
        % interpolation basis
        B     = full(kron(speye(rold),oned.mass_L') \ Q);
    otherwise
        T     = kron(speye(rold), oned.node2basis) * [B, Fu];
        [Q,R] = qr(T,0);
        if size(T,1) < size(T,2)
            Q = Q(:,1:size(T,1));
        end
        % basis evaluation multiplied by coefficients
        B     = full(kron(speye(rold), oned.basis2node)*Q);
end

%{
rnew     = size(B,2);
indf     = point_selection(options.int_method, B);
interp_x = local_index(oned, direction, interp_xold, indf);
    
Qr      = local_truncate(options.loc_err_tol*1E-5, 1, options.kick_rank, oned, Fr, rrold);
% rold is the leading dimension in the point selection
indr    = point_selection(options.int_method, Qr);
res_x   = local_index(oned, direction, res_xold, indr);

interp_atx  = B(indf,:);
if cond(interp_atx) > 1E5
    disp('warning: poor condition number in the interpolation')
end
core    = B/interp_atx;
%}

rnew        = size(B,2);
[indf,core] = point_selection(options.int_method, B);
interp_x    = local_index(oned, direction, interp_xold, indf);

interp_atx  = B(indf,:);
if cond(interp_atx) > 1E5
    disp('warning: poor condition number in the interpolation')
end

Qr      = local_truncate(options.loc_err_tol*eps, 1, options.kick_rank, oned, Fr, rrold);
% rold is the leading dimension in the point selection
indr    = point_selection(options.int_method, Qr);
res_x   = local_index(oned, direction, res_xold, indr);

couple  = reshape(interp_atx(:,1:r)*(R(1:r,1:r)*A), rnew, rnext, m);
if direction > 0
    % from left >k is integrated
    % get the index for transpose A(:,:,j) for each j
    core = permute(reshape(core,nnodes,rold,rnew), [2,1,3]);
    %the right factor is r x (rn1+enrich) x m, only push the r x rn1 x m 
    %block to the next block, first permute the block to m x r x rn1
    core_next = reshape(permute(couple, [3,1,2]), [], rnext) * reshape(core_next, rnext, []);
    core_next = permute(reshape(core_next, m, rnew, nn, rn2), [2, 3, 4, 1]);
    
    tmp     = reshape(permute(reshape(res_w_l*reshape(core, rold, nnodes*rnew),rrold,nnodes,rnew),[2,1,3]),[],rnew);
    res_w   = tmp(indr,:);
else
    % <k is integrated
    % unfold the 3D tensor
    core    = permute(reshape(core,nnodes,rold,rnew), [3,1,2]);
    %the right factor is r x (rn2+enrich) x m, only push the r x rn2 x m 
    %block to the next block, first permute the block to rn2 x r x m
    core_next = reshape(core_next, [], rnext) * reshape(permute(couple, [2,1,3]), rnext, []);
    core_next = reshape(core_next, rn1, nn, rnew, m);
    
    tmp     = reshape(reshape(reshape(core, rnew*nnodes, rold)*res_w_r,rnew,nnodes,rrold),rnew,[]);
    res_w   = tmp(:,indr);
end

end