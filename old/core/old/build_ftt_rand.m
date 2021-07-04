function ftt = build_ftt_rand(func, ftt, options, sample_x, debug_x)
%Construct function tensor train using random enrichment
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
    % at the head, update the random enrichment set
    if size(sample_x, 2) < options.kick_rank
        sample_x = [];        
        enrich  = zeros(d, options.kick_rank);
        for k = 1:d
            enrich(k,:) = sample_oned_domain(ftt.oneds{k}, options.kick_rank);
        end
    else
        enrich = sample_x(:,1:options.kick_rank);
        sample_x(:,1:options.kick_rank) = [];
    end

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
            else
                Jx_left = ftt.interp_x{k-1};
            end
            
            % evaluate interpolant function at x_k nodes
            [F, nf] = local_block(ftt.ng_flag, ftt.oneds{k}, Jx_left, ftt.interp_x{k+1}, func);
            [Fe,ne] = local_block(ftt.ng_flag, ftt.oneds{k}, Jx_left, enrich(k+1:d,:),   func);
            % update relative error and increase counter
            errs(k) = local_error(ftt.cores{k}, F);
            f_evals = f_evals + nf + ne;
            
            % Schmidt operator and interpolation, k - 1 is the previous index
            % push couple matrix to the right
            [ftt.cores{k}, ftt.interp_x{k}, ftt.cores{k+1}] = build_basis(options, ...
                ftt.oneds{k}, ftt.direction, Jx_left, ftt.cores{k+1}, cat(3,F,Fe));
            rs(k) = size(ftt.cores{k}, 3);
            
        else
            if k == d
                Jx_right = [];
            else
                Jx_right = ftt.interp_x{k+1};
            end
            % evaluate interpolant function at x_k nodes
            [F, nf] = local_block(ftt.ng_flag, ftt.oneds{k}, ftt.interp_x{k-1}, Jx_right, func);
            % different from left iteration, k + 1 is the previous index
            [Fe,ne] = local_block(ftt.ng_flag, ftt.oneds{k}, enrich(1:k-1,:),   Jx_right, func);
            % update relative error and increase counter
            errs(k) = local_error(ftt.cores{k}, F);
            f_evals = f_evals + nf + ne;
            
            % Schmidt operator and interpolation, k + 1 is the previous index
            % push couple matrix to the right
            [ftt.cores{k}, ftt.interp_x{k}, ftt.cores{k-1}] = build_basis(options, ...
                ftt.oneds{k}, ftt.direction, Jx_right, ftt.cores{k-1}, cat(1,F,Fe));
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

function [core, interp_x, core_next] = build_basis(options, oned, direction, interp_xold, core_next, F)
% find the eigen function of the Schmidt operator
% cross case, nileft x niright kernels,

% dimension of the refinement block
nbleft  = size(F, 1);
nbright = size(F, 3);
n_nodes = oned.num_nodes;
m       = size(F, 4);

% dimension of the next core
rn1 = size(core_next, 1);
nn  = size(core_next, 2);
rn2 = size(core_next, 3);

if direction > 0
    % from left >k is integrated
    % T*T' give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F, nbleft, n_nodes, nbright, []), [2,1,3,4]), n_nodes*nbleft,  []); %(nright+enrich)*m
    rold    = nbleft;
    %rnext   = rn1;
else
    % <k is integrated
    % T'*T give the schmidt operator evaluated at the quadrature points
    F       = reshape(permute( reshape(F, nbleft, n_nodes, nbright, []), [2,3,1,4]), n_nodes*nbright, []); %(nbleft+enrich)*m
    rold    = nbright;
    %rnext   = rn2;
end


[B,A,r]     = local_truncate(options.loc_err_tol, 1, options.max_rank, oned, F, rold);
[ind,core]  = point_selection(options.int_method, B);

interp_x    = local_index(oned, direction, interp_xold, ind);
interp_atx  = B(ind,:);
if cond(interp_atx) > 1E5
    disp('warning: poor condition number in the interpolation')
end

tmp     = reshape(interp_atx*A, r, [], m);
if direction > 0
    % from left >k is integrated
    % get the index for transpose A(:,:,j) for each j
    core = permute(reshape(core,n_nodes,rold,r), [2,1,3]);
    %the right factor is r x (rn1+enrich) x m, only push the r x rn1 x m 
    %block to the next block, first permute the block to m x r x rn1
    core_next = reshape(permute(tmp(:,1:rn1,:), [3,1,2]), [], rn1) * reshape(core_next, rn1, []);
    core_next = permute(reshape(core_next, m, r, nn, rn2), [2, 3, 4, 1]);
else
    % <k is integrated
    % unfold the 3D tensor
    core    = permute(reshape(core,n_nodes,rold,r), [3,1,2]);
    %the right factor is r x (rn2+enrich) x m, only push the r x rn2 x m 
    %block to the next block, first permute the block to rn2 x r x m
    core_next = reshape(core_next, [], rn2) * reshape(permute(tmp(:,1:rn2,:), [2,1,3]), rn2, []);
    core_next = reshape(core_next, rn1, nn, r, m);
end

end
