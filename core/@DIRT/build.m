function obj = build(obj, func, d, beta, pol, opt)

% this can be QMC sets from i.i.d. Gaussian
if obj.qmc_flag
    l = ceil(log2(opts.nsamples));
    opts.nsamples = 2^l;
    zs=qmcnodes(opts.np,l); 
    samples = erfinv(2*zs-1)*sqrt(2);
else
    samples = randn(opts.np, opts.nsamples);
end

beta = 0;
irtprev = 4;
for iter = 1:opts.max_iter
    beta_p = beta;
    disp(['iteration ' num2str(iter)])
    
    [us,lfs,Jr_uv] = eval_dirt(liss, irts, iter-1, samples, grids);
    
    mllkds  = zeros(1, opts.nsamples);
    gmllkds = zeros(opts.np, opts.nsamples);
    for i = 1:opts.nsamples
        [~,mllkds(i),~,gmllkds(:,i)] = opts.mlpost(us(:,i));
    end
    %
    data_sets{iter,1} = mllkds;
    data_sets{iter,2} = gmllkds;
    
    %
    % adapt on the temperature, increase temperature to 1
    beta = max(beta_p,opts.min_temp);
    % compute ess over sample size
    ess = ess_ratio(beta_p, beta, mllkds);
    while ess > opts.ess_tol
        beta = beta*opts.beta_fac;
        ess = ess_ratio(beta_p, beta, mllkds);
    end
    beta = min(1, beta);
    ess = ess_ratio(beta_p, beta, mllkds);
    betas(iter) = beta;
    % log of normalising const between two adjacent layers
    log_weight = (beta_p - beta)*mllkds;
    ref_weight = max(log_weight);
    %
    disp(['    temperature ' num2str(beta) ', ess ' num2str(ess)])
        
    if iter > 1
        if opts.debug
            lwd = -beta_p*mllkds-lfs;
            rwd = max(lwd);
            ess2 = ( sum(exp(lwd-rwd)).^2/sum(exp(2*(lwd-rwd))) ) / opts.nsamples;
            %
            disp(['N/ess ' num2str(1/ess2) ', ess ' num2str(ess2)])
        end
        log_weight = -beta*mllkds-lfs;
        ref_weight = max(log_weight);
    end
    
    % build LIS
    g = (beta_p - beta)*gmllkds;
    g = backtracking(liss, iter-1, Jr_uv, g);
    weights = exp(log_weight-ref_weight)/sum(exp(log_weight-ref_weight));
    [V,S,~] = svd( g.*sqrt(weights(:)'), 'econ' );
    s = diag(S);
    figure(1);
    semilogy(s);
    drawnow;
    cs = cumsum(s.^2);
    err = 0.5*sqrt(cs(end)-cs);
    r = max(sum(err>opts.tru_tol), opts.min_rank);
    r = opts.d;
    liss{iter}.V = V(:,1:r);
    liss{iter}.s = s;
    liss{iter}.r = r;
    %liss{k}.basis = liss{k}.V.*liss{k}.d(:)';
    %liss{k}.shift = 0;
    
    disp(['    LIS dimension ' num2str(r)])
    
    % Ratio function for current iteration
    f_sub = @(vr) ratio_fun(opts.mlpost, beta_p, beta, liss, irts, iter, liss{iter}.V, vr', grids)';
        
    if (strcmpi(opts.method, 'ftt'))
        % resampling
        ind = datasample(1:opts.nsamples, opts.nsamples, 'weights', weights);
        v_sub = liss{iter}.V'*samples(:,ind);
        % Build FTT
        ftts{iter} = build_ftt(f_sub, r, [], opt, 'sample_x', v_sub);
        % build IRT
        irts{iter} = build_irt(ftts{iter});
    else
        % build TT
        if (iter==1)
            grids{iter} = repmat({x0}, r, 1);
        else
            grids{iter} = repmat({xi}, r, 1);
        end
        ttgrids = tt_meshgrid_vert(cellfun(@(x)tt_tensor(x), grids{iter}, 'uni', 0));
        irts{iter} = amen_cross_s(ttgrids, f_sub, 0, 'nswp', 1, 'kickrank', 0, 'y0', irtprev);
        irtprev = irts{iter};
        figure(2);
        if (iter==1)
            surf(x0, x0, full(dot(tt_ones(numel(x0), r-2), irts{iter}, 3, r), numel(x0)*[1 1]), 'EdgeColor', 'none')
        else
            surf(xi, xi, full(dot(tt_ones(numel(xi), r-2), irts{iter}, 3, r), numel(xi)*[1 1]), 'EdgeColor', 'none')
        end
        shading interp;
        drawnow;
        irts{iter} = core2cell(irts{iter});
    end
    
    % stop
    if abs(beta - 1) < 1E-10
        break;
    end
end

if opts.debug && iter > 1
    [us,lfs] = eval_dirt(liss, irts , iter, samples, grids);
    mllkds = zeros(1, opts.nsamples);
    for i = 1:opts.nsamples
        [~,mllkds(i)] = opts.mlpost(us(:,i));
    end
    lwd = -mllkds-lfs;
    rwd = max(lwd);
    ess2 = ( sum(exp(lwd-rwd)).^2/sum(exp(2*(lwd-rwd))) ) / opts.nsamples;
    %
    disp(['N/ess ' num2str(1/ess2) ', ess ' num2str(ess2)])
end


end