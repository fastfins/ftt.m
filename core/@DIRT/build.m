function obj = build(obj, func)

% reference samples
%{
if obj.qmc_flag
    l = ceil(log2(obj.nsamples));
    opts.nsamples = 2^l;
    zs = qmcnodes(opts.np,l);
else
    zs = rand(obj.d, obj.nsamples);
end
%}

zs = rand(obj.d, obj.nsamples);
samples = eval_icdf(obj.diag, zs);

beta = 0;
obj.n_layers = 0;
while obj.n_layers < obj.max_iter
    %
    % evaluate the map and the target function
    [x,logfx] = eval_irt(obj, samples);
    [mllkds, mlps] = func(x);
    % determine the new temperature
    if obj.adapt_tempering
        beta_p = beta;
        beta = max(beta_p,obj.min_temperature);
        switch obj.method
            case {'Aratio'}
                % compute ess over sample size
                ess = ess_ratio((beta_p-beta)*mllkds);
                while ess > ess_tol
                    beta = beta*obj.tempering_factor;
                    ess = ess_ratio((beta_p-beta)*mllkds);
                end
                beta = min(1, beta);
                ess = ess_ratio((beta_p-beta)*mllkds);
            case {'Eratio'}
                % compute ess over sample size
                ess = ess_ratio(-beta*mllkds-mlps-logfx);
                while ess > ess_tol
                    beta = beta*obj.tempering_factor;
                    ess = ess_ratio(-beta*mllkds-mlps-logfx);
                end
                beta = min(1, beta);
                ess = ess_ratio(-beta*mllkds-mlps-logfx);
        end
        obj.temperatures(obj.n_layers) = beta;
    else
        beta_p = beta;
        beta = obj.temperatures(obj.n_layers);
    end
    
    % squared Hellinger error between logfx and (-beta_p*mllkds-mlps)
    if obj.n_layers > 1
        [~,dh2,~] = f_divergence(logfx, -beta_p*mllkds-mlps);
    else
        dh2 = 1;
    end
    %
    fprintf('iteration=%2d, Hellinger2 error=%3.3e, ess ratio=%3.3e cum#fevals=%3.3e\n', ...
        obj.n_layers-1, dh2, ess, f_eval);
    
    log_weights = -beta*mllkds-mlps-logfx;
    log_weights = log_weights - max(log_weights);
    %
    ind = datasample(1:obj.nsamples, obj.nsamples, 'weights', exp(log_weights));
    %
    % Ratio function for current iteration
    newf = @(z) ratio_fun(obj, func, beta_p, beta, z);
    %
    if obj.debug_flag
        obj.irts{obj.n_layers+1} =  SIRT(newf, d, pol, opt, 'debug_x', dx, 'sample_x', samples(:,ind));
    else
        obj.irts{obj.n_layers+1} =  SIRT(newf, d, pol, opt, 'sample_x', samples(:,ind));
    end
    obj.n_layers = obj.n_layers + 1;
    
    % stop
    if abs(beta - 1) < 1E-10
        break;
    end
end

end