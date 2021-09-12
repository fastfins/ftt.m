function obj = build(obj, func, oneds, sirt_opt)

beta_factor = 1.1;

%{

    if obj.n_layers > 0 % not in the first layer
        beta_p = obj.betas(obj.n_layers);
        beta = obj.betas(obj.n_layers+1);
        switch obj.method
            case {'Aratio'}
                logf =  (beta_p - beta)*mllkds + logfz;
            case {'Eratio'}
                logf =  - beta*mllkds - mlps - logft + logfz;
        end
    else
        beta = obj.betas(obj.n_layers+1);
        logf = - beta*mllkds - mlps;
    end
%}

% reference samples
%{
if obj.qmc_flag
    l = ceil(log2(obj.n_samples));
    opts.nsamples = 2^l;
    zs = qmcnodes(opts.np,l);
else
    zs = rand(obj.d, obj.n_samples);
end
%}

ns = obj.n_samples+obj.n_debugs;
zs = rand(obj.d, ns);
samples = eval_icdf(obj.diag, zs);

beta = 0;
obj.n_evals = 0;
obj.n_layers = 0;
obj.irts = cell(obj.max_layers, 1);
obj.logz = 0;
poly_counter = 1;
while obj.n_layers < obj.max_layers
    %
    % evaluate the map and the target function
    [x,logfx] = eval_irt(obj, samples);
    [mllkds, mlps] = func(x);
    obj.n_evals = obj.n_evals + ns;
    % determine the new temperature
    if obj.adapt_beta
        beta_p = beta;
        beta = max(beta_p,obj.min_beta);
        switch obj.method
            case {'Aratio'}
                % compute ess over sample size
                ess = ess_ratio((beta_p-beta)*mllkds);
                while ess > obj.ess_tol
                    beta = beta*beta_factor;
                    ess = ess_ratio((beta_p-beta)*mllkds);
                end
                beta = min(1, beta/beta_factor);
                ess = ess_ratio((beta_p-beta)*mllkds);
            case {'Eratio'}
                % compute ess over sample size
                ess = ess_ratio(-beta*mllkds-mlps-logfx);
                while ess > obj.ess_tol
                    beta = beta*beta_factor;
                    ess = ess_ratio(-beta*mllkds-mlps-logfx);
                end
                beta = min(1, beta);
                ess = ess_ratio(-beta*mllkds-mlps-logfx);
        end
        obj.betas(obj.n_layers+1) = beta;
    else
        beta_p = beta;
        beta = obj.betas(obj.n_layers+1);
        ess = ess_ratio((beta_p-beta)*mllkds);
    end
    
    % squared Hellinger error between logfx and (-beta_p*mllkds-mlps)
    if obj.n_layers > 0
        [~,dh2,~] = f_divergence(logfx, -beta_p*mllkds-mlps);
    else
        dh2 = 1;
    end
    %
    fprintf('\n\niteration=%2d, Hell err=%3.3e, cum#fevals=%3.3e, next beta=%3.3e, ess ratio=%3.3e \n', ...
        obj.n_layers, sqrt(dh2), obj.n_evals, beta, ess);
    
    log_weights = -beta*mllkds-mlps-logfx;
    log_weights = log_weights - max(log_weights);
    %
    ind = datasample(1:ns, ns, 'weights', exp(log_weights));
    %
    % Ratio function for current iteration
    % if FTT does not use sqrt_flag, then we should compute sqrt of the
    % ratio function here
    newf = @(z) ratio_fun(obj, func, z, ~sirt_opt.sqrt_flag);
    %
    if obj.n_debugs > 0
        ind1 = ind(1:obj.n_samples);
        ind2 = ind((1:obj.n_debugs)+obj.n_samples);
        %obj.irts{obj.n_layers+1} =  SIRT(newf, obj.d, oneds{poly_counter}, sirt_opt, 'debug_x', samples(:,ind2), 'sample_x', samples(:,ind1));
        
        if obj.n_layers > 1
            obj.irts{obj.n_layers+1} =  SIRT(newf, obj.d, obj.irts{obj.n_layers}, sirt_opt, 'debug_x', samples(:,ind2), 'sample_x', samples(:,ind1));
        else
            obj.irts{obj.n_layers+1} =  SIRT(newf, obj.d, oneds{poly_counter}, sirt_opt, 'debug_x', samples(:,ind2), 'sample_x', samples(:,ind1));
        end
    else
        if obj.n_layers > 1
            obj.irts{obj.n_layers+1} =  SIRT(newf, obj.d, obj.irts{obj.n_layers}, sirt_opt, 'sample_x', samples(:,ind));
        else
            obj.irts{obj.n_layers+1} =  SIRT(newf, obj.d, oneds{poly_counter}, sirt_opt, 'sample_x', samples(:,ind));
        end
    end
    obj.logz = obj.logz + log(obj.irts{obj.n_layers+1}.z);
    obj.n_evals = obj.n_evals + obj.irts{obj.n_layers+1}.n_evals;
    obj.n_layers = obj.n_layers + 1;
    if poly_counter < length(oneds)
        poly_counter = poly_counter + 1;
    else
        poly_counter = length(oneds);
    end
    
    % stop
    if abs(beta - 1) < 1E-10
        break;
    end
end

% evaluate the map and the target function
[x,logfx] = eval_irt(obj, samples);
[mllkds, mlps] = func(x);

% squared Hellinger error between logfx and (-mllkds-mlps)
[~,dh2,~] = f_divergence(logfx, -mllkds-mlps);

%
fprintf('\n\niteration=%2d, Hell error=%3.3e, cum#fevals=%3.3e\n', ...
    obj.n_layers, sqrt(dh2), obj.n_evals);
%
end