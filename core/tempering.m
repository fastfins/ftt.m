classdef tempering
    
    properties
        adapt_flag
        betas
        n_layers
        max_layers
        qmc_flag
        n_samples
        ess_tol
        beta_factor
        max_beta_factor
    end
    
    properties (Access = private, Constant = true)
        defaultSQRTFlag = false;
        defaultMaxALS   = 4;
        defaultALSTol   = 1E-2;
        defaultInitRank = 10;
        defaultKickRank = 2;
        defaultMaxRank  = 20;
        defaultLocTol   = 1E-10;
        defaultCDFTol   = 1E-10;
        defaultTTMethod = 'amen';
        expectedTTMethod  = {'random','amen'};
        defaultIntM     = 'MaxVol';
        expectedIntM    = {'QDEIM','WDEIM','MaxVol'};
    end
    
    methods (Static)
        function ess2n = ess_ratio(beta_p, beta, mllkds)
            log_weight = (beta_p - beta)*mllkds;
            log_weight = log_weight - max(log_weight);
            ess2n = ( sum(exp(log_weight)).^2/sum(exp(2*log_weight)) ) / length(mllkds(:));
        end
    end
    
    methods
        function obj = tempering(varargin)
            obj.n_layers = 1;
            if obj.adapt_flag
                obj.betas(obj.n_layers) = beta_p;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function beta = current_temperature(obj)
            beta = obj.betas(obj.n_layers);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [beta,ess] = new_temperature(obj, beta_p, mllkds)
            if obj.adapt_flag 
                if nargin < 3
                    error('needs minus log-likelihood')
                end
                beta = beta_p;
                % compute ess over sample size
                ess = 1;
                while ess > obj.ess_tol && beta <= obj.max_beta_factor*beta_p
                    beta = beta*obj.beta_factor;
                    ess = adaptTempering.ess_ratio(beta_p, beta, mllkds);
                end
                beta = min(1, beta);
                ess = adaptTempering.ess_ratio(beta_p, beta, mllkds);
                obj.n_layers = obj.n_layers+1;
                obj.betas(obj.n_layers) = beta;
                if obj.n_layers > obj.max_layers
                    error('exceeds the maximum number of layers')
                end
            else
                obj.n_layers = obj.n_layers+1;
                beta = obj.betas(obj.n_layers);
                ess = 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [beta,ess] = new_temperature_double(obj, beta_p, mllkds, mlps)
            if obj.adapt_flag
                if nargin < 4
                    error('needs minus log-likelihood and minus log-prior')
                end
                beta = beta_p;
                % compute ess over sample size
                ess = 1;
                while ess > obj.ess_tol && beta <= obj.max_beta_factor*beta_p
                    beta = beta*obj.beta_factor;
                    ess = adaptTempering.ess_ratio(beta_p, beta, mllkds);
                end
                beta = min(1, beta);
                ess = adaptTempering.ess_ratio(beta_p, beta, mllkds);
                obj.n_layers = obj.n_layers+1;
                obj.betas(obj.n_layers) = beta;
                if obj.n_layers > obj.max_layers
                    error('exceeds the maximum number of layers')
                end
            else
                obj.n_layers = obj.n_layers+1;
                beta = obj.betas(obj.n_layers);
                ess = 1;
            end
        end
    end
end