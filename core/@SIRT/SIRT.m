classdef SIRT < FTT
    
    properties
        marginal_direction
        z
        ys
        ms
        oned_cdfs
    end
    
    methods (Static)
        function T = eval_oned_core_213(oned, core, x)
            rkm = size(core, 1);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval(oned, reshape(permute(core, [2,1,3]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_213_deri(oned, core, x)
            rkm = size(core, 1);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval_deri(oned, reshape(permute(core, [2,1,3]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_231(oned, core, x)
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval(oned, reshape(permute(core, [2,3,1]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        function T = eval_oned_core_231_deri(oned, core, x)
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = evalderi(oned, reshape(permute(core, [2,3,1]), nn, []), x(:));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        % [z,ys,ms] = cumint(ftt, dir)
        % Marginalise the pdf represented by ftt dimension by dimension
        
    end
    
    methods
        z = eval_rt(obj, r)
        % Evaluate squared RT z = T(r), where z is uniform and r is target r.v.

        [r,f] = eval_irt(obj, z)
        % Evaluate squared IRT r = T(z), where z is uniform

        [r,f] = eval_cond_irt(obj, x, z)
        % Using SIRT to draw conditional samples from the target pdf 
        % approximated by FTT

        J = eval_rt_jac(firt, r, z)
        % Evaluate the jacobian of the squared RT z = T(r), where z is 
        % uniform and r is target r.v.

        fx = eval_marginal_pdf(obj, x)
        % Evaluate the marginalise pdf represented by ftt

        obj = marginalise(obj, dir) 
        % Marginalise the pdf represented by ftt dimension by dimension
        
        function obj = SIRT(func, d, arg, varargin)
            % Setup data structure used for IRT
            
            defaultErrTol = 1E-8;
            %
            obj@FTT(func, d, arg, varargin{:})
            %
            if strcmp(obj.opt.tt_method, 'amen')
                obj = round(obj);
            end
            obj.oned_cdfs = cell(size(obj.oneds));
            for i = 1:length(obj.oneds)
                obj.oned_cdfs{i} = CDFconstructor(obj.oneds{i}, defaultErrTol);
            end
        end
    end
    
end