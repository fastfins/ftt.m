classdef SIRT < FTT
    
    properties
        irt_dir
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
            tmp = eval(oned, x(:), reshape(permute(core, [2,1,3]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_213_deri(oned, core, x)
            rkm = size(core, 1);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval_deri(oned, x(:), reshape(permute(core, [2,1,3]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rkm, []), [2,1,3]), rkm*nx, []);
        end
        
        function T = eval_oned_core_231(oned, core, x)
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = eval(oned, x(:), reshape(permute(core, [2,3,1]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        function T = eval_oned_core_231_deri(oned, core, x)
            rk  = size(core, 3);
            nn  = oned.num_nodes;
            nx  = length(x);
            % evaluate the updated basis function
            tmp = evalderi(oned, x(:), reshape(permute(core, [2,3,1]), nn, []));
            T   = reshape(permute(reshape(tmp, nx, rk, []), [2,1,3]), rk*nx, []);
        end
        
        [z,ys,ms] = cumint(ftt, dir)
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

        function obj = IRT(ftt, varargin)
            % Setup data structure used for IRT
            %
            %   ftt - A given function tensor train
            %
            %   err_tol - error tolerance for evaluating the IRT
            
            obj = ftt;
            %
            [obj.z, obj.ys, obj.ms] = SIRT.cumint(obj, obj.direction);
            obj.oned_cdfs = cell(size(obj.oneds));
            for i = 1:length(obj.oneds)
                obj.oned_cdfs{i} = CDFconstructor(obj.oneds{i}, varargin{:});
            end
        end
    end
    
end