classdef FTT
    % FTT abstract class, need to use either FTTamen or FTTrand
    %
    % FTT Properties:
    %
    %   opt         - FTToption
    %
    %   direction   - ALS direction, >0: built from left to right
    %                                <0: built from right to left
    %
    %   oneds{k}    - data structure containing information for building
    %                 one dimensional polynomial representation
    %
    %   cores{k}    - nodal values or coefficent tensors of oned functions
    %                 the dimension of the current core is organised as
    %                 previous rank x oned{k}.num_nodes x current rank
    %
    %   interp_x{x} - interpolation coordinates
    %
    %
    % FTToption Methods:
    %
    %   eval        - evaluate fTT. The output is horizontally aligned.
    %
    %   eval_block  - evaluate fTT for either the first or last k variables.
    %                 The output is horizontally aligned.
    %
    %   round       - round the TT cores
    %
    %   int         - integrate the entire TT
    %
    %   int_block   - integrate a block of TT cores
    %
    %   size        - size of the TT
    %
    %%%%%%%%%%%%%%%%%
    %
    % Example (vector function outputs, m = 2):
    %
    % poly = setup_oned(5, 'type', 'Lagrange', 'lag_elems', 10, 'ghost_size', 1E-10);
    % func = @(x) [sqrt(1./sum(1E-5+x.^2,1)); sqrt(1./sum(1E-2+x.^2,1))];
    % d    = 10;
    %
    % debug_size = 1E4;
    % debug_x = zeros(d, debug_size);
    % for k = 1:d
    %     debug_x(k,:) = sample_oned_domain(poly, debug_size);
    % end
    %
    % % use alternating energy enrichment (AMEN)
    % opt1 = ftt_options('method', 'AMEN', 'oned_ref', poly, ...
    %     'err_tol', 1E-8, 'loc_err_tol', 1E-10, 'max_rank', 50);
    % ftt1 = build_ftt(func, d, [], opt1, 'debug_x', debug_x);
    %
    % % use random enrichment
    % opt2 = ftt_options('method', 'Random', 'oned_ref', poly, ...
    %     'err_tol', 1E-8, 'loc_err_tol', 1E-10, 'max_rank', 50);
    % ftt2 = build_ftt(func, d, [], opt2, 'debug_x', debug_x);
    %
    % % evaluate the function and its factorisations
    % exact   = func(debug_x);
    % approx1 = eval_ftt(ftt1, debug_x);
    % approx2 = eval_ftt(ftt2, debug_x);
    % % plot the error
    % figure
    % plot(exact(:) - approx1(:), 'x')
    % hold on
    % plot(exact(:) - approx2(:), '.')
    %
    %
    
    properties
        opt FTToption
        cores
        oneds
        interp_x
        direction
        res_x
        res_w
    end
        
    methods (Static)
        [core, interp_x, res_w, res_x, core_next] = build_basis_amen(oned, ...
            interp_xold, res_xold, res_w_l, res_w_r, core_next, F, Fu, Fr, ...
            dir, int_method, loc_err_tol, max_rank, kick_rank)
        % build tt core by amen

        [core,interp_x,core_next] = build_basis_svd(oned, ...
            interp_xold, core_next, F, ...
            dir, int_method, loc_err_tol, max_rank)
        % build tt core by svd

        [B,A,r] = local_truncate(loc_err_tol, min_rank, max_rank, oned, F)
        % truncate local core

        [f, f_evals] = local_block(sqrt_flag, oned, xleft, xright, func)
        % evaluate function for a coordinate

        interp_x = local_index(oned, direction, interp_xold, ind)
        % index selection

        rerr = local_error(core, f)
        % error of a core
    end
    
    methods
        
        obj = cross(obj, func, d, sample_x, debug_x)
        % cross iterations

        obj = round(obj, thres)
        % round the TT cores
        
        z = int(obj)
        % Integrate the entire TT
        
        ftt = int_block(obj, ind)
        % Integrate a block of TT cores
        
        fx = eval(obj, x, varargin)
        % Evaluate the ftt function. The output is horizontally aligned.
        
        fx = eval_block(obj, x, dir)
        % Evaluate the fTT for either the first or last k variables.
        
        [d,rs,ns] = size(obj)
        % size of ftt

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj = FTT(func, d, arg, varargin)
            % Construct tensor train for a function mapping from R^d to R^m.
            %
            % Inputs:
            %
            %   func - a function (R^d to R^m) that take inputs as a dxn
            %          matrix and returns mxn vector
            %
            %   d    - dimension of the input variable
            %
            %   arg  - either an existing ftt used as the initial guess, or
            %          a set of one dimensional bases for discretising the
            %          function. If only one set of basis is supplied, each
            %          coordinate is discretised using the same basis
            %
            %   opt  - FTT options
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % parsing inputs
            defaultOption    = FTToption();
            defaultSampleSet = [];
            defaultDebugSet  = [];
            
            p = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            %
            addRequired(p,'func',@(x) isa(x, 'function_handle'));
            addRequired(p,'d',validScalarPosNum);
            addRequired(p,'arg');
            addOptional(p,'option',defaultOption);
            %
            addParameter(p,'sample_x', defaultSampleSet);
            addParameter(p,'debug_x',  defaultDebugSet);
            %
            p.KeepUnmatched = false;
            parse(p,func,d,arg,varargin{:});
            %
            obj.opt = p.Results.option;
            sample_x    = p.Results.sample_x;
            debug_x     = p.Results.debug_x;
            %
            if isa(arg, 'FTT')
                obj = arg;
            else
                obj.oneds = cell(d,1);
                if isa(arg, 'cell')
                    for k = 1:d
                        if ~isa(arg{k}, 'oned')
                            error('wrong type of argument')
                        end
                        obj.oneds{k} = arg{k};
                    end
                    obj.cores = [];
                    obj.res_x = [];
                    obj.res_w = [];
                elseif isa(arg, 'oned')
                    for k = 1:d
                        obj.oneds{k} = arg;
                    end
                    obj.cores = [];
                    obj.res_x = [];
                    obj.res_w = [];
                else
                    error('wrong type of argument')
                end
            end
            % build ftt
            obj = cross(obj, func, d, sample_x, debug_x);
        end
        
    end
end