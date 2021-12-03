classdef SymmetricReference < Reference
    
    properties
        mu
        sigma
        domain
        left
        right
    end
    
    methods (Abstract)
        [u,f] = eval_ref_cdf(obj, z)
        [f,g] = eval_ref_pdf(obj, z)
        z = invert_ref_cdf(obj, u)
        %
        [f,g] = log_joint_ref_pdf(obj, z)
    end
    
    methods
        
        function x = invert_cdf(obj, u)
            u(isinf(u)) = 1-eps;
            u(isnan(u)) = eps;
            u(u>(1-eps)) = 1-eps;
            u(u<eps) = eps;
            z = invert_ref_cdf(obj, u);
            x = z*obj.sigma + obj.mu;
        end
        
        function [u,f] = eval_cdf(obj, x)
            z = (x-obj.mu)/obj.sigma;
            [u,f] = eval_ref_cdf(obj, z);
            u(u > 1-eps) = 1-eps;
            u(u < eps) = eps;
            f = f/obj.sigma;
        end
        
        function[f,g] = eval_pdf(obj, x)
            z = (x-obj.mu)/obj.sigma;
            [f,g] = eval_ref_pdf(obj, z);
            f = f/obj.sigma;
            g = g/obj.sigma^2;
        end
        
        function[f,g] = log_joint_pdf(obj, x)
            z = (x-obj.mu)/obj.sigma;
            [f,g] = log_joint_ref_pdf(obj, z);
            f = f - log(obj.sigma)*size(x,1);
            g = g/obj.sigma;
        end
        
        function obj = set_domain(obj, domain)
            obj.domain = [min(domain), max(domain)];
            obj.left  = eval_cdf(obj, obj.domain(1));
            obj.right = eval_cdf(obj, obj.domain(2));
        end
        
        
        function z = random(obj, d, n)
            % pseudo random samples
            u = rand(d, n);
            u = u*(obj.right-obj.left) + obj.left;
            z = invert_cdf(obj, u);
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
        end
        
        function z = sobol(obj, d, n)
            % QMC samples using Sobol sequence
            S = sobolset(d);
            u = net(S,n);
            u = u'*(obj.right-obj.left) + obj.left;
            z = invert_cdf(obj, u);
            z(z >= obj.domain(2)) = obj.domain(2);
            z(z <= obj.domain(1)) = obj.domain(1);
        end
        
        function debug(obj, n)
            xs = linspace(-5, 5, n);
            [c,f] = eval_cdf(obj, xs);
            fd = diff(c)/(10/n);
            figure;
            subplot(3,2,1)
            plot(xs(1:n-1),fd);hold on;plot(xs,f)
            legend('FD', 'pdf')
            title('pdf')
            subplot(3,2,2)
            plot(xs,c)
            title('cdf')
            %
            disp(norm(invert_cdf(obj, c)-xs));
            %
            disp(sum(f)*(10/n))
            %
            [f,g] = eval_pdf(obj, xs);
            gd = diff(f)/(10/n);
            subplot(3,2,3)
            plot(xs(1:n-1),gd);hold on;plot(xs,g)
            legend('GD', 'grad')
            title('grad')
            subplot(3,2,4)
            plot(xs,f)
            title('pdf')
            
            [lf,lg] = log_joint_pdf(obj, xs);
            lgd = diff(lf)/(10/n);
            subplot(3,2,5)
            plot(xs(1:n-1),lgd);hold on;plot(xs,lg)
            legend('GD', 'grad')
            title('log grad')
            subplot(3,2,6)
            plot(xs,lf)
            title('log pdf')
        end
        
        function obj = SymmetricReference(varargin)
            defaultMu = 0;
            defaultSigma = 1;
            defaultDomain = [-4,4];
            %
            p = inputParser;
            addOptional(p,'mu',defaultMu,@(x) isnumeric(x) && isscalar(x));
            addOptional(p,'sigma',defaultSigma,@(x) isnumeric(x) && isscalar(x) && (x>=0));
            addOptional(p,'domain',defaultDomain,@(x) isnumeric(x));
            p.KeepUnmatched = false;
            parse(p,varargin{:});
            %
            obj.mu = p.Results.mu;
            obj.sigma = p.Results.sigma;
            %
            obj = set_domain(obj, p.Results.domain);
        end
    end
    
end