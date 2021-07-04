function data = pdf2cdf(oned_cdf, pdf)
%
%Define data structure for solving the inverse transform for given pdf(s)
%
%%%%%
%Inputs:
%
%oned_cdf:
%  The input oned polynomial for defining cdf
%
%pdf:
%  Pdf values, d x n, where n is the number of pdfs need to be sampled
%
%%%%%
%Output:
%
%data:
%  A data structure contains:
%
%  coef:        coefficients for defining the CDFs
%  cdf_nodes:   cdf function evaluated at nodes, used for root finding
%  cdf_grid:    only for Lagrange basis, cdf values at the end points of
%               each element
%  base:        base value of the intergrals for working out cdf 
%  norm:        normalisation of the integration
%  
%Tiangang Cui, August, 2019

if (sum(pdf(:)<0)>0)
    disp(['negative pdf ' num2str(sum(pdf(:)<0))])
end

data.size = size(pdf,2);

%if sum(pdf(:)<=0) > 1
%    pdf = pdf - min(pdf, [], 1) + 1E-5*max(pdf, [], 1);
%    disp('Warning: negative pdf values found')
%end

switch oned_cdf.type
    case{'Lagrange'}
        if data.size > 1
            data.pdf_left  = pdf(1,:);
            data.pdf_right = pdf(end,:);
            % 1st coord: local, 2nd coord: elems, 3rd: pdfs
            local_pdf   = reshape(pdf(oned_cdf.global2local,:), oned_cdf.local.num_nodes, oned_cdf.num_elems, data.size);
            % permute: 1st coord: local, 2nd coord: pdfs, 3rd: elems
            local_pdf   = permute(local_pdf, [1,3,2]);
            
            data.coef   = reshape(oned_cdf.local.node2basis*reshape(local_pdf, oned_cdf.local.num_nodes, []), ...
                oned_cdf.local.num_nodes, data.size, oned_cdf.num_elems);
            
            data.cdf_grid   = zeros(oned_cdf.num_elems+1, data.size);
            data.cdf_nodes  = zeros(oned_cdf.num_nodes, data.size);
            data.base       = zeros(oned_cdf.num_elems, data.size);
            switch oned_cdf.bc
                case{'Dirichlet'}
                    data.cdf_grid(1,:) = 0.5*data.pdf_left*oned_cdf.gs^2;
                otherwise
                    %bug line below
                    %data.cdf_grid(2,:) = data.pdf_left*oned_cdf.gs;
                    data.cdf_grid(1,:) = data.pdf_left*oned_cdf.gs;
            end
            ind = reshape(oned_cdf.global2local, [], oned_cdf.num_elems);
            for i = 1:oned_cdf.num_elems
                tmp             = (oned_cdf.cdf_basis2node*oned_cdf.jac(i))*data.coef(:,:,i);
                data.base(i,:)  = tmp(1,:);
                data.cdf_nodes(ind(:,i),:)  = tmp - data.base(i,:) + data.cdf_grid(i,:);
                data.cdf_grid(i+1,:)        = data.cdf_grid(i,:) + tmp(end,:) - data.base(i,:);
            end
            switch oned_cdf.bc
                case{'Dirichlet'}
                    data.norm = data.cdf_grid(end,:) + 0.5*data.pdf_right*oned_cdf.gs^2;
                otherwise
                    data.norm = data.cdf_grid(end,:) + data.pdf_right*oned_cdf.gs;
            end
            data.cdf_nodes = data.cdf_nodes./data.norm;
            data.cdf_grid  = data.cdf_grid ./data.norm;
        else
            data.pdf_left  = pdf(1);
            data.pdf_right = pdf(end);
            % 1st coord: local, 2nd coord: elems
            local_pdf   = reshape(pdf(oned_cdf.global2local), oned_cdf.local.num_nodes, oned_cdf.num_elems);
            data.coef   = oned_cdf.local.node2basis*local_pdf;
            %
            data.cdf_grid   = zeros(oned_cdf.num_elems+1, 1);
            data.cdf_nodes  = zeros(oned_cdf.num_nodes, 1);
            data.base       = zeros(oned_cdf.num_elems, 1);
            switch oned_cdf.bc
                case{'Dirichlet'}
                    data.cdf_grid(1) = 0.5*data.pdf_left*oned_cdf.gs^2;
                otherwise
                    data.cdf_grid(1) = data.pdf_left*oned_cdf.gs;
            end
            ind = reshape(oned_cdf.global2local, [], oned_cdf.num_elems);
            for i = 1:oned_cdf.num_elems
                tmp             = (oned_cdf.cdf_basis2node*oned_cdf.jac(i))*data.coef(:,i);
                data.base(i)    = tmp(1);
                data.cdf_nodes(ind(:,i)) = tmp - data.base(i) + data.cdf_grid(i);
                data.cdf_grid(i+1) = data.cdf_grid(i) + tmp(end) - data.base(i);
            end
            switch oned_cdf.bc
                case{'Dirichlet'}
                    data.norm = data.cdf_grid(end) + 0.5*data.pdf_right*oned_cdf.gs^2;
                otherwise
                    data.norm = data.cdf_grid(end) + data.pdf_right*oned_cdf.gs;
            end
            data.cdf_nodes = data.cdf_nodes/data.norm;
            data.cdf_grid  = data.cdf_grid /data.norm;
        end
        data.grid_pts = oned_cdf.grid_pts;
    otherwise
        data.coef       = oned_cdf.node2basis*pdf;
        data.cdf_nodes  = oned_cdf.cdf_basis2node*data.coef;
        data.base       = data.cdf_nodes(1,:);
        data.norm       = data.cdf_nodes(end,:) - data.base;
        data.cdf_nodes  = (data.cdf_nodes - data.base)./data.norm;
        
end

end