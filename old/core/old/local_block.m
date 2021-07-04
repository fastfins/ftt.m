function [f, f_evals] = local_block(ng_flag, oned, xleft, xright, func)
%
%Evaluate the user function at cross indices and qudrature/intepolation points
%weighed by the interplation matrix
%
%Tiangang Cui, August, 2019

if isempty(xleft)
    % left boudnary
    niright = size(xright, 2);
    % f       = zeros(niright, oned_def.nquad);
    % define parameters, each xright binding with a block of xquad
    param   = [ repmat(oned.nodes(:)', 1, niright); ...
        reshape(repmat(xright, oned.num_nodes, 1), size(xright,1), niright*oned.num_nodes) ];
    % return of func is a vector, reshape to a matrix
    % 1st index: xquad, 2nd: xright
    f       = reshape( func(param)', 1, oned.num_nodes, niright, []);
elseif isempty(xright)
    % right boundary
    nileft  = size(xleft, 2);
    % f       = zeros(nileft, oned_def.nquad);
    % define parameters, each xquad binding with a block of xleft
    param   = [ repmat(xleft, 1, oned.num_nodes); ...
        reshape( repmat(oned.nodes(:)', nileft, 1), 1, nileft*oned.num_nodes ) ];
    % return of func is a vector, reshape to a matrix
    % 1st index: xleft, 2nd: xquad
    f       = reshape( func(param)', nileft, oned.num_nodes, 1, []);
else
    %
    nileft  = size(xleft, 2);
    niright = size(xright, 2);
    % f       = zeros(nileft, niright, oned_def.nquad);
    
    tmp_lq  = [ repmat(xleft, 1, oned.num_nodes); ...
        reshape( repmat(oned.nodes(:)', nileft, 1), 1, nileft*oned.num_nodes ) ];
    
    param   = [ repmat(tmp_lq, 1, niright); ...
        reshape(repmat(xright, oned.num_nodes*nileft, 1), size(xright,1), nileft*oned.num_nodes*niright) ];
    
    % return of func is a vector, reshape to tensor
    f       = reshape( func(param)', nileft, oned.num_nodes, niright, []);
end

if ng_flag
    if any(f(:)<0)
        error('Error: negative function evaluations')
    end
    f   = f.^(0.5);
end

f_evals = size(param, 2);

end