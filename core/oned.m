classdef oned
    
    properties
        domain(1,2)
        order
        %
        num_nodes
        nodes(:,1)
    end
    
    methods (Abstract)
        eval(obj)
        eval_deri(obj)
        eval_basis(obj)
        eval_basis_deri(obj)
        %
        sample_domain(obj)
        mass_r(obj)
        integral(obj)
        get_name(obj)
        %
        point_selection(obj)
    end
    
end