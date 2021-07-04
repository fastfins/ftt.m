classdef IRTmarginal < IRT
    methods
        function obj = IRTmarginal(ftt, ind, varargin)
            ftt_new = int_block(ftt, ind);
            obj@IRT(ftt_new, varargin{:});
        end
    end
end
