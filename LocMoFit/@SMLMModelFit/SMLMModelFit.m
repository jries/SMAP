classdef SMLMModelFit<LocMoFit
    % This class is obselete. Please use LocMoFit instead. SMLMModelFit is
    % now defined as a subclass of LocMoFit during the transition phase.
    methods
        function obj = SMLMModelFit(varargin)
            obj@LocMoFit(varargin{:});
        end
    end
end