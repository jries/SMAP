function output = getSolverOption(solverName, varargin)
%% GETSOLVEROPTION Configuration of the solver options
switch solverName
    case 'fminsearchbnd'
        options = {'Display','FunValCheck','MaxFunEvals','MaxIter','TolFun','TolX','PlotFcn','OutputFcn'};
        if length(varargin)==0
            output = options;
        else
        end
    case 'fmincon'
        options = {'Display','PlotFcn'};
        if length(varargin)==0
            output = options;
        else
        end
    case 'particleswarm'
        options = {'Display','UseParallel','UseVectorized','FunctionTolerance','FunValCheck','HybridFcn','MaxIterations','MaxStallIterations','MaxStallTime','MaxTime','MinNeighborsFraction','ObjectiveLimit','SelfAdjustmentWeight','SocialAdjustmentWeight','SwarmSize','OutputFcn'};
        if length(varargin)==0
            output = options;
        else
        end
end
end