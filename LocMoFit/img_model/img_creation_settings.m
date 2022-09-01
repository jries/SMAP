function out = img_creation_settings(modelName)
    if nargin == 0
        out = {'arc2D'};
        return
    else
        out.SMAP = [];          % general parameters
        out.LocMoFit = [];      % model sepecific parameters
    end
    out.SMAP = [];
    out.LocMoFit.parID = {};
    out.LocMoFit.value = [];
    out.LocMoFit.internalSettings = {};
    out.LocMoFit.internalSettings_val = [];
    switch modelName
        case 'arc2D'
            out.LocMoFit.parID = {'m1.mPar.radius', 'm1.mPar.theta'};
            out.LocMoFit.value = [50 270];
    end
end