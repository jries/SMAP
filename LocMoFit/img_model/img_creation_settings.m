function out = img_creation_settings(modelName)
    if nargin == 0
        out = {'sphericalCap3D_surfaceArea','arc2D'};
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
        case 'sphericalCap3D_surfaceArea'
            out.LocMoFit.parID = {'m1.mPar.surfaceArea', 'm1.mPar.closeAngle', 'm91.sim.numOfMol'};
            out.LocMoFit.value = [5 120 1000];
            out.SMAP.general.view = {'xy', 'xz'};
            out.SMAP.settings.se_sitefov = 250;
            out.SMAP.settings.se_siteroi = 250;
    end
end