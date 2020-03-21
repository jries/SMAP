classdef uTrack<interfaces.DialogProcessor
    % calls the uTrack single particle tracking software from the Danuser
    % Lab
    %Please install from: https://github.com/DanuserLab/u-track
    %add directory here in the plugin.
%     Published in: 1.Jaqaman, K., Loerke, D., Mettlen, M., Kuwata, H., Grinstein, S.,
%     Schmid, S. L. & Danuser, G. Robust single-particle tracking in
%     live-cell time-lapse sequences. Nat Methods 5, 695–702 (2008).

    methods
        function obj=uTrack(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson'};
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=uTracki(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end
function out=uTracki(obj,p)
if ~isdeployed
addpath(genpath('External/u-track/software'))
addpath(genpath('/Users/jonasries/Documents/MATLAB/u-track/software'))
end

%prepare input parameters

locs=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','phot','frame'},'layer',1,'position','roi');
[sortf,indsort]=sort(locs.frame);
ind1=1;
for ff=1:sortf(end)
    while ind1<=length(sortf)&&sortf(ind1)<ff
        ind1=ind1+1;
    end
    ind2=ind1;
    while ind2<=length(sortf)&&sortf(ind2)==ff
        ind2=ind2+1;
    end   
    ind2=ind2-1;
    rangeh=ind1:ind2;
    ind1=ind2+1;

    movieInfo(ff).xCoord=horzcat(locs.xnm(indsort(rangeh)), locs.locprecnm(indsort(rangeh)))/100;
    movieInfo(ff).yCoord=horzcat(locs.ynm(indsort(rangeh)), locs.locprecnm(indsort(rangeh)))/100;
    if ~isempty(locs.znm)
    movieInfo(ff).zCoord=horzcat(locs.znm(indsort(rangeh)), locs.locprecznm(indsort(rangeh)))/100;
    end
    movieInfo(ff).amp=horzcat(locs.phot(indsort(rangeh)), zeros(size(locs.phot(indsort(rangeh)))));
end

global tracksFinal
[gapCloseParam,costMatrices,kalmanFunctions]=getDefaultParameters;
%% additional input

%saveResults
% saveResults.dir = 'C:\kjData\test\140825_Sungsoo\'; %directory where to save input and output
% saveResults.filename = 'tracksTest4.mat'; %name of file where input and output are saved
saveResults = 0; %don't save results

%verbose state
verbose = 1;

%problem dimension
probDim = 2;
[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose)


plotTracks2D(tracksFinal)
out=[];

end
         

function [gapCloseParam,costMatrices,kalmanFunctions]=getDefaultParameters
%% general gap closing parameters
%
% Copyright (C) 2016, Danuser Lab - UTSouthwestern 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
gapCloseParam.timeWindow = 2; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 1; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

%parameters

parameters.linearMotion = 0; %use linear motion Kalman filter.
parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = 2; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
% parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.

%optional input
parameters.diagnostics = []; %if you want to plot the histogram of linking distances up to certain frames, indicate their numbers; 0 or empty otherwise. Does not work for the first or last frame of a movie.

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

%parameters

%needed all the time
parameters.linearMotion = 0; %use linear motion Kalman filter.

parameters.minSearchRadius = 2; %minimum allowed search radius.
parameters.maxSearchRadius = 2; %maximum allowed search radius.
parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.

parameters.brownScaling = [0 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
% parameters.timeReachConfB = 3; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
parameters.timeReachConfB = gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).

parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.

parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.

parameters.useLocalDensity = 0; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.linStdMult = 1*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.

parameters.linScaling = [0.25 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
% parameters.timeReachConfL = 4; %similar to timeReachConfB, but for the linear part of the motion.
parameters.timeReachConfL = gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';


end
             
function pard=guidef
pard.t1.object=struct('String','u-Track from Danuser lab','Style','text');
pard.t1.position=[1,2];
pard.t1.Width=4;

pard.t2.object=struct('String','time window (f)','Style','text');
pard.t2.position=[2,1];
pard.t2.Width=.75;
pard.t2.TooltipString='Gap closing time window (10 frames). Depends on SNR and fluorophore blinking. Critical if too small or too large. Robust if in proper range.';

pard.timeWindow.object=struct('String','10','Style','edit');
pard.timeWindow.position=[2,1.75];
pard.timeWindow.Width=.25;
pard.timeWindow.TooltipString=pard.t2.TooltipString;

pard.mergeSplit.object=struct('String','Merge and Split','Style','checkbox');
pard.mergeSplit.position=[2,2];
pard.mergeSplit.Width=1;
pard.mergeSplit.TooltipString='Flag for merging and splitting (1). 1 to allow merging and splitting between track segments, 0 otherwise';

pard.minTrackLent.object=struct('String','time window (f)','Style','text');
pard.minTrackLent.position=[2,3];
pard.minTrackLent.Width=.75;
pard.minTrackLent.TooltipString='Minimum track segment length used in the gap closing, merging and splitting step (2 frames). Since false positives during detection generally lead to tracks with short lifetimes, this parameter excludes very short tracks from participating in the gap closing, merging and splitting step. May need adjustment for very fast movies.';

pard.minTrackLen.object=struct('String','2','Style','edit');
pard.minTrackLen.position=[2,3.75];
pard.minTrackLen.Width=.25;
pard.minTrackLen.TooltipString=pard.minTrackLent.TooltipString;

pard.linearMotion.object=struct('String','Linear motion','Style','checkbox');
pard.linearMotion.position=[3,1];
pard.linearMotion.Width=1;
pard.linearMotion.TooltipString='Flag for linear motion (0). 1 for position propagation prior to linking using a linear motion model (as described in the text), 0 otherwise (hence assuming random motion).';


pard.plugininfo.description=sprintf('uTrack');
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.Name='uTrack';
end

