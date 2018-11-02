classdef boundaryFinder<interfaces.SEEvaluationProcessor
    % This plug-in depends on the BALM_fibril_growth.
    % Green line is the original boundary
    % White line is the refined boundary

    properties
        boundary
    end
    methods
        function obj=boundaryFinder(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, inp)
            %% Data import
            % Import the kimograph
            kimo = obj.site.evaluation.BALM_fibril_growth.kimograph;
            kimo(kimo>=1)=1;

            slideStep = inp.slideStep;
            gridRange = inp.gridRange;
            adjM = inp.adjM;

            % extend the kimo on x to compensate the sliding window
            oriSize = size(kimo);
            expand = zeros(oriSize(1),gridRange*slideStep);
%             kimo = [kimo, expand];

            % Binarize the kimograph
            % X means time, Y means space
            [kimoRefx, kimoRefy] = find(kimo);
            [kimoNZx, kimoNZy] = find(kimo);
            
            %% Remove noise based on euclidean distances between data points
            % use euclidean distances to filter out noises
            K = round(length(kimoRefx)*5/100);

            ED = pdist2([kimoNZx kimoNZy],[kimoNZx kimoNZy],'euclidean');
            sortedED = sort(ED,2);
            sortedED = sortedED(:,1:K);
            rmean = mean(sortedED,2);

            EDcutoff = quantile(rmean, 0.95);

            kimo(sub2ind(size(kimo), kimoNZx([rmean>=EDcutoff]), kimoNZy([rmean>=EDcutoff]))) = 0;
            
            %% Get a rough boundary based on density
            % Get cordinates of the kimograph
                        
            Size = size(kimo);
            kimo = [kimo kimo(:,Size(2):-1:1)];                     % generate the mirrored kimograph
%             kimo = [kimo(Size(2):-1:1,:); kimo];
%             kimo(1:Size(1),Size(2)+1:end) = 0;
%             figure; imagesc(kimo)
            
            bwKimo = kimo>0;
            h = fspecial('gaussian', inp.std*2,inp.std);
            KDM = filter2(h, bwKimo);
%             allVal = KDM(:);
%             allVal(allVal==0)=[];
%             maxVal = prctile(allVal,80);
%             KDM(KDM>maxVal)=maxVal;
            KDM = KDM(:,end:-1:1)';
%             KDM = getKernelMatrix(Size*2,[kimoFZx, kimoFZy]);
            
            
%             figure; imagesc(KDM)                                  % show the gaussian-filtered kimograph
            tempPlot = figure('Name','Temp');
            ax1 = axes(tempPlot);
            
            [~,ct] = contour(ax1, KDM, 100);
            
            tl = ct.LevelList(inp.contourLevel);
            close(tempPlot)
%             tl2 = ct.LevelList(27);
%             tl3 = ct.LevelList(18);
%             figure(999)
%             cla
%             imagesc(kimo)
%             
%             figure(1000)
%             cla
%             contour(KDM, 200);
%             hold on
            [theLevel,~] = contour(KDM, [tl tl],'LineWidth',2);
%             [theLevel2,~] = contour(KDM, [tl2 tl2],'LineWidth',2);
%             [theLevel3,~] = contour(KDM, [tl3 tl3],'LineWidth',2);
            rBoundary = round(theLevel)';
            rBoundary(rBoundary(:,2)>Size(2),:)=[];
%             rBoundary(rBoundary(:,1)<=Size(1),:)=[];
%             rBoundary(:,1)=rBoundary(:,1)-Size(1)+1;
%             rBoundary2 = round(theLevel2)';
%             rBoundary3 = round(theLevel3)';
%             ulimX = prctile(rBoundary(:,2), 99);
%             rmRBoundary2 = rBoundary2(:,2) <= ulimX;
%             rBoundary2(rmRBoundary2,:)=[];
%             rBoundary = [rBoundary;rBoundary2];
%             ulimX = prctile(rBoundary(:,2), 99);
%             rmRBoundary3 = rBoundary3(:,2) <= ulimX;
%             rBoundary3(rmRBoundary3,:)=[];
%             rBoundary = [rBoundary;rBoundary3];
            %scatter(kimoRefx, kimoRefy)
            %hold off
            boundary = zeros(Size(1),1);
            CxFParent = 1:Size(1);
            for i=CxFParent                                         % for each x
                idx = rBoundary(:,1)==i;
                if sum(idx)>0
                    boundary(i) = max(rBoundary(idx,2));
                end
            end
            CyFParent = cummax(boundary);
            CxFParent = CxFParent';
%             figure(893); plot(CxFParent, CyFParent) 
%             figure(894); plot(rBoundary(:,1), rBoundary(:,2),' o') 
            %% refinement of the boundary
            Mref = calMeasurement(CxFParent,CyFParent,kimoNZx,kimoNZy,Size);
            MrefO = Mref;
            CyFParentA = CyFParent;
            for i=size(CxFParent):-1:1
               thisFrameY = kimoNZy(kimoNZx == CxFParent(i));
               thisFrameY = thisFrameY(thisFrameY>CyFParent(i));
               ii = 1;
               futherStep = 0;
               while ii <= length(thisFrameY)
                   CyFParentO = CyFParentA;

                   % New boundary
                   CyFParentA(i)=thisFrameY(ii);

                   CyFParentA = cummax(CyFParentA);

                   % fit to new boundary
                   Mnew = calMeasurement(CxFParent,CyFParentA,kimoNZx,kimoNZy,Size);
                   if Mnew >= MrefO*adjM && Mnew >= Mref %#Par
                       MrefO = Mnew;
                       futherStep = 0;
                   else
                       if futherStep <= 5
                        CyFParentA = CyFParentO;
                        Mnew = MrefO;
                        futherStep = futherStep+1;
                       else
                        CyFParentA = CyFParentO;
                        Mnew = MrefO;
                        break
                       end
                   end
                   ii=ii+1;
               end
                if isinteger(i/10)
                    figure(51)
                    scatter(kimoNZx, kimoNZy)
                    hold on
                    plot(CxFParent, CyFParentA)
                    hold off
                end
            end
            
            KimoPar = obj.site.evaluation.BALM_fibril_growth;
            fig = KimoPar.kimograph;
            co=quantile(fig(:),0.999);
            fig(fig>co)=co;
            out.boundary = [CxFParent CyFParentA];

            
            %% Merge small steps
            steps = inp.mergeSteps;
            for i=1:length(steps)
                switch steps(i)
                    case 1
                    % Merge steps by move
                        idxRef = [];
                        while 1
                            [newboundary, idx] =  mergeStallsByMove(out.boundary, inp.minWidth);
                            if isequal(idxRef, idx)
                                break
                            end
                            idxRef = idx;
                            out.boundary =  newboundary;
                            out.boundary = rmRedundant(out.boundary);
                        end
                    case 2
                    % Merge steps by move
                    	idxRef = [];
                        while 1
                            [newboundary, idx] =  mergeByMove(out.boundary, inp.minWidth);
                            if isequal(idxRef, idx)
                                break
                            end
                            idxRef = idx;
                            out.boundary =  newboundary;
                            out.boundary = rmRedundant(out.boundary);
                        end
                    case 3
                    % Merge steps by time
                    out.boundary = mergeByTime(out.boundary, inp.minT);
                    out.boundary = rmRedundant(out.boundary);
                    case 0
                        % Remove redundant points (in stalled region)
                    out.boundary = rmRedundant(out.boundary);
                end
            end
            
            % Remove redundant points (in stalled region)
            stepWidth = diff(out.boundary(:,2));
            indStepWidthNon0 = find(stepWidth);
            indPoint2Keep = [indStepWidthNon0; indStepWidthNon0+1];
            indPoint2Keep = unique(indPoint2Keep);
            out.boundary = out.boundary([1; indPoint2Keep; end],:);
            
            out.stallTime = calStallTime(out.boundary);
            out.stepWidth = calStepWidth(out.boundary);
            out.growthTime = calGrowthTime(out.boundary);
            out.avgRate = calAvgRate(out.boundary);
            
            % Visualize the boundary
            h=obj.setoutput('kimograph');
            imagesc(h,(fig))
            hold(h,'on')
            plot(h,CyFParentA,CxFParent, 'LineWidth', 1, 'Color', 'g')
            plot(h,out.boundary(:,2),out.boundary(:,1), 'LineWidth', 1.5, 'Color', 'w')
            hold(h,'off')            
            

            %xlabel(h,'xnm')
            %ylabel(h,'frame')
            %xticklabels(h, KimoPar.xn)
            %yticklabels(h, KimoPar.fr)


            
            if 0
                out.avgRate = (CyFParentA(end)-CyFParentA(1))/(CxFParent(end)-CxFParent(1));
                out.stepWidth = stepWidthNon0;
                out.stallTime = stallTime;

                h2=obj.setoutput('statistics');
                axes(h2);
                ax1 = subplot(1,2,1);
                histogram(ax1, stepWidthNon0, 'BinWidth', 1);
                title(ax1,'Rate');
                ax2 = subplot(1,2,2);
                histogram(ax2, stallTime, 'BinWidth', 1);
                title(ax2,'Stall time');
            end
        end
     
        function pard=guidef(obj)
            pard=guidef;
        end
    end

end



function pard=guidef
pard.t_gridRange.object=struct('Style','text','String','Bin size of the grids');
pard.t_gridRange.position=[1,1];
pard.t_gridRange.Width=2;

pard.gridRange.object=struct('Style','edit','String',5);
pard.gridRange.position=[1,3];
pard.gridRange.TooltipString = 'If you set it as 5, it means before the density comparison, every grid will be set to cover a 5-by-5 area in the original coordinates';
            
pard.t_slideStep.object=struct('Style','text','String','Slide step(s)');
pard.t_slideStep.position=[2,1];
pard.t_slideStep.Width=2;

pard.slideStep.object=struct('Style','edit','String',5);
pard.slideStep.position=[2,3];
pard.slideStep.TooltipString = 'If you set it as 5, it means during the density comparison, every grid value will be campared to its following 4 (5 minus 1, which means the reference grid itself) right neighbors';

pard.t_adjM.object=struct('Style','text','String','Adjustment of M');
pard.t_adjM.position=[3,1];
pard.t_adjM.Width=2;

pard.adjM.object=struct('Style','edit','String',1.0003);
pard.adjM.position=[3,3];
pard.adjM.TooltipString = 'If you set it as 1.0003, it means during the optimization, if the measurment of current step (Mcur) is 0.0003-time worse than the measurment of the previous step (Mpre), Mcur will still be considered as a good result. The measurment, which defines the boundary is good or not, can be definde by users.';

pard.t_mergeSteps.object = struct('Style','text','String','Order of merging steps');
pard.t_mergeSteps.position=[4,1];
pard.t_mergeSteps.Width = 2;

pard.mergeSteps.object = struct('Style','edit','String','0 3 1 2 3 1 2 0');
pard.mergeSteps.position=[4,3];
pard.mergeSteps.TooltipString = '0 = clean up; 1 = by space(ambiguous); 2 = by space(small), 3 = by time';

pard.t_minT.object = struct('Style','text','String','Minimum time');
pard.t_minT.position=[5,1];
pard.t_minT.Width = 2;

pard.minT.object = struct('Style','edit','String',5);
pard.minT.position=[5,3];
pard.minT.TooltipString = 'minimum time of a step (arbitrary unit)';

pard.t_minWidth.object = struct('Style','text','String','Minimum width');
pard.t_minWidth.position=[6,1];
pard.t_minWidth.Width = 2;

pard.minWidth.object = struct('Style','edit','String',2);
pard.minWidth.position=[6,3];
pard.minWidth.TooltipString = 'minimum width of a step (arbitrary unit)';

pard.t_std.object = struct('Style','text','String','Std.');
pard.t_std.position=[7,1];
pard.t_std.Width = 0.5;

pard.std.object = struct('Style','edit','String',100);
pard.std.position=[7,1.5];
pard.std.TooltipString = 'minimum width of a step (arbitrary unit)';
pard.std.Width = 0.5;

pard.t_contourLevel.object = struct('Style','text','String','Level');
pard.t_contourLevel.position=[7,2];
pard.t_contourLevel.Width = 0.5;

pard.contourLevel.object = struct('Style','edit','String',90);
pard.contourLevel.position=[7,2.5];
pard.contourLevel.TooltipString = 'The level (out of 100) chosen as the outline of the rough boundary';
pard.contourLevel.Width = 0.5;

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';

end

function M = calMeasurement(x,y,qx,qy,Size)
    ref = interp1(x,y,qx);
    rightIdx = qy > ref;
    leftIdx = ~rightIdx;
    leftA = sum(y);
    rightA = Size(1)*Size(2)-sum(y);
    leftD = sum(leftIdx)/leftA;
    rightD = sum(rightIdx)/rightA;
    M = 0-((leftD-1)^2 + (rightD-0)^2)^(1/2);
end

function [newBoundary, idx] = mergeStallsByMove(boundary, minWidth)
    % Merge steps by move
    move = diff(boundary(:,2));
    stalled = move==0;
    stalled = [boundary(find(stalled), 1) boundary(find(stalled)+1, 1)];
    indConStalled = find(diff(stalled,1,2)>1);
    conStalled = stalled(indConStalled,:);
    conStallTime = diff(conStalled,1,2);

    smallMove = move <= minWidth;
    smallMove = [boundary(find(smallMove), 1) boundary(find(smallMove)+1, 1)];

    afterStall = ismember(smallMove(:,1),conStalled(:,2));
    beforeStall = ismember(smallMove(:,2), conStalled(:,1));
    abStallSM = smallMove(afterStall+beforeStall==2,:);
    lPosition2change = [];
    for i = 1:size(abStallSM,1)
    upstreamStall = conStalled(:,2) == abStallSM(i,1);
    downstreamStall = conStalled(:,1) == abStallSM(i,2);
        if conStallTime(upstreamStall) >= conStallTime(downstreamStall);
            lPosition2change = ismember(boundary(:,1), conStalled(downstreamStall,:)');
            lPositionValue = ismember(boundary(:,1), conStalled(upstreamStall,:)');
            boundary(lPosition2change,2) = boundary(lPositionValue,2);
        else
            lPosition2change = ismember(boundary(:,1), conStalled(upstreamStall,:)');
            lPositionValue = ismember(boundary(:,1), conStalled(downstreamStall,:)');
            boundary(lPosition2change,2) = boundary(lPositionValue,2);
        end
    end
    newBoundary = boundary;
    idx=lPosition2change;
end
function [newBoundary, idx] = mergeByMove(boundary, minWidth)
    move = diff(boundary(:,2));
    point2drag = move <= minWidth;
    indPoint2drag = find(point2drag)+1;
    boundary(indPoint2drag,2) = boundary(indPoint2drag-1,2);
    newBoundary = boundary;
    idx = indPoint2drag;
end
function newBoundary = mergeByTime(boundary, minT)
    point2keepCp = [];
    while 1
        stallTime = diff(boundary(:,1));
        point2keep = stallTime >= minT;
        if isequal(point2keep, point2keepCp)
            break
        end
        point2keepCp = point2keep;
        indPoint2keep=unique([find(boundary(:,1)==1); find(point2keep); find(point2keep)+1]);
        point2keep(indPoint2keep)=1;
        boundary = boundary(point2keep,:);
        newBoundary = boundary;
    end
end

function salltime = calStallTime(boundary)
    move = diff(boundary(:,2));
    stalled = move==0;
    stalled = [boundary(find(stalled), 1) boundary(find(stalled)+1, 1)];
    salltime = diff(stalled,1,2);    
end

function growthTime = calGrowthTime(boundary)
    move = diff(boundary(:,2));
    growth = move>0;
    growth = [boundary(find(growth), 1) boundary(find(growth)+1, 1)];
    growthTime = diff(growth,1,2);    
end

function stepWidth = calStepWidth(boundary)
    move = diff(boundary(:,2));
    stepWidth = move(move>0);
end

function avgRate = calAvgRate(boundary)
    startNEnd = boundary([1 end], :);
    diffTimeNPosition = diff(startNEnd, 1, 1);
    avgRate = diffTimeNPosition(2)/diffTimeNPosition(1);
end

function [z, bw] = getKernelMatrix(matrixSize, sVrot, varargin) 
    [meshx, meshy] = meshgrid(0:matrixSize(1), 0:matrixSize(2)); % make full grid
    q = [meshx(:), meshy(:)]; % convert the grid to positions
    if length(varargin)==0
        [sVrotx,xy,bw]=ksdensity(sVrot, q); % get kernel density estimation
    else
        [sVrotx,xy,bw]=ksdensity(sVrot, q, varargin{1,1}, varargin{1,2}); % get kernel density estimation    
    end
    
    % converted into an image
    idx = sub2ind(matrixSize(end:-1:1)+1, xy(:,2)+1, xy(:,1)+1); 
    z = zeros(matrixSize(end:-1:1)+1);
    z(idx) = sVrotx;
end

function newBoundary = rmRedundant(boundary)
    % Remove redundant points (in stalled region)
    stepWidth = diff(boundary(:,2));
    indStepWidthNon0 = find(stepWidth);
    indPoint2Keep = [indStepWidthNon0; indStepWidthNon0+1];
    indPoint2Keep = unique(indPoint2Keep);
    newBoundary = boundary([1; indPoint2Keep; end],:);
end 