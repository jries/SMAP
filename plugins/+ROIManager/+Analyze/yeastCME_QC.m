classdef yeastCME_QC<interfaces.DialogProcessor&interfaces.SEProcessor
%     Segements yeast cells with labeled endocytic sites. Foucs on top or bottom pole.
    methods
        function obj=yeastCME_QC(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
          cellLevelQC(obj,p)
          siteLevelQC(obj,p)
          out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard=[];
pard.t1.object=struct('String','cutoff','Style','text');
pard.t1.position=[1,1];
pard.cutoff.object=struct('String','1','Style','edit');
pard.cutoff.position=[1,2];

pard.preview.object=struct('String','preview','Style','checkbox','Value',1);
pard.preview.position=[3,1];
pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='Segements yeast cells with labeled endocytic sites. Foucs on top or bottom pole.';

% pard.t2.object=struct('String','sigmaNMS','Style','text');
% pard.t2.position=[2,1];
% pard.sigmaNMS.object=struct('String','5','Style','edit');
% pard.sigmaNMS.position=[2,2];
% pard.t3.object=struct('String','diameterNPC','Style','text');
% pard.t3.position=[3,1];
% pard.diameterNPC.object=struct('String','110','Style','edit');
% pard.diameterNPC.position=[3,2];
% pard.t4.object=struct('String','rim','Style','text');
% pard.t4.position=[4,1];
% pard.rim.object=struct('String','20','Style','edit');
% pard.rim.position=[4,2];
% 
% pard.saveon.object=struct('String','saveon','Style','checkbox');
% pard.saveon.position=[1,3];
% 
% pard.getmask.object=struct('String','getmask','Style','checkbox');
% pard.getmask.position=[2,3];
end


function cellLevelQC(obj,p)
    p.width_cellRimWidth = 150;

    se = obj.locData.SE;
    % Shaking
    allFileNumber = getFieldAsVector(obj.locData.files.file, 'number');
    for m = 1:length(allFileNumber)
        fileNumber = allFileNumber(m);
        locs = obj.locData.getloc({'xnm','ynm','groupindex','cell_DBSCAN'}, 'layer', 1, 'filenumber', fileNumber, 'position','all');

        file_cell = getFieldAsVector(se.cells, 'info.filenumber')';
        l_inTheCell = file_cell == fileNumber;
        cellSubset = se.cells(l_inTheCell);
        
        for c = 1:length(cellSubset)
            currentCell = cellSubset(c);
            lLocs = locs.cell_DBSCAN == currentCell.ID;
            xy = [locs.xnm(lLocs) locs.ynm(lLocs)];
            txy_bound = currentCell.evaluation.cellBoundary.txy;
            xy_bound = txy_bound(:,2:3);
            arcLen = arclength(xy_bound(:,1), xy_bound(:,2));
            [~, dist, t]= distance2curve(xy_bound,xy);
            l_kept = dist<=p.width_cellRimWidth;
            xy_insideCell = xy(l_kept,:);
            t_rim = t(l_kept)*arcLen;
            
            n = round(arcLen./10);
            
            theta_profile = histcounts(t_rim,n);
            r = rotAutoXorr(theta_profile,150);
            r = r(1:end-1);
            r = r-prctile(r(50:end),25);
            r = r./max(r);
            score = find(r<0.5,1, 'first');
            currentCell.evaluation.yeastCME_QC.score = score;
            if c == 1
                figure; plot(r)
            end
            hold on;
            plot(r)
            
        end
        allscore = getFieldAsVector(cellSubset, 'evaluation.yeastCME_QC.score')';
        figure; plot(allscore)
    end
    % Non-permeablized

end

function siteLevelQC(obj,p)

end

function r = rotAutoXorr(v, n)
    len = length(v);
    len = min([n len]);
    r = [];
    for k = 1:len
        r(k) = sum(v.*[v(end-k+1:end) v(1:end-k)]);
    end
end