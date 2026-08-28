%% This pluging finds trajectories and list to the ROI manager. Modified from Segment_cotracks.m by Lasse Drengk.
%% additional functions: duplicate filter, skipping the first n locs when determining the length

classdef Segment_singletracks<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    methods
        function obj = Segment_singletracks(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults = true;
        end
        function out = run(obj,p)
            
            out = runi(obj,p);
        end
        function pard = guidef(obj)
            pard = guidef(obj);
        end
       
    end
end

function out = runi(obj,p)
out = [];

% Parameters
minlenlocs = p.minlenlocs;
lennmstartend = p.lennmstartend;
skipN = p.skipNLocs;
layers = find(obj.getPar('sr_layerson'));

% Get localization data
[locs,~] = obj.locData.getloc({'xnm','ynm','time','tid','filenumber','thi'},'layer',layers,'position','all','grouping','ungrouped','removefilter','filenumber');

% Remove duplicates safely
locMat = [locs.xnm(:), locs.ynm(:), locs.time(:), locs.tid(:), locs.filenumber(:), locs.thi(:)];
[~, ia] = unique(locMat, 'rows', 'stable'); % Keep first occurrences
fields = fieldnames(locs);
for f = 1:numel(fields)
    locs.(fields{f}) = locs.(fields{f})(ia);
end
% end of removing duplicates. This code part may be removed if wished

SE = obj.locData.SE;

% Initialize cell grid if empty
if isempty(SE.cells)
    cg = ROIManager.Segment.makeCellGrid;
    cg.attachLocData(obj.locData);
    cg.attachPar(obj.P);
    p = cg.getAllParameters;
    cg.run(p);
end

cellp = vertcat(SE.cells.pos);

for c = length(SE.cells):-1:1
    cellfn(c,1) = SE.cells(c).info.filenumber;
end

fnums = unique(locs.filenumber);

for ff = 1:length(fnums)
    fhh = fnums(ff);

    % get tracks in the channel of interest, e.g. thi==0
    tids0 = locs.tid(locs.thi == 0 & locs.filenumber == fhh);
    numloc = histcounts(tids0,1:max(tids0)+1);
    id0 = find(numloc > minlenlocs);

    for k = 1:length(id0)
        idh = locs.tid==id0(k) & locs.filenumber == fhh;
        xh = locs.xnm(idh); 
        yh = locs.ynm(idh);
        time = locs.time(idh);

        % Filtering out the first N locs for lenght of track
        xh = xh((skipN+1):end);
        yh = yh((skipN+1):end);
        time = time((skipN+1):end);

        len = sqrt((xh(end)-xh(1)).^2+(yh(end)-yh(1)).^2);
        if len < lennmstartend
            continue
        end

        currentsite = interfaces.SEsites;
        currentsite.pos = [mean(xh), mean(yh), 0];  

        thisf = ~(cellfn == fhh);
        [~,cind] = min(sum((currentsite.pos(1:2)-cellp).^2,2)+thisf*1e9);

        currentsite.info.cell = SE.cells(cind).ID;
        currentsite.info.filenumber = fhh;
        currentsite.annotation.comments = num2str(time([1 end])/1000);

        SE.addSite(currentsite);
    end
end

SE.processors.preview.updateSitelist;
end
             
function pard = guidef(obj)

pard.minlenlocst.object = struct('String','Min length track (locs)','Style','text');
pard.minlenlocst.position = [2,1];
pard.minlenlocst.Width = 1.5;

pard.minlenlocs.object = struct('String','100','Style','edit');
pard.minlenlocs.position = [2,2.5];
pard.minlenlocs.Width = .5;

pard.lennmstartendt.object = struct('String','start-end (nm) >','Style','text');
pard.lennmstartendt.position = [2,3.5];
pard.lennmstartendt.Width = 1.;

pard.lennmstartend.object = struct('String','50','Style','edit');
pard.lennmstartend.position = [2,4.5];
pard.lennmstartend.Width = .5;

pard.skipNLocst.object = struct('String','Skip first N locs','Style','text');
pard.skipNLocst.position = [3,1];
pard.skipNLocst.Width = 1;

pard.skipNLocs.object = struct('String','0','Style','edit');
pard.skipNLocs.position = [3,2.5];
pard.skipNLocs.Width = 0.5;

pard.plugininfo.description = sprintf('single-path tracking analysis');
pard.plugininfo.type = 'ProcessorPlugin';
end