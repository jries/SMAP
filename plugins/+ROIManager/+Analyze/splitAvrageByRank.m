classdef splitAvrageByRank<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=splitAvrageByRank(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer', 'se_siteroi', 'se_sitefov'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            %% Initialization
            % this should be added as a plugin of SMAP later
            betweenAvgDistance = p.betweenBins;       % distance between bins
            windowSize = p.windowSize;                % # of sites that form one bin
            gapBetweenAvg = p.siteRange(2);            % # of sites between the middle sites in neighbouring bins
            lastSite = p.siteRange(3);
            startingSite = p.siteRange(1);            % starting sites
            offset = p.fromStack;                  % offset (distance) from the original average
            trim = p.trim;
            
            %% Main
            % back up the original pos of locs
            if ~isfield(obj.locData.loc,'xnmori')
                obj.locData.loc.xnmori = obj.locData.loc.xnm;
                obj.locData.loc.ynmori = obj.locData.loc.ynm;
            end
            % get the pos from the backup
            obj.locData.loc.xnm = obj.locData.loc.xnmori;
            obj.locData.loc.ynm = obj.locData.loc.ynmori;
            
            lMiddlePart_x = obj.locData.loc.xnm>=p.se_siteroi-p.se_siteroi/2+trim(1)&obj.locData.loc.xnm<=(p.se_siteroi+p.se_siteroi/2-trim(1));
            lMiddlePart_y = obj.locData.loc.ynm>=p.se_siteroi-p.se_siteroi/2+trim(2)&obj.locData.loc.ynm<=(p.se_siteroi+p.se_siteroi/2-trim(2));
            lMiddlePart = lMiddlePart_x&lMiddlePart_y;
            % get the order of sites based on the previous sorting
            siteOrder = obj.locData.loc.(p.rankOption.selection);

            % the site number of the centers of bins
            centralValue = (startingSite:gapBetweenAvg:lastSite);
            if isEven(windowSize)
                windowLb = round(windowSize/2);
                windowUb = round(windowSize/2);
            else
                windowLb = round(windowSize/2);
                windowUb = round(windowSize/2)-1;
            end
            centralValue = centralValue+windowLb;
            centralValue = centralValue(centralValue<lastSite);
            % shift the bins
            for k = length(centralValue):-1:1
                dx = k*betweenAvgDistance+offset;
                lb = centralValue(k)-windowLb;
                ub = centralValue(k)+windowUb;
                lselectedLocs = siteOrder>=lb&siteOrder<ub;
                lShift = lselectedLocs&lMiddlePart;
                obj.locData.loc.xnm(lShift) = obj.locData.loc.xnmori(lShift) + dx;
            end
            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function recover_callback(a,b,obj)
    % get the pos from the backup
    obj.locData.loc.xnm = obj.locData.loc.xnmori;
    obj.locData.loc.ynm = obj.locData.loc.ynmori;
end

function rank_callback(a,b,obj)
    % get the pos from the backup
    fn = fieldnames(obj.locData.loc);
    a.String = fn';
    a.Enable = 'on';
end


function pard=guidef(obj)
pard.t1.object=struct('String','Between bins (nm)','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.betweenBins.object=struct('Style','edit', 'String','240');
pard.betweenBins.position=[1,2];
pard.betweenBins.Width=1;

pard.t2.object=struct('String','Range','Style','text');
pard.t2.position=[2,1];
pard.t2.Width=1;

pard.siteRange.object=struct('String','Start binSize end','Style','edit');
pard.siteRange.position=[2,2];
pard.siteRange.Width=1;

pard.t3.object=struct('String','Window size','Style','text');
pard.t3.position=[3,1];
pard.t3.Width=1;

pard.windowSize.object=struct('Style','edit', 'String','10');
pard.windowSize.position=[3,2];
pard.windowSize.Width=1;

pard.t4.object=struct('String','From original stack','Style','text');
pard.t4.position=[4,1];
pard.t4.Width=1;

pard.fromStack.object=struct('Style','edit', 'String','1000');
pard.fromStack.position=[4,2];
pard.fromStack.Width=1;

pard.tTrim.object=struct('String','Trim (x and y)','Style','text');
pard.tTrim.position=[5,1];
pard.tTrim.Width=1;

pard.trim.object=struct('Style','edit', 'String','x y');
pard.trim.position=[5,2];
pard.trim.Width=1;

pard.t5.object=struct('String','Rank','Style','text');
pard.t5.position=[6,1];
pard.t5.Width=1;

fn = fieldnames(obj.locData.loc);
pard.rankOption.object=struct('Style', 'popupmenu', 'String',{fn'},'value',1);
pard.rankOption.position=[6,2];
pard.rankOption.Width=1;
pard.rankOption.ButtonDownFcn = {@rank_callback,obj};
pard.rankOption.Enable = 'Inactive';

pard.recover.object=struct('Style', 'pushbutton', 'String', 'Recover','Callback',{{@recover_callback, obj}});
pard.recover.position=[7,1];
pard.recover.Width=1;

pard.plugininfo.type='ROI_Analyze';

end