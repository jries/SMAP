classdef Localizationnumber<interfaces.DialogProcessor&interfaces.SEProcessor
    properties
    end
    methods
        function obj=Localizationnumber(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            tic
            out=runintern(obj,p);
            toc
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

    end
end


function out=runintern(obj,p)
out=[];
tid=obj.locData.loc.tid;
time=obj.locData.loc.time;
file=obj.locData.loc.filenumber;
allid=unique(tid);
allfiles=unique(file);
locintrack=0*tid;
for f=1:length(allfiles)
    for id=1:length(allid)
        ixh=allfiles(f)==file & allid(id)==tid;
        numloc=sum(ixh);
        locintrack(ixh)=1:numloc;
    end
end
obj.locData.loc.locintrack=locintrack;
obj.locData.regroup;
end


function pard=guidef(obj)
pard.sourcet.object=struct('String','Calculate localization number in track','Style','text');
pard.sourcet.position=[1,1];
pard.sourcet.Width=4;
% pard.source.object=struct('String',{{'All filtered','ROI manager'}},'Style','popupmenu');
% pard.source.position=[1,2];
% 
% pard.skipt.object=struct('String','skip first','Style','text');
% pard.skipt.position=[3,1];
% pard.skip.object=struct('String','20','Style','edit');
% pard.skip.position=[3,2];
% 
% pard.minlt.object=struct('String','min length','Style','text');
% pard.minlt.position=[3,3];
% pard.minlen.object=struct('String','50','Style','edit');
% pard.minlen.position=[3,4];
% 
% 
% 
% pard.prect.object=struct('String','Convergence: ','Style','text');
% pard.prect.position=[4,1];
% pard.precmode.object=struct('String',{{'first d < ', 'last d > '}},'Style','popupmenu');
% pard.precmode.position=[4,2];
% pard.prec.object=struct('String','2','Style','edit');
% pard.prec.position=[4,3];
% pard.prec.Width=0.5;
% pard.prect2.object=struct('String','nm','Style','text');
% pard.prect2.position=[4,3.5];
% 
% 
% 
% pard.maxplottailst.object=struct('String','max number of plots','Style','text');
% pard.maxplottailst.position=[5,1];
% pard.maxplottails.object=struct('String','100','Style','edit');
% pard.maxplottails.position=[5,2];
% 
% pard.maxlt.object=struct('String','max length analysis','Style','text');
% pard.maxlt.position=[5,3];
% pard.maxlen.object=struct('String','50','Style','edit');
% pard.maxlen.position=[5,4];

% pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description=' ';

end