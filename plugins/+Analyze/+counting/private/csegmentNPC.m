classdef csegmentNPC<recgui.DialogProcessor
    methods
        function obj=csegmentNPC(varargin)        
                obj@recgui.DialogProcessor;
%                 if ~isempty(parameter)
%                     obj.setParameters(parameter)
%                 end
            if nargin>0
                obj.handle=varargin{1};
            
%             if nargin>1
%                 obj.makeGui(guidef);
            end
%             end   

        end
        
        function out=run(obj,p)
           
            locs=obj.locData.getloc({'frame','xnm','ynm','phot','bg','PSFxnm','locprecnm'},'layer',1,'position','roi');

            obj.locData.guiData.counting.clusters=cluster_segment_NPC(locs,p);
%             obj.setGuiParameters(pout);
%             out=0;
        end
        
        function refit_callback(obj)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        

    end
    methods(Static)
        function info=info(obj)
            info.name='2b. segment NPC';
            info.class=@csegmentNPC;
            info.tag='csegmentNPC';
%             obj.info=info;
        end

    end
end




function pard=guidef(obj)
% 
pard.text2.object=struct('String','length scale nm','Style','text');
pard.text2.position=[1,1];

pard.lengthscale.object=struct('String','15','Style','edit');
pard.lengthscale.position=[1,2];
pard.lengthscale.isnumeric=1;


pard.text3.object=struct('String','maxPSF focus','Style','text');
pard.text3.position=[2,1];
pard.segment_maxPSF.object=struct('String','130','Style','edit');
pard.segment_maxPSF.position=[2,2];
pard.segment_maxPSF.isnumeric=1;
% pard.text1.object=struct('String','parameters','Style','text');
% pard.text1.position=[1,1];
% 
% pard.N0_fit.object=struct('String','N0','Style','radiobutton');
% pard.N0_fit.position=[2,2];
% 
% pard.N0_v.object=struct('String','10','Style','edit');
% pard.N0_v.position=[2,3];
% pard.N0_v.isnumeric=1;
% 
% 
% pard.pmature_fit.object=struct('String','p mature','Style','radiobutton');
% pard.pmature_fit.position=[3,2];
% 
% pard.pmature_v.object=struct('String','.5','Style','edit');
% pard.pmature_v.position=[3,3];
% pard.pmature_v.isnumeric=1;
% 
% 
% pard.pblink_fit.object=struct('String','p blink','Style','radiobutton');
% pard.pblink_fit.position=[4,2];
% 
% pard.pblink_v.object=struct('String','.2','Style','edit');
% pard.pblink_v.position=[4,3];
% pard.pblink_v.isnumeric=1;
% 
% 
% pard.monomer_fit.object=struct('String','monomer fraction','Style','radiobutton');
% pard.monomer_fit.position=[5,2];
% 
% pard.monomer_v.object=struct('String','.2','Style','edit');
% pard.monomer_v.position=[5,3];
% pard.monomer_v.isnumeric=1;
% 
% 
% 
% 
% pard.ncluster_fit.object=struct('String','n in cluster','Style','radiobutton');
% pard.ncluster_fit.position=[6,2];
% 
% pard.ncluster_v.object=struct('String','0','Style','edit');
% pard.ncluster_v.position=[6,3];
% pard.ncluster_v.isnumeric=1;
% 
% pard.blinkmode.object=struct('Style','popupmenu','String','Poisson|Exponential|Noblink');
% pard.blinkmode.position=[7,2];

% 
% pard.results1='results 1';
% pard.results3='other results';
% 
% pard.refit.object=struct('String','length scale nm','Style','pushbutton','Callback',{@obj.refit_callback});
% pard.refit.position=[4,4];

end