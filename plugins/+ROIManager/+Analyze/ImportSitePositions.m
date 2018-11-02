classdef ImportSitePositions<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=ImportSitePositions(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'se_layerson'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
            out=[];
           l=load(p.oldfile);
           se=l.saveloc.siteexplorer;
           obj.locData.SE=copyfields(obj.locData.SE,se);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function load_callback(obj,a,b)
            f=obj.getSingleGuiParameter('oldfile');
            p=fileparts(f);
            [f,path]=uigetfile([p filesep '*.mat']);
            if f                
                obj.setGuiParameters(struct('oldfile',[path f]));
            end
                
        end
    end
end




function pard=guidef(obj)

pard.t1.object=struct('String','Import sites and cells from other file. ','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;
pard.t1.object.TooltipString='If multiple files were loaded, make sure the corresponding ones are loaded here and have the same order';


pard.oldfile.object=struct('String','','Style','edit');
pard.oldfile.position=[2,1];
pard.oldfile.Width=3;


pard.loadfile.object=struct('String','load','Style','pushbutton','Callback',@obj.load_callback);
pard.loadfile.position=[2,4];
pard.loadfile.Width=1;



pard.plugininfo.type='ROI_Analyze';


end