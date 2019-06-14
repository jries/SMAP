classdef RemoveDoubleSites<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=RemoveDoubleSites(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)  
            out=[];
            sites=obj.locData.SE.sites;
            x=getFieldAsVectorInd(sites,'pos',1);
            y=getFieldAsVectorInd(sites,'pos',2);
            ind=1;
            indices=1:length(x);
            while ind<=length(indices)
                xyh=sites(indices(ind)).pos;
                d2=(x(indices)-xyh(1)).^2+(y(indices)-xyh(2)).^2;
                indb= d2<p.mind^2 & d2>0;
                indices(indb)=[];
                ind=ind+1;
            end
            obj.locData.SE.sites=sites(indices);
            
          
        
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.t1.object=struct('String','min distance (nm): ','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1.5;

pard.mind.object=struct('String','500','Style','edit');
pard.mind.position=[1,2.5];
pard.mind.Width=0.5;


pard.plugininfo.type='ROI_Analyze';


end