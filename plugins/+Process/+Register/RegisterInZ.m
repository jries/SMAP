classdef RegisterInZ<interfaces.DialogProcessor
    properties
    end
    methods
        function obj=RegisterInZ(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_pixrec'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            out=[];
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function savebutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            [f,path]=uiputfile(fn,'Select transformation file _T.mat');
            if f
                obj.guihandles.Tfile.String=[path f];
                transformation=obj.transformation;
                save([path,f],'transformation');
            end      
        end 
    end
end




function pard=guidef(obj)

pard.texta.object=struct('String','target','Style','text');
pard.texta.position=[1,1];
pard.textb.object=struct('String','reference','Style','text');
pard.textb.position=[1,2];
pard.targetselect.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5');
pard.targetselect.position=[2,1];
pard.targetselect.load=false;
pard.refselect.object=struct('Style','popupmenu','String','layer1|layer2|layer3|layer4|layer5');
pard.refselect.position=[2,2];
pard.refselect.load=false;


pard.t1.object=struct('String','Width of stripes (nm)','Style','text');
pard.t1.position=[3,1];
pard.t1.Width=2;

pard.stripewidth.object=struct('Style','edit','String','150');
pard.stripewidth.position=[3,3];

pard.t2.object=struct('Style','text','String','Bin Width z (nm)');
pard.t2.position=[4,1];


pard.binwidth.object=struct('Style','edit','String','5');
pard.binwidth.position=[4,2];

pard.plugininfo.description=sprintf(['']);
pard.plugininfo.name='Register In Z';
pard.plugininfo.type='ProcessorPlugin';
end