classdef FilterField<interfaces.DialogProcessor&interfaces.SEProcessor
%    Filters a localization attribute. This is mainly used for the BatchAnalysis
%     plugin.
    methods
        function obj=FilterField(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)  
            out=[];
            sf={p.field,p.fieldmin,p.fieldmax, p.fieldfilter};
            l1=find(obj.getPar('sr_layerson'));
            for k=1:length(l1)
                obj.setPar('selectedField',sf,'layer',l1(k))
            end
%             obj.locData.filter;
            obj.locData.regroup;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function fieldchange_callback(a,b,obj)
field=obj.getSingleGuiParameter('fieldselect').selection;
obj.setGuiParameters(struct('field',field));
l1=find(obj.getPar('sr_layerson'));
ft=obj.getPar('filtertable','layer',l1(1));
ind=strcmp(ft(:,1),field);
obj.setGuiParameters(struct('fieldmin',ft{ind,2}))
obj.setGuiParameters(struct('fieldmax',ft{ind,6}))
obj.setGuiParameters(struct('fieldfilter',ft{ind,7}))
end


function pard=guidef(obj)

pard.t1.object=struct('String','Filter a field. Mostly used for batch analysis.','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=1;

pard.t1m.object=struct('String','min','Style','text');
pard.t1m.position=[1,3];
pard.t1m.Width=.5;
pard.t1mm.object=struct('String','max','Style','text');
pard.t1mm.position=[1,3.5];
pard.t1mm.Width=.5;

pard.fieldt.object=struct('String','Field','Style','text');
pard.fieldt.position=[2,1];
pard.fieldt.Width=1;

pard.field.object=struct('String','','Style','edit');
pard.field.position=[2,2];
pard.field.Width=1;

pard.fieldselect.object=struct('String','','Style','popupmenu','Callback',{{@fieldchange_callback,obj}});
pard.fieldselect.position=[3,2];
pard.fieldselect.Width=1;

pard.fieldmin.object=struct('String','0','Style','edit');
pard.fieldmin.position=[2,3];
pard.fieldmin.Width=.5;
pard.fieldmax.object=struct('String','0','Style','edit');
pard.fieldmax.position=[2,3.5];
pard.fieldmax.Width=.5;

pard.fieldfilter.object=struct('String','filter','Style','checkbox','Value',true);
pard.fieldfilter.position=[2,4];
pard.fieldfilter.Width=1;

pard.syncParameters={{'locFields','fieldselect',{'String'}}}; %,{@fieldchange_callback,obj}

pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description=' Filters a localization attribute. This is mainly used for the BatchAnalysis plugin.';
end