classdef Intensity2Channel<interfaces.DialogProcessor
    methods
        function obj=Intensity2Channel(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;  
            obj.history=false;
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            out=[];
            if obj.processorgui==false || p.assignfield1.Value==1%run from WF
                p.assignfield1.selection='fit_n2';
                p.assignfield2.selection='fit_n1';
            end
            [ll,ind]=obj.locData.getloc({'xnm','inungrouped'},'Position','roi','layer',find(obj.getPar('sr_layerson')),'removeFilter','channel','grouping','ungrouped');
%             loc=get_intensity2ch(obj.locData.loc,p,ll.inungrouped);
            loc=get_intensity2ch(obj.locData.loc,p,ind);
            if p.combineunassigned
                loc.channel(loc.channel==3)=1;
                loc.channel(loc.channel==4)=2;
            end
            obj.locData.loc=copyfields(obj.locData.loc,loc);
            obj.locData.regroup;
%             
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
            obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
            obj.updateLocFields;
        end

        function updateLocFields(obj)
            excludefields={'frame','x','xnm','y','ynm','z','znm','locprecxnm','locprecynm',...
                'locprecznm','channel','PSFxnm','PSFynm','peakfindxnm','peakfindynm','filenumber',...
                'groupindex','numberInGroup','locprecnm'};
            if ~isempty(obj.locData.loc)
            fnpresent=fieldnames(obj.locData.loc);
            showfields=setdiff(fnpresent,excludefields);
            obj.guihandles.assignfield1.String=showfields;
            obj.guihandles.assignfield2.String=showfields;
            obj.guihandles.assignfield1.Value=min(obj.guihandles.assignfield1.Value,length(obj.guihandles.assignfield1.String));
            obj.guihandles.assignfield2.Value=min(obj.guihandles.assignfield2.Value,length(obj.guihandles.assignfield2.String));
            end
        end
    end
end




function pard=guidef
pard.texta.object=struct('String','dataset','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];

pard.dataselect.object.TooltipString='choose localization file data set';

pard.textb.object=struct('String','slope','Style','text');
pard.textb.position=[3,3];
pard.split_slope.object=struct('Style','edit','String','1');
pard.split_slope.position=[3,4];

pard.split_slope.object.TooltipString='slope to split. give two numbers for upper and lower slope, data in between is excluded.';
pard.textb.object.TooltipString=pard.split_slope.object.TooltipString;


pard.textc.object=struct('String','edge [1 2]','Style','text');
pard.textc.position=[2,3];
pard.split_edge.object=struct('Style','edit','String','.5');
pard.split_edge.position=[2,4];
pard.split_edge.object.TooltipString='additional part around slopes to be excluded. One number or two for asymmetric edge.';
pard.textc.object.TooltipString=pard.split_edge.object.TooltipString;

pard.texto.object=struct('String','offset (int 1)','Style','text');
pard.texto.position=[1,3];
pard.split_offset.object=struct('Style','edit','String','0 0');
pard.split_offset.position=[1,4];
pard.split_offset.object.TooltipString='Offset for split line.';
pard.texto.object.TooltipString=pard.split_offset.object.TooltipString;


pard.textd.object=struct('String','int min [1 2]','Style','text');
pard.textd.position=[4,3];
pard.split_intmin.object=struct('Style','edit','String','0 0');
pard.split_intmin.position=[4,4];

pard.split_intmin.object.TooltipString='minimum intensity required for split.';
pard.textd.object.TooltipString=pard.split_intmin.object.TooltipString;


% pard.assignmode.object=struct('Style','popupmenu','String','field|N/A: two rois');
% pard.assignmode.position=[3,1];
% pard.assignmode.object.TooltipString='NA';
% 


pard.t1.object=struct('String','value 1','Style','text');
pard.t1.position=[4,1];
pard.t2.object=struct('String','value 2','Style','text');
pard.t2.position=[4,2];

pard.assignfield1.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1.position=[5,1];
pard.assignfield2.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield2.position=[5,2];

pard.logscale.object=struct('Style','checkbox','String','log scale','Value',1);
pard.logscale.position=[5,4];



pard.assignfield1.object.TooltipString='choose which field to use for splitting';
pard.assignfield2.object.TooltipString=pard.assignfield1.object.TooltipString;

pard.combineunassigned.object=struct('String','Associate unassiged locs','Style','checkbox','Value',true);
pard.combineunassigned.position=[7,1];
pard.combineunassigned.Width=2;
pard.plugininfo.type='ProcessorPlugin';
end