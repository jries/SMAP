classdef SALMInt2z<interfaces.DialogProcessor
    methods
        function obj=SALMInt2z(varargin)  
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            f1=p.assignfield1.selection;
            f2=p.assignfield2.selection;
            
            i1=obj.locData.loc.(f1);
            i2=obj.locData.loc.(f2);
%             
%             r5050=1.146;
%             b5050=1/(r5050+1);a5050=1-b5050;
%             aBS=a5050;bBS=b5050;
%             cal.a=0.657;
%             cal.b=4.116;
%             cal.c=-0.2024;
%             IsIu=1./(i2/aBS./(i1/bBS-i1/aBS));
            
            
            cal.a=0.5483;
            cal.b=5.054;
            cal.c=0.7466;
            IsIu=i1./i2;

            z=(-1/cal.b*log((IsIu-cal.c)/cal.a))*1000;
            figure(23)
            hist(real(z),500)
            z(isnan(z))=-2;
            z(imag(z)~=0)=-2;
            
            obj.locData.loc.znm=z;
            
            obj.locData.regroup;
%             
        end
        function pard=guidef(obj)
            pard=guidef;
        end
%         function attachLocData(obj,locData)
%             attachLocData@recgui.GuiProcessor(obj,locData);
% %             addlistener(obj.locData,'synchronizeGui',@obj.synchronizeGui);
%             addlistener(obj.locData,'loaded',@obj.updateGui);
%             addlistener(obj.locData,'updateGui',@obj.updateGui);
%         end
        function updateGui(obj,event,data)
            if isfield(obj.locData.files,'file')
            ff=obj.locData.files.file;
            str={};
            for k=1:length(ff)
                if 0 %ff(k).istiff
                    str{end+1}=' ';
                else
                    str{end+1}=ff(k).name;
                end
            end
            obj.guihandles.dataselect.String=str; 
            obj.guihandles.assignmode.Callback={@changemode,obj};
            changemode(0,0,obj);
            end
        end
    end
end


function changemode(a,b,obj)
excludefields={'frame','x','xnm','y','ynm','z','znm','locprecxnm','locprecynm',...
    'locprecznm','channel','PSFxnm','PSFynm','peakfindxnm','peakfindynm','filenumber',...
    'groupindex','numberInGroup','locprecnm'};
if ~isempty(obj.locData.loc)
fnpresent=fieldnames(obj.locData.loc);
showfields=setdiff(fnpresent,excludefields);
obj.guihandles.assignfield1.String=showfields;
obj.guihandles.assignfield2.String=showfields;
end
end

function pard=guidef
pard.texta.object=struct('String','dataset','Style','text');
pard.texta.position=[1,1];

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];

pard.dataselect.object.TooltipString='choose localization file data set';

pard.textb.object=struct('String','slope','Style','text');
pard.textb.position=[1,3];
pard.split_slope.object=struct('Style','edit','String','.5 .51');
pard.split_slope.position=[1,4];

pard.split_slope.object.TooltipString='slope to split. give two numbers for upper and lower slope, data in between is excluded.';
pard.textb.object.TooltipString=pard.split_slope.object.TooltipString;


pard.textc.object=struct('String','edge [1 2]','Style','text');
pard.textc.position=[2,3];
pard.split_edge.object=struct('Style','edit','String','0 0');
pard.split_edge.position=[2,4];
pard.split_edge.object.TooltipString='additional part around slopes to be excluded. One number or two for asymmetric edge.';
pard.textc.object.TooltipString=pard.split_edge.object.TooltipString;



pard.textd.object=struct('String','int min [1 2]','Style','text');
pard.textd.position=[3,3];
pard.split_intmin.object=struct('Style','edit','String','0 0');
pard.split_intmin.position=[3,4];

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


pard.assignfield1.object.TooltipString='choose which field to use for splitting';
pard.assignfield2.object.TooltipString=pard.assignfield1.object.TooltipString;
pard.plugininfo.type='ProcessorPlugin';
end