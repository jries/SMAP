classdef QuantifyChannelDrift<interfaces.DialogProcessor
%     Calculates the  difference in x and y between two chanels and plots
%     it vs the frame to test for drift between channels over time.
    properties
        isz=0;
        transformation=[];
        register_parameters;
    end
    methods
        function obj=QuantifyChannelDrift(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.showresults=true;
%             obj.history=true;
        end
        
        function out=run(obj,p)
          f1a=p.assignfield1a.selection;
          f2a=p.assignfield2a.selection;
          f1b=p.assignfield1b.selection;
          f2b=p.assignfield2b.selection;
          
          locs=obj.locData.getloc({f1a,f1b,f2a,f2b,'frame'},'layer',find(obj.getPar('sr_layerson')),'position','roi');
          da=locs.(f1a)-locs.(f2a);
          db=locs.(f1b)-locs.(f2b);
          ax=obj.initaxis('drift vs frame');
          fmax=max(locs.frame);
          df=10^floor(log10(fmax)-2);
          ff=0:df:fmax;
          fff=0:df*10:fmax;
          mode=p.mode.selection;
          dab=bindata(locs.frame,da,ff,mode);
          dbb=bindata(locs.frame,db,ff,mode);
          dabf=bindata(locs.frame,da,fff,mode);
          dbbf=bindata(locs.frame,db,fff,mode);
          
          plot(ax,ff,dab,'b',ff,dbb,'m',fff,dabf,'k',fff,dbbf,'r',fff,mean(dabf)+0*fff,'b',fff,mean(dbbf)+0*fff,'r',fff,0*fff,'k')
          xlabel('frame')
          ylabel([f1a ', ' f1b])
          legend(f1a,f1b)
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
%         function makeGui(obj)
%             makeGui@interfaces.DialogProcessor(obj);
%         end
        function initGui(obj)
            obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
            obj.updateLocFields;
        end
        function updateLocFields(obj)
            excludefields={'frame','channel','peakfindxnm','peakfindynm','filenumber',...
                'groupindex','numberInGroup'};
            if ~isempty(obj.locData.loc)
                fnpresent=fieldnames(obj.locData.loc);
                showfields=setdiff(fnpresent,excludefields);
                obj.guihandles.assignfield1a.String=showfields;
                obj.guihandles.assignfield2a.String=showfields;
                obj.guihandles.assignfield1a.Value=min(obj.guihandles.assignfield1a.Value,length(obj.guihandles.assignfield1a.String));
                obj.guihandles.assignfield2a.Value=min(obj.guihandles.assignfield2a.Value,length(obj.guihandles.assignfield2a.String));
                obj.guihandles.assignfield1b.String=showfields;
                obj.guihandles.assignfield2b.String=showfields;
                obj.guihandles.assignfield1b.Value=min(obj.guihandles.assignfield1a.Value,length(obj.guihandles.assignfield1a.String));
                obj.guihandles.assignfield2b.Value=min(obj.guihandles.assignfield2a.Value,length(obj.guihandles.assignfield2a.String));

            end
        end
    end
end



function pard=guidef
pard.texta.object=struct('String','target:','Style','text');
pard.texta.position=[2,2.05];
pard.assignfield1a.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1a.position=[2,1];
pard.assignfield2a.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield2a.position=[2,2];
pard.assignfield1b.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1b.position=[3,1];
pard.assignfield2b.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield2b.position=[3,2];
pard.mode.object=struct('Style','popupmenu','String',{{'mean','median'}});
pard.mode.position=[1,3];

pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description='Calculates the  difference in x and y between two chanels and plots it vs the frame to test for drift between channels over time.';
end