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
          
          locs=obj.locData.getloc({f1a,f1b,f2a,f2b,'frame',[f1a 'err'],[f2a 'err'],[f1b 'err'],[f2b 'err']},'layer',find(obj.getPar('sr_layerson')),'position','roi');
          indgood=~isnan(locs.(f1a)) & ~isnan(locs.(f1b)) & ~isnan(locs.(f2a)) & ~isnan(locs.(f2b)) ;
          da=locs.(f1a)(indgood)-locs.(f2a)(indgood);
          db=locs.(f1b)(indgood)-locs.(f2b)(indgood);
          ax=obj.initaxis('drift vs frame');
          fmax=max(locs.frame(indgood));
          df=10^floor(log10(fmax)-2);
          ff=0:df:fmax;
          fff=0:df*10:fmax;
          mode=p.mode.selection;
          

          if strcmp(mode,'weighted mean')
              wx=1./(locs.([f1a 'err'])(indgood).^2+locs.([f2a 'err'])(indgood).^2);  %weights: It is variance! see Wikipedia
              wy=1./(locs.([f1b 'err'])(indgood).^2+locs.([f2b 'err'])(indgood).^2);
              [fita,dab]=getsmoothcurve(locs.frame(indgood),da,wx,ff');
              [fitb,dbb]=getsmoothcurve(locs.frame(indgood),db,wy,ff');
              fff=ff;
              dabf=fita(ff);
              dbbf=fitb(ff);
              
          else
              dab=bindata(locs.frame(indgood),da,ff,mode);
              dbb=bindata(locs.frame(indgood),db,ff,mode);
              dabf=bindata(locs.frame(indgood),da,fff,mode);
              dbbf=bindata(locs.frame(indgood),db,fff,mode);
          end
          
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

function [fitx,dxb]=getsmoothcurve(frames,dx,wx,ff)      
if nargin<4
    dframe=100;
    ff=(min(frames): dframe:max(frames))';
else
    dframe=ff(2)-ff(1);
end

if nargin<3 || isempty(wx)
    wx=[];
end
[dxb,sb]=bindatamean(frames,dx,ff,wx);
[fitx,gof,p]=fit(ff,dxb,'smoothingspline','Weights',sb,'SmoothingParam',1e-11*dframe);
end
            


function pard=guidef
pard.textx.object=struct('String','xfield','Style','text');
pard.textx.position=[2,1];
pard.texty.object=struct('String','yfield','Style','text');
pard.texty.position=[3,1];
pard.textr.object=struct('String','ch1','Style','text');
pard.textr.position=[1,2];
pard.textt.object=struct('String','ch2','Style','text');
pard.textt.position=[1,3];

pard.assignfield1a.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1a.position=[2,2];
pard.assignfield2a.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield2a.position=[2,3];
pard.assignfield1b.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1b.position=[3,2];
pard.assignfield2b.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield2b.position=[3,3];
pard.mode.object=struct('Style','popupmenu','String',{{'weighted mean','mean','median'}});
pard.mode.position=[1,4];

pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description='Calculates the  difference in x and y between two chanels and plots it vs the frame to test for drift between channels over time.';
end