classdef intensity_ratios_channels<interfaces.DialogProcessor
% Calibrates field-dependent intensity ratio for ratiometric dual-camera
% imaging. Used for salvaged fluorescence.
    properties
    end
    methods
        function obj=intensity_ratios_channels(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;  
            obj.showresults=true;
        end
        function out=run(obj,p)
            out=[];
            dc=0.01;
            df=0.002;
            window=2;
            if obj.processorgui==false || p.assignfield1.Value==1%run from WF
                p.assignfield1.selection='fit_n2';
                p.assignfield2.selection='fit_n1';
            end
            field1=p.assignfield1.selection;
            field2=p.assignfield2.selection;
            
            locs=obj.locData.getloc({'xnm','ynm',field1,field2},'layer',find(obj.getPar('sr_layeron')),'Position','roi');
            ratio=locs.(field2)./locs.(field1);
            g=(ratio-1)./(ratio+1);
            n=-1:dc:1;
            weights=(locs.(field1)+locs.(field2));
            weights=weights/mean(weights);
            [wh,nc]=weightedhistcounts(g,weights,n);
            ax=obj.initaxis('ratio');
            hold(ax,'off')
            plot(ax,nc,wh)
            hold(ax,'on')
            ind=find(islocalmax(wh));
            [maxsort,inds]=sort(wh(ind),'descend');
            peakpos=nc(ind(inds(1:p.numcolors)));
            gpos=zeros(1,length(peakpos));
            for k=1:length(peakpos)
                nh=peakpos(k)-window*dc:df:peakpos(k)+window*dc;
                [hih,nch]=weightedhistcounts(g,weights,nh);
                plot(ax,nch,hih*dc/df,'g');
                fitp=fit(nch',hih'*dc/df,'poly2');
                plot(ax,nch,fitp(nch),'r')
                gpos(k)=-fitp.p2/2/fitp.p1;
                plot(ax,[gpos(k) gpos(k)],[0 fitp(gpos(k))],'k')
            end
            rpos=(1+gpos)./(1-gpos);
            rpos=sort(rpos,'descend');
            title(num2str(rpos))
            clipboard('copy',rpos)
%             histogram(ax,g)
%             bad=isnan(ratio)|isinf(ratio) | ratio<0;
%             xnm=locs.xnm(~bad);ynm=locs.ynm(~bad);ratio=ratio(~bad);
%             ax=obj.initaxis('ratio');
%             plot(ax,xnm,ratio,'.')
%             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
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
                setdefaultfields(0,0,obj)
            end
        end

    end
end

function setdefaultfields(a,b,obj)
f1='phot1';
f2='phot2';
fs1=obj.getSingleGuiParameter('assignfield1');
fs2=obj.getSingleGuiParameter('assignfield2');

ind1=find(strcmp(fs1.String,f1),1,'first');
if ~isempty(ind1)
    fs1.Value=ind1;
    fs1.selection=f1;
end
ind2=find(strcmp(fs2.String,f2),1,'first');
if ~isempty(ind2)
    fs2.Value=ind2;
    fs2.selection=f2;
end
obj.setGuiParameters(struct('assignfield1',fs1,'assignfield2',fs2));
end

function pard=guidef(obj)
pard.t1.object=struct('String','fields','Style','text');
pard.t1.position=[1,1];

pard.assignfield1.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield1.position=[2,1];
pard.assignfield2.object=struct('Style','popupmenu','String','n1|n2');
pard.assignfield2.position=[3,1];
pard.setdefault.object=struct('String','default','Style','pushbutton','Callback',{{@setdefaultfields,obj}});
pard.setdefault.position=[2,2];
pard.setdefault.Width=0.6;


pard.numcolorst.object=struct('String','Number of colors','Style','text');
pard.numcolorst.position=[1,3];
pard.numcolorst.Width=1;
pard.numcolors.object=struct('String','2','Style','edit');
pard.numcolors.position=[1,4];

pard.assignfield1.object.TooltipString='choose which field to use for splitting';
pard.assignfield2.object.TooltipString=pard.assignfield1.object.TooltipString;


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Calibrates intensity ratios among channels for global multi-channel fitting.';
end

function [out,nc]=weightedhistcounts(in,weights,n)
[ins,inds]=sort(in);
ws=weights(inds);
out=zeros(1,length(n)-1);
for k=1:length(n)-1
    indh=ins>n(k) & ins <= n(k+1);
    out(k)=sum(ws(indh));
end
nc=n(1:end-1)+(n(2)-n(1))/2;

end