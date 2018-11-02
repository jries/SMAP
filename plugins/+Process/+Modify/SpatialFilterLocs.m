classdef SpatialFilterLocs<interfaces.DialogProcessor
    methods
        function obj=SpatialFilterLocs(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','SpatialFilterLocs');
            notify(obj.P,'backup4undo');
            if p.dataselectall
                filenumber=unique(obj.locData.loc.filenumber);
            else
                filenumber=p.dataselect.Value;
            end
%             locs = obj.locData.loc;
            for fi=1:length(filenumber)
            [locs,indinb]=obj.locData.getloc({'xnm','ynm',p.locfield.selection},'layer',1,'position','all','filenumber',filenumber(fi),'grouping','ungrouped');
%             indin=find(locs.filenumber==p.dataselect.Value);
            indin=find(indinb);
            [xnm,sortind]=sort(locs.xnm);
            ynm=locs.ynm(sortind);
            value=locs.(p.locfield.selection)(sortind);
%                         [xnm,sortind]=sort(locs.xnm(indin));
%             ynm=locs.ynm(indin(sortind));
%             value=locs.(p.locfield.selection)(indin(sortind));
            
            dist=p.spatialscale;
            ind1=1;
            ind2=1;
            outval=zeros(size(indinb));
            
            switch p.filter.selection
                case 'mean'
                    locfilter=@mean;
                case 'median'
                    locfilter=@median;
                case 'min'
                    locfilter=@min;
                case 'max'
                    locfilter=@max;
                case 'std'
                    locfilter=@std;
                case '5-95% percentile'
                    locfilter=@perc9_95;
            end
            lx=length(sortind);
            
            timerx=tic;
            for k=1:lx
                
                xh=xnm(k);
                while ind1<lx && xnm(ind1)<xh-dist
                    ind1=ind1+1;
                end
%                 ind2=ind1;
                while ind2<=lx && xnm(ind2)<xh+dist;
                    ind2=ind2+1;
                end
                
                indrange=ind1:ind2-1;
                indregion=(xnm(indrange)-xnm(k)).^2+(ynm(indrange)-ynm(k)).^2<dist^2;
%                 indin(sortind(k))
                vfilt=locfilter(value(indrange(indregion)));
                outval(indin(sortind(k)))=vfilt;
                
                if toc(timerx)>1
                    timerx=tic;
                    obj.status([num2str(100*k/lx,3) '%']);
                    drawnow
                end
                
            end
            
            


                if isfield(obj.locData.loc,p.resultfield)
                    obj.locData.loc.(p.resultfield)(indin)=outval(indin);
                else
                 obj.locData.setloc(p.resultfield,outval);
                end
            end
                 obj.locData.filter(p.resultfield)
                 obj.locData.regroup;
                 obj.setPar('locFields',fieldnames(obj.locData.loc))

             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function locfield_callback(obj,a,b)
            p=obj.getGuiParameters;
            field=p.locfield.selection;
            obj.guihandles.resultfield.String=[field '_' p.filter.selection];
            
        end
    end
end

function pard=guidef(obj)
pard.t1.object=struct('String','Only localizations displayed in layer 1 are used','Style','text');
pard.t1.position=[1,1];
pard.t1.Width=4;

pard.resultfield.object=struct('String','','Style','edit');
pard.resultfield.position=[3,1];
pard.resultfield.Width=0.8;


pard.t2.object = struct('String','=','Style','text');
pard.t2.position=[3,1.8];
pard.t2.Width=0.2;

pard.locfield.object=struct('String','','Style','popupmenu','Callback',{{@obj.locfield_callback}});
pard.locfield.position=[3,2];
pard.locfield.Width=1;

pard.filter.object=struct('String',{{'median','mean','max','min','std','5-95% percentile'}},'Style','popupmenu','Callback',{{@obj.locfield_callback}});
pard.filter.position=[3,3];
pard.filter.Width=1;

pard.t3.object = struct('String','Spatial scale (nm)','Style','text');
pard.t3.position=[4,1];
pard.t3.Width=1;

pard.spatialscale.object = struct('String','100','Style','edit');
pard.spatialscale.position=[4,2];
pard.spatialscale.Width=1;
% pard.t3.object=struct('String','Equation. Use fieldnames (e.g. xnm, phot) as variables','Style','text');
% pard.t3.position=[2,2];
% pard.t3.Width=3;
% 
% pard.equation.object = struct('String','(locprecnm<25 & PSFxnm>100) | numberInGroup>1','Style','edit');
% pard.equation.position=[3,2];
% pard.equation.Width=3;


pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];
pard.dataselect.object.TooltipString='choose localization file data set';


pard.dataselectall.object=struct('Style','checkbox','String','all','Value',1);
pard.dataselectall.position=[2,2];
pard.dataselectall.object.TooltipString='apply on all';


pard.syncParameters={{'filelist_short','dataselect',{'String'}},{'locFields','locfield',{'String'}}};

pard.plugininfo.type='ProcessorPlugin';
end


function out=perc9_95(in)
q=myquantilefast(in,[0.05 0.95]);
out=q(2)-q(1);
end