classdef VersatileRenderer<interfaces.DialogProcessor
%     Render any field vs any other using any 3rd field for coloring
    methods
        function obj=VersatileRenderer(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_layerson','numberOfLayers'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            %additional field for sx, sy. checkbox: set
            out=[];
            ax=initaxis(p.resultstabgroup,'image');
            
            
            lp=obj.getLayerParameters;
           
            for k=length(lp):-1:1
                plotfields{k}=lp{k}.renderfield.selection;
            end
           
            [~,inroi]=obj.locData.getloc('xnm','position','roi');
             lochere=obj.locData.copy(unique({plotfields{:} p.assignfield1.selection p.assignfield2.selection 'locprecnm' 'znm' 'xnm' 'ynm' 'phot' 'numberInGroup'}),inroi);
%             lochere.removelocs(~inroi);
%             lochere.regroup;
            v1=lochere.loc.(p.assignfield1.selection)/p.pixelsize1;
            v2=lochere.loc.(p.assignfield2.selection)/p.pixelsize2;
            lochere.loc.x=v1;
            lochere.loc.y=v2;

            v1=lochere.grouploc.(p.assignfield1.selection)/p.pixelsize1;
            v2=lochere.grouploc.(p.assignfield2.selection)/p.pixelsize2;
            lochere.grouploc.x=v1;
            lochere.grouploc.y=v2;
            
            pall=obj.getLayerParameters;
%             if p.setsigma1
                    lochere.loc.sx=lochere.loc.x*0+p.sigma1;
                    lochere.grouploc.sx=lochere.grouploc.x*0+p.sigma1;
%             else
%                     lochere.loc.sx=lochere.loc.(p.sigmafield1.selection)/p.pixelsize1;
%                     lochere.grouploc.sx=lochere.grouploc.(p.sigmafield1.selection)/p.pixelsize1;
%           
%             end
%             if p.setsigma2
                    lochere.loc.sy=lochere.loc.x*0+p.sigma2;  
                    lochere.grouploc.sy=lochere.grouploc.x*0+p.sigma2;
                    
%             else   
%                     lochere.loc.sy=lochere.loc.(p.sigmafield2.selection)/p.pixelsize2;
%                     lochere.grouploc.sy=lochere.grouploc.(p.sigmafield2.selection)/p.pixelsize2;
%             end
%             phere.gaussfac=1;
            phere.sr_pixrec=1;
            phere.rangex=[p.min1 p.max1]/p.pixelsize1;
            phere.rangey=[p.min2 p.max2]/p.pixelsize2;
%             phere.mingausspix=p.filtersize;
%             phere.gaussfac=.01;
%             phere.mingausspix=0;
            phere.mingaussnm=0;
            
            phere.sr_axes=[];
            phere.addscalebar=false;
            for k=1:length(pall)
                pall{k}=copyfields(pall{k},p);
                pall{k}=copyfields(pall{k},phere);
            end
            
            img=TotalRender(lochere,pall,{'xnm','ynm'});
            imagesc(ax, [p.min1 p.max1],[p.min2 p.max2],img.image)
            xlabel(ax,p.assignfield1.selection)
            ylabel(ax,p.assignfield2.selection)
            axis(ax,'xy')
            
            linemode=p.plotline.selection;
            if ~strcmp(linemode,'none')
                rangex=p.min1:p.pixelsize1:p.max1;
                layers=find(obj.getPar('sr_layerson'));
                 hold(ax,'on')
         
                for k=1:length(layers)
                    lp{k}.lut.selection
                    lut=mymakelut(lp{k}.lut.selection);
                    if lp{k}.lutinv
                        lut=1-lut;
                    end                    
%                     lutm=mean(lut,1);
                    lutm=lut(end/2,:);
                    lutm=lutm/max(lutm);

                    loch=obj.locData.getloc({p.assignfield1.selection, p.assignfield2.selection },'position','roi','layer',layers(k));
                    v1=loch.(p.assignfield1.selection);
                    v2=loch.(p.assignfield2.selection);
                    linedat=bindata(v1,v2,rangex,linemode);
                    plot(ax,rangex,linedat,'k','LineWidth',3)
                    plot(ax,rangex,linedat,'Color',lutm)
                end
                hold(ax,'off')
            end
            
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end
function assignfield_callback(object,a,obj,axis)
    field=object.String{object.Value};
    v=obj.locData.getloc(field,'layer',find(obj.getPar('sr_layerson')),'position','roi').(field);
    
    v=v(~isnan(v) & ~isinf(v));
    q=myquantilefast(v,[0.001,0.5,0.999]);

    dx=(q(3)-q(1))/100;
    if dx==0
        dx=1;
    end
    dx=10^floor(log10((dx)));
    obj.guihandles.(['pixelsize' num2str(axis)]).String=num2str(dx);
    minv=floor(q(1)/dx)*dx-5*abs(dx);
    maxv=ceil(q(3)/dx)*dx+5*abs(dx);
    obj.guihandles.(['min' num2str(axis)]).String=num2str(minv);
    obj.guihandles.(['max' num2str(axis)]).String=num2str(maxv);
end
function pard=guidef(obj)
pard.t1.object=struct('String','horizontal axis','Style','text');
pard.t1.position=[1,2];
pard.t2.object=struct('String','vertical axis','Style','text');
pard.t2.position=[1,3];

pard.assignfield1.object=struct('Style','popupmenu','String',{{'n1,n2'}},'Callback',{{@assignfield_callback,obj,1}});
pard.assignfield1.position=[2,2];
pard.assignfield2.object=struct('Style','popupmenu','String',{{'n1,n2'}},'Callback',{{@assignfield_callback,obj,2}});
pard.assignfield2.position=[2,3];
pard.assignfield1.object.TooltipString='choose which field to render';
pard.assignfield2.object.TooltipString=pard.assignfield1.object.TooltipString;

pard.t3.object=struct('String','pixelsize','Style','text');
pard.t3.position=[3,1];
pard.t4.object=struct('String','min','Style','text');
pard.t4.position=[4,1];
pard.t5.object=struct('String','max','Style','text');
pard.t5.position=[5,1];
% pard.t7.object=struct('String','sigma','Style','text');
% pard.t7.position=[6,1];


% pard.sigmafield1.object=struct('Style','popupmenu','String',{{'n1,n2'}});
% pard.sigmafield1.position=[6,2];
% pard.sigmafield2.object=struct('Style','popupmenu','String',{{'n1,n2'}});
% pard.sigmafield2.position=[6,3];

pard.setsigma.object=struct('String','sigma','Style','text');
pard.setsigma.position=[6,1];
% pard.setsigma1.object=struct('String','','Style','checkbox','Value',1);
% pard.setsigma1.position=[7,2];
% pard.setsigma1.Width=.2;

% pard.setsigma2.object=struct('String','','Style','checkbox','Value',1);
% pard.setsigma2.position=[7,3];
% pard.setsigma2.Width=.2;

pard.sigma1.object=struct('String','1','Style','edit');
pard.sigma1.position=[6,2];
pard.sigma1.Width=1;
pard.sigma2.object=struct('String','1','Style','edit');
pard.sigma2.position=[6,3];
pard.sigma2.Width=1;

pard.pixelsize1.object=struct('String','','Style','edit');
pard.pixelsize1.position=[3,2];
pard.pixelsize2.object=struct('String','','Style','edit');
pard.pixelsize2.position=[3,3];

pard.min1.object=struct('String','','Style','edit');
pard.min1.position=[4,2];
pard.min2.object=struct('String','','Style','edit');
pard.min2.position=[4,3];

pard.max1.object=struct('String','','Style','edit');
pard.max1.position=[5,2];
pard.max2.object=struct('String','','Style','edit');
pard.max2.position=[5,3];

pard.plotlinet.object=struct('String','Plot line:','Style','text');
pard.plotlinet.position=[7,1];
pard.plotline.object=struct('String',{{'none','mean','median','std','robustmean'}},'Style','popupmenu');
pard.plotline.position=[7,2];

% pard.t6.object=struct('String','sigma rec (pix)','Style','text');
% pard.t6.position=[3,4];
% pard.filtersize.object=struct('String','1','Style','edit');
% pard.filtersize.position=[2,4];

pard.syncParameters={{'locFields','assignfield1',{'String'}},{'locFields','assignfield2',{'String'}}};%,{'locFields','sigmafield2',{'String'}},{'locFields','sigmafield1',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Render any field vs any other using any 3rd field for coloring';
end

