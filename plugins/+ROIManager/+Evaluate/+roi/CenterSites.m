classdef CenterSites<interfaces.SEEvaluationProcessor
%     Calculates the center coordinates of structres in a ROI based on
%     median, mean, circular fit or cross-correlation with a mask and
%     shifts thte ROI to center the structure.
    properties
        
    end
    methods
        function obj=CenterSites(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function makeGui(obj)
            makeGui@interfaces.SEEvaluationProcessor(obj);
%             obj.guihandles.saveimagesb.Callback={@saveimagesb_callback,obj};
        end
        function out=run(obj,p)
            try
            out=runintern(obj,p);
            catch err
                err
                out=[];
            end
         
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)
p(1).value=[1,2];p(1).on={'cz'}; p(1).off={'radiust','radius'}; %mean, median
p(2).value=3;p(2).on={'radiust','radius','blurrt','blurr'}; p(2).off={'cz'};
p(3).value=4;p(3).on={'radiust','radius'}; p(3).off={'blurrt','blurr','cz'};

pard.centermode.object=struct('Style','popupmenu','String',{{'median','mean','mask','fitring'}},'Callback',{{@obj.switchvisible,p}});
pard.centermode.position=[1,1];
pard.centermode.Width=2;

pard.radiust.object=struct('Style','text','String','R (nm)','Visible','off');
pard.radiust.position=[2,1];
pard.radiust.Width=1;
pard.radius.object=struct('Style','edit','String','50','Visible','off');
pard.radius.position=[2,2];
pard.radius.Width=1;

pard.blurrt.object=struct('Style','text','String','sigma (nm)','Visible','off');
pard.blurrt.position=[3,1];
pard.blurrt.Width=1;
pard.blurr.object=struct('Style','edit','String','20','Visible','off');
pard.blurr.position=[3,2];
pard.blurr.Width=1;


pard.iterationst.object=struct('Style','text','String','iterations:');
pard.iterationst.position=[1,3];
pard.iterationst.Width=1;

pard.iterations.object=struct('Style','edit','String','1');
pard.iterations.position=[1,4];
pard.iterations.Width=.5;

pard.cxy.object=struct('Style','checkbox','String','correct x,y');
pard.cxy.position=[2,3];
pard.cxy.Width=2;

pard.cz.object=struct('Style','checkbox','String','correct z');
pard.cz.position=[3,3];
pard.cz.Width=2;

pard.plugininfo.type='ROI_Evaluate';
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','layer1_','layer2_','se_sitepixelsize'};
pard.plugininfo.description='Calculates the center coordinates of structres in a ROI based on median, mean, circular fit or cross-correlation with a mask and shifts thte ROI to center the structure.';
end


function out=runintern(obj,p)
for k=1:p.iterations
%     obj.site.pos(1)
    locs=obj.getLocs({'xnm','ynm','znm','locprecznm','locprecnm'},'layer',find(obj.getPar('sr_layerson')),'size',p.se_siteroi/2);
    if obj.display
    ax=obj.setoutput('images');
    end
    switch p.centermode.selection
        case 'median'
            dx=median(locs.xnm)-obj.site.pos(1);
            dy=median(locs.ynm)-obj.site.pos(2);
            if ~isempty(locs.znm)
            dz=median(locs.znm)-obj.site.pos(3);
            else
                dz=0;
            end
            obj.site.pos=obj.site.pos+[dx*p.cxy dy*p.cxy dz*p.cz];
        case 'mean'
            dx=mean(locs.xnm)-obj.site.pos(1);
            dy=mean(locs.ynm)-obj.site.pos(2);
            if ~isempty(locs.znm)
            dz=mean(locs.znm)-obj.site.pos(3);
            else
                dz=0;
            end
            obj.site.pos=obj.site.pos+[dx*p.cxy dy*p.cxy dz*p.cz];
        case 'mask'
            img=obj.site.image.image;
            
            sim=size(img);
            pm=floor((sim(1)-1)/2);
            ns=(-pm:pm)*p.se_sitepixelsize;
            ns2=(-pm:0.25:pm)*p.se_sitepixelsize;
            
%             ns=round(-p.se_sitefov/2:p.se_sitepixelsize:p.se_sitefov/2);
            [Xs,Ys]=meshgrid(ns,ns);
            mask=Xs.^2+Ys.^2<(p.se_siteroi/2)^2;

            pixrec=p.se_sitepixelsize;
            rim=p.blurr;
            diameterNPC=p.radius*2;
            if isempty(diameterNPC)
                diameterNPC=p.se_siteroi/2;
                disp('please specify radius')
            end
            n=-1*diameterNPC:pixrec:1*diameterNPC;
            [X,Y]=meshgrid(n,n);

            h=exp(-abs(X.^2+Y.^2-(diameterNPC/2)^2)/4/rim^2);
            h=h/sum(h(:));
            imsq=sqrt(sum(img,3));
            imsq(~mask)=0;
            imf=filter2(h,imsq);
            imff=imresize(imf,4,'Method','bicubic');
            [~,mind]=max(imff(:));
            [my,mx]=ind2sub(size(imff),mind);
            dx=ns2(mx);
            dy=ns2(my);   
            dz=0;
            
            %plot 
            se=obj.locData.SE;
            obj.site.pos=obj.site.pos+[dx*p.cxy dy*p.cxy dz*p.cz];
            obj.site.image=[];
            se.plotsite(obj.site);
            if obj.display
                imagesc(ax,[ imf/max(imf(:)),imsq/max(imsq(:))])
                title(ax,[dx dy])
                ax2=obj.setoutput('mask');
                imagesc(ax2,h);
                xlabel(ax2,'x (nm)')
                ylabel(ax2,'y (nm)')
                xlabel(ax,'x (nm)')
                ylabel(ax,'y (nm)')
            end
        case 'fitring'
            posold=obj.site.pos;
            [x0,y0,R]=fitposring(locs.xnm,locs.ynm,p.radius);
            if p.cxy
                obj.site.pos(1:2)=[x0,y0];
            end
            if obj.display
                plot(ax,locs.xnm,locs.ynm,'.')
                hold(ax,'on')
                circle(x0,y0,R,'Parent',ax)
                hold(ax,'off')
                dx=x0-posold(1);
                dy=y0-posold(2);
                dz=0;
                title(ax,[dx dy R])
                xlabel(ax,'x (nm)')
                ylabel(ax,'y (nm)')
                axis(ax,'equal')
            end
    end 
end

out.dx=dx;
out.dy=dy;
out.dz=dz;
end

