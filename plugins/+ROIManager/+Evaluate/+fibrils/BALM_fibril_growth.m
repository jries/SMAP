classdef BALM_fibril_growth<interfaces.SEEvaluationProcessor
    properties
        poly
        hpoly
    end
    methods
        function obj=BALM_fibril_growth(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj,p)
            out=[];
            lw=p.width/2;
            locs=obj.getLocs({'locprecnm','PSFxnm','xnm','ynm','frame'},'layer',1,'size',p.se_siteroi);  
            pol=obj.site.annotation.rotationpos.pos;
            
            dpol=pol(2,:)-pol(1,:);
            alpha=-atan2(dpol(1),dpol(2));

            len=sqrt(sum(dpol.^2))*1000/2;
            if len==0
                out=[];
                return
            end
            midp=mean(pol,1)*1000;
            [xr,yr]=rotcoord(locs.xnm-midp(1),locs.ynm-midp(2),alpha);

            indb=abs(xr)>lw|abs(yr)>len;
            indg=~indb;
            ynmline=xr(indg);xnmline=yr(indg);
            frame=locs.frame(indg);
            
            %segment filament, only keep locs in filament
            pixrec=10; %nm
            xrange=-len-pixrec:pixrec:len+pixrec;
            yrange=-lw-pixrec:pixrec:lw+pixrec;
            him=histcounts2(xnmline, ynmline,xrange,yrange);
            sigma=0.5;
            h=fspecial('gaussian',12,sigma);
            max1=mean(him(:));
            
            himf=fibermetric(him,round(p.fiberw/pixrec),'StructureSensitivity',max1/3);
            himf=filter2(h,himf);
             max1=mean(himf(:));
%             max1=1/sigma^2/pi/2;
%             max1=4*max(h(:));
            
            imbw=himf>max1;
%             figure(89);subplot(2,2,1);imagesc(himf);
            imbw=bwareaopen(imbw,round(sum(imbw(:))/4));
%              subplot(2,2,2);imagesc(imbw);
            
             seb=strel('disk',3,4);
%              ses=strel('disk',1,4);
%              imbw=imerode(imbw,ses);
%              imbw=bwareaopen(imbw,round(sum(imbw(:))/4));
            imbw=imdilate(imbw,seb);
            imbw=imerode(imbw,seb);
            imbw=imdilate(imbw,seb);
             h=obj.setoutput('images');
             himp=him/max(him(:));
             imagesc(h,horzcat(himp,himf,double(imbw)))
%             
%            imbw=imdilate(imbw,ses);
%             imbw=imerode(imbw,ses);
%             imbw=imopen(imbw,seb);
%              imbw=imerode(imbw,ses);
%               imbw=imerode(imbw,ses);

           
           
         
%             imbw=imsegfmm(himf,imbw,0.01);
%             imbw=imdilate(imbw,ones(3));
%             bwac=activecontour(himf,true(size(himf)));
%             subplot(2,2,3);imagesc(imbw)
            xr2=round((xnmline+len)/pixrec)+1;
            yr2=round((ynmline+lw)/pixrec)+1;
            indlin=sub2ind(size(imbw),xr2,yr2);
            indgood=imbw(indlin);
            out.mask=imbw;
            
            %
%             h2=fspecial('gaussian',12,1.5);
            
 
            xn=-len:p.dx:len;
            fr=1:p.df:max(frame(indgood));
            him=histcounts2(frame(indgood),xnmline(indgood),fr,xn);
%             hxf=filter2(h2,him);
%             mf=max(h2(:));
%             hxfb=hxf>mf*1.1;
            
            out.kimograph=him;
            out.fr=fr;
            out.xn=xn;
            co=quantile(him(:),0.999);
            him(him>co)=co;
            h=obj.setoutput('kimograph');
            imagesc(h,xn,fr,(him))
            xlabel(h,'xnm')
            ylabel(h,'frame')
            
            uicontrol(h.Parent,'Position',[10,10,40,20],'String','Add Line','Callback',{@obj.addline,h});
            
            if length(obj.poly)<obj.site.ID || isempty(obj.poly{obj.site.ID})
                valh=obj.site.evaluation.(obj.modulename);
                if isfield(valh,'poly')
                   %obj.hpoly=impoly(h,valh.poly,'Closed',false,'PositionConstraintFcn',@obj.polyconstrain);
                   obj.poly{obj.site.ID}=valh.poly;
                end
            end
            if length(obj.poly)>=obj.site.ID &&~isempty(obj.poly{obj.site.ID})
                obj.hpoly=impoly(h,obj.poly{obj.site.ID},'Closed',false,'PositionConstraintFcn',@obj.polyconstrain);
            
            
                 setConstrainedPosition(obj.hpoly,obj.hpoly.getPosition);
                hi=addNewPositionCallback(obj.hpoly,@obj.polycallback);
    %             obj.hpoly.wait;
                obj.poly{obj.site.ID}=obj.hpoly.getPosition;
                out.poly=obj.hpoly.getPosition;
                obj.polycallback(obj.hpoly.getPosition);
            end
           

        end
        
        function addline(obj,a,b,h)
           % if length(obj.poly)<obj.site.ID || isempty(obj.poly{obj.site.ID}) 
                obj.hpoly=impoly(h,'Closed',false,'PositionConstraintFcn',@obj.polyconstrain);
           % else %add vortex to end
           %     pos=obj.hpoly.getPosition;
           %     pos(end+1,:)=pos(end,:);
           %     obj.hpoly.setPosition(vertcat(pos(1,:),pos));
                 
          %  end
             setConstrainedPosition(obj.hpoly,obj.hpoly.getPosition);
            hi=addNewPositionCallback(obj.hpoly,@obj.polycallback);
            obj.polycallback(obj.hpoly.getPosition);
        end
        
        function outp=polyconstrain(obj,inp)
            
            outp=vertcat(inp(1,:),max(inp(2:end,:),inp(1:end-1,:)));
            
        end
        function polycallback(obj,inp)
%             setConstrainedPosition(obj.hpoly,inp);
            obj.poly{obj.site.ID}=inp;
            obj.site.evaluation.(obj.name).poly=inp;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end

end

function pard=guidef
pard.dxt.object=struct('Style','text','String','dx (nm), dt (frames)');
pard.dxt.position=[1,1];
pard.dxt.Width=2;
pard.dx.object=struct('Style','edit','String','10');
pard.dx.position=[1,3];
pard.df.object=struct('Style','edit','String','1000');
pard.df.position=[1,4];

pard.widtht.object=struct('Style','text','String','width (nm)');
pard.widtht.position=[2,1];
pard.widtht.Width=2;
pard.width.object=struct('Style','edit','String','400');
pard.width.position=[2,3];
pard.width.Width=2;

pard.fiberwt.object=struct('Style','text','String','fiber width (nm)');
pard.fiberwt.position=[3,1];
pard.fiberwt.Width=2;
pard.fiberw.object=struct('Style','edit','String','150');
pard.fiberw.position=[3,3];
pard.fiberw.Width=2;

% pard.dxt.Width=3;
pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi'};
pard.plugininfo.type='ROI_Evaluate';
end