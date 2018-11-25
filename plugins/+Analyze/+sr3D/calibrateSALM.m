classdef calibrateSALM<interfaces.DialogProcessor
    properties
        
    end
    methods
        function obj=calibrateSALM(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;

        end
        function out=run(obj,p)
            locs=obj.locData.getloc({'filenumber','xnm','ynm','frame','znm','phot','numberInGroup','groupindex','PSFxpixGauss','PSFypixGauss'},...
                'layer', find(obj.getPar('sr_layeron')),'grouping','ungrouped','position','roi',...
                'removeFilter',{'filenumber'});
           dz=50;
           maxzabs=500;
           minframes=maxzabs/dz * 0.5;

           indgood=abs(locs.znm)<maxzabs; %only close to fitted z, here other filters could apply
           ngroup=1:max(locs.groupindex);
           hcb= histcounts(locs.groupindex(indgood),ngroup);
           beadgroup=ngroup(hcb>=minframes);
            
%            symmetricPSF=isempty(locs.PSFypixGauss);
           symmetricPSF=true;
           
           for k=length(beadgroup):-1:1
               indh=(locs.groupindex==beadgroup(k)) & indgood;
               PSFx=locs.PSFxpixGauss(indh);
               frames=locs.frame(indh);
               if symmetricPSF
                   f0(k)=getf02D(frames,PSFx);
               else
                   f0(k)=getf03D(frames,PSFx,PSFy);
               end
               filenumber(k)=locs.filenumber(find(indh,1,'first'));
           end
           figure(99);histogram(f0,floor(min(f0)):1:max(f0))
              

        end
        
        function initGui(obj)           
        end
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function f0=getf02D(frames,PSFx)
wind=3;
[~,ind]=min(PSFx);
r=max(1,ind-wind):min(length(frames),ind+wind);
fitp=fit(frames(r),PSFx(r),'poly2');
f0=-fitp.p2/2/fitp.p1;
end

function [Zint]=getZinterp(beads,Zintold,p,zaxis)
if nargin<4
    zaxis='zglass';
end
 %make big array for interpolation
 zglassall=[];z0relativeall=[];zfitall=[];idall=[];zplot=[];dzerrall=[];z0glassall=[];
 
            for k=1:length(beads)
                
                zglassall=double(vertcat(zglassall,beads(k).loc.zglass));
                z0glassall=double(vertcat(z0glassall,beads(k).loc.z0glass));
                z0relativeall=double(vertcat(z0relativeall,beads(k).loc.z0relative));
%                 z0all=double(vertcat(z0all,beads(k).f0*p.dz*ones(length(beads(k).loc.zglass),1)));

                zfitall=double(vertcat(zfitall,beads(k).loc.znm));
%                 fileall=double(vertcat(fileall,double(beads(k).filenumber)*ones(length(beads(k).loc.zglass),1)));
                idall=double(vertcat(idall,k*ones(length(beads(k).loc.zglass),1)));
                if 0%beads(k).f0glass<50
                    zplot=double(vertcat(zplot,beads(k).loc.dzcorr*0));
                else
                    zplot=double(vertcat(zplot,beads(k).loc.dzcorr));
                end
                
%                 if ~isempty(dzerr)
%                     dzerrall=double(vertcat(dzerrall,dzerr{k}));
%                 else
%                     dzerrall=double(vertcat(dzerrall,0*beads(k).loc.znm));
%                 end
            end
            
            if strcmp(zaxis,'z0glass')
                zzax=z0glassall;
            else
                zzax=zglassall;
            end
            
            if 0%p.setzero
                zplot(z0glassall<100)=0;
            end
%             zglassall=z0glassall;
%             zglassall=zglassall-
            qzfit=myquantile(zfitall,[0.05,0.95]);
            qzfit(1)=qzfit(1)+p.dz;qzfit(2)=qzfit(2)-p.dz;
%             qzfit(1)=max(qzfit(1),p.zrangeuse(1));
%             qzfit(2)=min(qzfit(2),p.zrangeuse(2));
            %for now don't follow up, later: remove single wrong
            %localizations.
            
            %determine range/outliers
%             inz=z0relativeall>p.zrangeuse(1)&z0relativeall<p.zrangeuse(2);
            inz=abs(z0relativeall)<800;
%             inz=inz&z0glassall>-50;
            inz=inz&(zfitall)<qzfit(2)&(zfitall)>qzfit(1);
%             figure(78);scatter3(zglassall,zfitall,z0all)
%             zplot=z0relativeall-zfitall;
            inz=inz&abs(zplot)<800;
            
            if ~isempty(Zintold)  
                dz=Zintold.interp(zzax,zfitall)-zplot;

                inz=inz&abs(dz)<p.cutoffrefine;
                h=histcounts(idall(inz),(1:max(idall)+1))';
%                 minpoints=(p.zrangeuse(2)-p.zrangeuse(1))/p.dz/2; %ONLY relevant use of zrangeuse. remove!
                minpoints=800/p.dz;
                innump=h(idall)>minpoints;
                inz=inz&innump;
            end
            if p.setzero
                    
                    zfitx=(qzfit(1):p.dz:qzfit(2))';
                    zfitallh=vertcat(zfitall(inz),zfitx,zfitx,zfitx);
                    zploth=vertcat(zplot(inz),0*zfitx,0*zfitx,0*zfitx);
%                     zploth(z0glassall(inz)<50)=0;
                    zzaxh=vertcat(zzax(inz),min(zzax)+0*zfitx,min(zfitx)+0*zfitx,mean([min(zfitx),min(zzax)])+0*zfitx);
                
            else
                zfitallh=zfitall(inz);
                zploth=zplot(inz);
                zzaxh=zzax(inz);

            end            
            
            xrange=min(zzax(inz)):100:max(zzax(inz));
            
            yrange=[qzfit(1):10: qzfit(2)];
            [X,Y]=meshgrid(xrange,yrange);  
            %interpolation

            idallh=idall(inz);
            idallh(length(zzaxh))=0;
            Z=RegularizeData3D(zzaxh,zfitallh,zploth,xrange,yrange,'smoothness',p.smoothing);
             Zint.interp=griddedInterpolant(X',Y',Z');
             Zint.zaxis=zaxis;
           
            
%              Z=RegularizeData3D(zglassall(inz),zfitall(inz),zplot(inz),xrange,yrange,'smoothness',[.1,.001]);
%             
%             Zint=griddedInterpolant(X',Y',Z');
            
            if ~isempty(p.axhere)
                scatter3(p.axhere,zzax(inz),zfitall(inz),zplot(inz),[],idall(inz))

%                 scatter3(p.axhere,zzaxh,zfitallh,zploth,[],idallh)
                xlabel(p.axhere,'zglass');ylabel(p.axhere,'zfit'); zlabel(p.axhere,'z0r');
                colormap(p.axhere,'lines')
                hold(p.axhere,'on')
                mesh(p.axhere,X,Y,Zint.interp(X,Y),'FaceAlpha',0.2)
            end
end
function  f0glass=getf0glass(beads,p)
if isempty(p.axhere)
    f=figure;ax=gca;
else
    ax=p.axhere;
end
for k=1: max([beads(:).filenumber]) 
    indf=[beads(:).filenumber]==k;
      f0=[beads(indf).f0];
      phot=[beads(indf).phot];
      dzh=50/p.dz;
      induse=f0<dzh*60;
%       induse=induse&phot>median(phot);
%       mean(phot)
    
    f0=f0(induse); %only look in the first 3 um
    
    range=0:dzh:1000;
    h=histogram(ax,f0,range);


    [mh]=max(h.Values);
    ind=find(h.Values>mh*.4,1,'first');
    f0h=range(ind);
    ind=find(f0>f0h-2*dzh&f0<f0h+2*dzh);
    f0glass(k)=mean(f0(ind));
    hold(ax, 'on');
    
end
if isempty(p.axhere)
    close(f)
else
    plot(ax,f0glass,ones(size(f0glass)),'k*')
end
end
function [err1,dzerr]=geterrors(beads,Zint,p)
%   figure(99)
%             hold off
            xrange=Zint.interp.GridVectors{1};
            yrange=Zint.interp.GridVectors{2};

  for k=1:length(beads)
        zh=double(beads(k).loc.znm);
        zglass=beads(k).loc.zglass;
        z0f=beads(k).loc.dzcorr;
        inz=abs(zh<300) & abs(z0f)<300 & (zh)<yrange(end) & (zh)>yrange(1);
%                 inz= (zh)<qzfit(2) & (zh)>qzfit(1);
        dz=Zint.interp(zglass(inz),zh(inz))-z0f(inz);
% 
%         plot(zh(inz),dz)
%         hold on
  if p.setzero&&beads(k).f0glass*p.dz<2*p.dz %on glass
      factor=0.2;
  else
      factor=1;
  end        
        err1(k)=sqrt(mean(dz.^2))*factor;
        err2(k)=mean(abs(dz))*factor;
        err3(k)=std(dz)*factor;
%         dzerr{k}=ones(size(beads(k).loc.znm))+NaN;
         dzerr{k}=Zint.interp(zglass,zh)-z0f;
   end
end
function getcoords(a,b,obj)
locsall=obj.locData.getloc({'frame','xnm','ynm','PSFxnm','PSFynm','filenumber','phot'},'position','all','layer',1,'removeFilter','filenumber');
p=obj.getAllParameters;
x=locsall.xnm/p.cam_pixelsize_nm;
y=locsall.ynm/p.cam_pixelsize_nm;
img=myhist2(x,y,1,1,[0 512.5],[0 512.5]);
f1=figure;
imagesc(img');
h=imrect;
pos=wait(h);% [x y wx wy]
delete(f1);
answ=inputdlg({'number of rows','number of columns'},'set tiles', 1,{'2','1'});
if isempty(answ)
    return;
end
nx=str2double(answ{2});
ny=str2double(answ{1});

% pos=pos([2 1 4 3]);
p.Xmin=round(pos(1)); p.Ymin=round(pos(2)); p.Xmax=round(pos(1)+pos(3)); p.Ymax=round(pos(2)+pos(4));
p.Xd=floor(pos(3)/nx);
p.Yd=floor(pos(4)/ny);
obj.setGuiParameters(p)

end
function pard=guidef(obj)


pard.t1.object=struct('String','value 1','Style','text');
pard.t1.position=[4,1];
pard.t2.object=struct('String','value 2','Style','text');
pard.t2.position=[4,2];

pard.assignfield1.object=struct('Style','popupmenu','String','');
pard.assignfield1.position=[5,1];
pard.assignfield2.object=struct('Style','popupmenu','String','');
pard.assignfield2.position=[5,2];


 pard.syncParameters={{'locFields','assignfield1',{'String'}},{'locFields','assignfield2',{'String'}}};
            
%             obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
%             obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
pard.plugininfo.type='ProcessorPlugin';


end
