classdef calibrate3Daberrations<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=calibrate3Daberrations(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
%              obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-1;obj.guiPar.Vrim=obj.guiPar.Vrim-20;
        end
        function out=run(obj,p)
            p.EMon=obj.locData.files.file(1).info.EMon;
            fp=obj.locData.history{1}.children.fitparamters;
            if isfield(fp,'fitterGUI')
                fp=fp.fitterGUI.children.MLE_GPU_Yiming;
            elseif isfield(fp,'MLE_GPU_Yiming')
                fp=fp.MLE_GPU_Yiming;
            end
            if fp.userefractive_index_mismatch
                p.RIM=fp.refractive_index_mismatch;
            else 
                p.RIM=1;
            end
            p.dz=p.dz*p.RIM;
            out=[];
            %             * Make 3D fitting model
            % * Fit bead stacks (gel, cells) with this model 
            %     * (refractive index mismatch =1) or rescale later
            % * extract bead positions (grouping)
            disp('get beads')
            if p.beadsource.Value==1 %ROImanager
                beads=sites2beads(obj,p);
            else % segment new
                beads=segmentb(obj,p);
            end
            
            %framerange
            frange=[min(obj.locData.loc.frame) max(obj.locData.loc.frame)];
            
            % get f0 for beads
            for k=length(beads):-1:1
                   beads(k).loc.znm=beads(k).loc.znm;
                [beads(k).f0]=getf0Z(beads(k).loc,p);
                beads(k).phot=beads(k).loc.phot(min(length(beads(k).loc.phot),max(1,round(beads(k).f0))));
            end
            f0all=([beads(:).f0]);
            badind=f0all<frange(1)|isnan(f0all); %|f0all>frange(2)
            beads(badind)=[];
            
            % * determine fglass, glass
            if p.setglass
                if isempty(p.glassframe)
                    beadst=beads;
                    for k=1:length(beadst)
                        beadst(k).filenumber=1;
                    end
                    p.axhere=obj.initaxis('z0positions');
                    p.glassframe=getf0glass(beadst,p);
                end
                f0glass=p.glassframe*ones(1,max([beads(:).filenumber]));
            else
                p.axhere=obj.initaxis('z0positions');
                f0glass=getf0glass(beads,p);
            end
            p.f0glass=f0glass;
            
            %calculate relevant other coordinates
            axh=obj.initaxis('z vs frame');
            hold off
            for k=1:length(beads)
                beads(k).loc.zglass=(beads(k).loc.frame-f0glass(beads(k).filenumber))*p.dz;
                beads(k).loc.z0relative=-(beads(k).f0-beads(k).loc.frame)*p.dz;
                beads(k).loc.dzcorr=beads(k).loc.z0relative-beads(k).loc.znm;
                beads(k).f0glass=beads(k).f0-f0glass(beads(k).filenumber);
                beads(k).loc.z0glass=beads(k).f0glass*p.dz+0*beads(k).loc.zglass;
                 indplot=abs(beads(k).loc.z0relative)<800;
%                 indplot=(beads(k).loc.z0relative)>p.zrangeuse(1)&(beads(k).loc.z0relative)<p.zrangeuse(2);
%                 indplot=true(size(beads(k).f0glass));
%                 indplot=(beads(k).loc.z0relative)>p.zrangeuse(1)&(beads(k).loc.z0relative)<p.zrangeuse(2);
%                 plot(axh,beads(k).loc.frame,beads(k).loc.znm,'.')
                plot(axh,beads(k).loc.zglass(indplot),beads(k).loc.znm(indplot),'.')
%                 beads(k).stddz=std(diff(beads(k).loc.znm(indplot)));
                hold on
            end
            
            p.axhere=[];
%             p.setzero=true;
            %get interpolation
            phere=p;
            phere.smoothing=[0.05 0.002];
           [ZcorrInterp]=getZinterp(beads,[],phere,'zglass');
           phere=p;
           phere.cutoffrefine=500;
           [ZcorrInterp]=getZinterp(beads,ZcorrInterp,phere);
%            p.cutoffrefine=150;
            %calculate errors
           [err1,dzerr]=geterrors(beads,ZcorrInterp,p);
           err0=err1;
           goodind=find(true(length(beads),1));
           beads2=beads;

           cutofffactor=4;
           while  1% length(beads2)>length(beads)/2
                cutoff=cutofffactor*nanmean(err1);
                badind=(err1>cutoff|isnan(err1));
                if sum(badind)==0
                    break
                end
                goodind=goodind(~badind);
                beads2=beads(goodind);
%                 dzerr2=dzerr(goodind);
                [ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);
                %calculate errors
                [err1,dzerrh]=geterrors(beads2,ZcorrInterp,p);    
                dzerr(goodind)=dzerrh;
           end
           
           %ploto utput
           p.axhere=obj.initaxis('Interpolation');
           [ZcorrInterp]=getZinterp(beads2,ZcorrInterp,p);
           
           axhere=obj.initaxis('error');
           n=1:length(beads);
           plot(n,err0,n(goodind),err1,'*-');
           %correct beads for testing
           
           ax1=obj.initaxis('validation');
           f=ax1.Parent;
           ax2=axes(f,'Position',[0.5 0 1 1]);
           subplot(1,2,1,ax1);
           subplot(1,2,2,ax2);
       
            minv=inf;
            maxv=-inf;
           for k=length(beads):-1:1

               dZ=ZcorrInterp.interp(beads(k).loc.zglass,beads(k).loc.znm);
%                 z0glass=beads(k).loc.zglass-(beads(k).loc.frame)*p.dz
%                 dZ=ZcorrInterp.interp(z0glass,beads(k).loc.znm);
               beads(k).loc.znmcorrected=beads(k).loc.znm+dZ;
               
               if any(goodind==k)
                   col='k.';
                   inr=abs(beads(k).loc.z0relative)<1000;
                   if sum(inr)>0
                   minv=min(min(minv,min(beads(k).loc.znm(inr))),min(beads(k).loc.znmcorrected(inr)));
                   maxv=max(max(maxv,max(beads(k).loc.znm(inr))),max(beads(k).loc.znmcorrected(inr)));
                   end
               else
                   col='r.';
               end
               plot(ax1,beads(k).loc.z0relative,beads(k).loc.znm,col)
               plot(ax2,beads(k).loc.z0relative,beads(k).loc.znmcorrected,col)
               hold(ax1,'on');
               hold(ax2,'on');
           end
           
           xlim(ax1,[-1000 1000])
           ylim(ax1,[minv maxv]);
           xlim(ax2,[-1000 1000])
           ylim(ax2,[minv maxv]);
           
           file=obj.getPar('lastSMLFile');
           file=strrep(file,'_sml.mat','.mat');
           file=strrep(file,'.mat','_3Dcorr.mat');
           [f,p]=uiputfile(file);
           
           if f
               save([p f],'ZcorrInterp')
           end
            %get image stacks if needed
            
               % * zfitted(z0-zglass,ZObjective)
            %     * also for xfitted, yfitted, 
            %     * also spatially resolved
            % * z0-zglass=ztrue(zfitted, zObjective): interpolated or lookup table
            % * use for correction
            % * save with 3Dcal.mat
            % * fitter: instead of refractive index mismatch choose this correction.
            
            %save: either select existing 3Dcal:then it is appended. Or
            %save new (correction then with plugin)
           

        end
        
        function initGui(obj)
%             setvisible(0,0,obj)
%             beaddistribution_callback(0,0,obj)           
        end
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
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
%





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
tp=3.6;tmin=4.1;td=4.4;tmax=4.7;
w=0.3;
wp=0.5;
wcb=1.;

pard.beadsource.object=struct('String',{{'RoiManager','Segment'}},'Style','popupmenu','Value',2);
pard.beadsource.position=[1,1];
pard.beadsource.Width=1;

pard.dzt.object=struct('Style','text','String','dz (nm)'); 
pard.dzt.position=[2,1];
pard.dzt.Width=.4;
pard.dz.object=struct('Style','edit','String','10'); 
pard.dz.position=[2,1.4];
pard.dz.Width=.35;


% pard.alignz.object=struct('Style','checkbox','String','Align in z with f0','Value',1); 
% pard.alignz.position=[7,tp];
% pard.alignz.Width=1.2;


% pard.zrangeuset.object=struct('String','zrange (nm)','Style','text');
% pard.zrangeuset.position=[3,1];
% pard.zrangeuset.Width=1.5;
% pard.zrangeuse.object=struct('String','-800 800','Style','edit');
% pard.zrangeuse.position=[3,2.5];
% pard.zrangeuse.Width=1;

pard.smoothingt.object=struct('String','Smoothing (frame, zfit)','Style','text');
pard.smoothingt.position=[4,1];
pard.smoothingt.Width=1.5;
pard.smoothing.object=struct('String','0.02 0.0001','Style','edit');
pard.smoothing.position=[4,2.5];
pard.smoothing.Width=1;

pard.cutoffrefinet.object=struct('String','Maximum distance for outliers (nm)','Style','text');
pard.cutoffrefinet.position=[5,1];
pard.cutoffrefinet.Width=1.5;
pard.cutoffrefine.object=struct('String','150','Style','edit');
pard.cutoffrefine.position=[5,2.5];
pard.cutoffrefine.Width=1;

pard.setzero.object=struct('String','Set dz on glass to zero','Style','checkbox','Value',1);
pard.setzero.position=[6,1];
pard.setzero.Width=1.5;


pard.setglass.object=struct('String','Set glass to frame (empty: automatic): ','Style','checkbox','Value',0);
pard.setglass.position=[7,1];
pard.setglass.Width=3;

pard.glassframe.object=struct('String','','Style','edit');
pard.glassframe.position=[7,4];
pard.glassframe.Width=0.5;


% pard.spatialcalibration.object=struct('Style','checkbox','String','Spatial calibration','Value',0,'Callback',{{@setvisible,obj}}); 
% pard.spatialcalibration.position=[1,tp];
% pard.spatialcalibration.Width=wcb+.2;
% 
% pard.getcoords.object=struct('String','select','Style','pushbutton','Callback',{{@getcoords,obj}});
% pard.getcoords.position=[1,tp+wcb];
% pard.getcoords.Width=3*w+wp-wcb;
% pard.getcoords.Height=1;

% pard.tt1.object=struct('String','grid','Style','text');
% pard.tt1.position=[2,tp];
% pard.tt1.Width=wp;
% pard.mint.object=struct('String','min','Style','text');
% pard.mint.position=[2,tmin];
% pard.mint.Width=w;
% pard.dxt.object=struct('String','delta','Style','text');
% pard.dxt.position=[2,td];
% pard.dxt.Width=w;
% pard.maxt.object=struct('String','max','Style','text');
% pard.maxt.position=[2,tmax];
% pard.maxt.Width=w;
% 
% pard.Xt.object=struct('String','X (pix)','Style','text');
% pard.Xt.position=[3,tp];
% pard.Xt.Width=wp;
% 
% pard.Xmin.object=struct('String','0','Style','edit');
% pard.Xmin.position=[3,tmin];
% pard.Xmin.Width=w;
% pard.Xd.object=struct('String','512','Style','edit');
% pard.Xd.position=[3,td];
% pard.Xd.Width=w;
% pard.Xmax.object=struct('String','512','Style','edit');
% pard.Xmax.position=[3,tmax];
% pard.Xmax.Width=w;
% 
% pard.Yt.object=struct('String','Y (pix)','Style','text');
% pard.Yt.position=[4,tp];
% pard.Yt.Width=wp;
% 
% pard.Ymin.object=struct('String','0','Style','edit');
% pard.Ymin.position=[4,tmin];
% pard.Ymin.Width=w;
% pard.Yd.object=struct('String','256','Style','edit');
% pard.Yd.position=[4,td];
% pard.Yd.Width=w;
% pard.Ymax.object=struct('String','512','Style','edit');
% pard.Ymax.position=[4,tmax];
% pard.Ymax.Width=w;
% 
% pard.xyoverlapt.object=struct('String','x,y overlap (pix)','Style','text');
% pard.xyoverlapt.position=[5,tp];
% pard.xyoverlapt.Width=1;
% pard.xyoverlap.object=struct('String','10','Style','edit');
% pard.xyoverlap.position=[5,tmax];
% pard.xyoverlap.Width=w;
% 
% pard.zcalc.object=struct('String','z-dependent calibration','Style','checkbox','Value',1,'Callback',{{@setvisible,obj}});
% pard.zcalc.position=[6,tp];
% pard.zcalc.Width=1.5;
% 
% pard.Zt.object=struct('String','Z vals (nm)','Style','text');
% pard.Zt.position=[7,tp];
% pard.Zt.Width=wp;
% pard.Zval.object=struct('String',' 0:1000:3000','Style','edit');
% pard.Zval.position=[7,tmin];
% pard.Zval.Width=5-tp-wp;
% 
% pard.framerangecombinet.object=struct('String','z beyond interval (nm)','Style','text');
% pard.framerangecombinet.position=[8,tp];
% pard.framerangecombinet.Width=1;
% pard.framerangecombine.object=struct('String','100','Style','edit');
% pard.framerangecombine.position=[8,tmax];
% pard.framerangecombine.Width=w;







pard.inputParameters={'cam_pixelsize_nm'};
pard.plugininfo.type='ProcessorPlugin';


end
