classdef SimulateCameraImages<interfaces.WorkflowModule
%     Plugin to simulate camera raw frames using a realistic fluorophore
%     blinking model, PSF model and camera noise model. Useful to test fitting
%     workflows.
    properties
        locs
        simulator
        par
        PSF
        file 
    end
    methods
        function obj=SimulateCameraImages(varargin)        
                obj@interfaces.WorkflowModule(varargin{:});
                obj.inputChannels=1;
                obj.isstartmodule=true;

        end
        function initGui(obj)
            simulation_callback(obj.guihandles.simulationsource, 0,obj)
            obj.setcampar;
            obj.makeinfobutton('ne');
%             psfpar_callback(0,0,obj,true)
        end
        function out=run(obj,data,p)  
                p=obj.par;   
                preview=obj.getPar('loc_preview');
                 if ~isfield(obj.locs,'bg')
                        obj.locs.bg=zeros(size(obj.locs.x));
                 end
                 if preview
                    allframes=max(1,obj.getPar('loc_previewframe'));  
                    indh=obj.locs.frame==allframes;

                    locgth=copystructReduce(obj.locs,indh,{'x','y','znm','bg','phot'});
                    locgt.x=(locgth.x-p.xrange(1))/p.pixelsize;
                    locgt.y=(locgth.y-p.yrange(1))/p.pixelsize;
                    locgt.z=locgth.znm;
                    locgt.N=locgth.phot;
                    locgt.bg=locgth.bg+p.background;
                    obj.setPar('loc_gt_preview',locgt)
                 else
                    allframes=max(1,p.frames(1)):min(p.frames(end),max(obj.locs.frame));
                 end
              allimgs=[];
              for k=1:length(allframes)
                  data=interfaces.WorkflowData;
                  data.frame=allframes(k);
                  data.ID=k;
                  [img,simulpar]=simulatecamera(obj.locs,p,allframes(k),obj.PSF); %obj.PSF not working
                  data.data=img;
                  if p.savetiffs &&~obj.getPar('loc_preview')
                      if isempty(allimgs)
                          sim=size(img);
                          allimgs=zeros(sim(1),sim(2),length(allframes),'single');
                      end
                      allimgs(:,:,k)=img;
                  end
                  obj.output(data);
                  obj.status(['simulating camera image frame ' num2str(k)]);
                  
              end
              data=interfaces.WorkflowData;
              data.eof=true;
              data.ID=k+1;
              data.frame=allframes(end)+1;
              obj.output(data);
              if obj.getPar('loc_preview')
                  locfit=obj.getPar('locdata_preview');
                  frameprev=locfit.frame(1);
                  locgt=copystructReduce(obj.locs,obj.locs.frame==frameprev);
              else
                  locfit=obj.locData.loc;
                  locgt=obj.locs;
              end
              locgt.bg=locgt.bg+p.background;
              %calculate CRLB, add to gt
              
              settings=obj.getPar('loc_cameraSettings');
              crlbfac=sqrt(settings.EMon+1);
              pixelsize=settings.cam_pixelsize_um*1000;
              psfmodel=obj.PSF;
              psfmodel.roisize=obj.getPar('loc_ROIsize');
              crlb=psfmodel.crlb(locgt.phot,locgt.bg,locgt.znm);  %why minus?
%                     crlb=psfmodel.crlb(lTn.phot,lTn.bg,-lTn.z);
                    locgt.xerr=sqrt(crlb(:,2))*pixelsize(1)*crlbfac;
                    locgt.yerr=sqrt(crlb(:,1))*pixelsize(end)*crlbfac;
                    locgt.zerr=sqrt(crlb(:,5))*crlbfac;
                    locgt.Nerr=sqrt(crlb(:,3))*crlbfac;
                    
              if ~obj.getPar('loc_preview')
                  figure
                simulationerror(locgt,locfit,1)
              end
              
              if p.savetiffs && ~obj.getPar('loc_preview')
                  saveim=uint16(allimgs);
                  pathh=fileparts(obj.getPar('loc_fileinfo').basefile);
                  [f,path]=uiputfile([pathh filesep 'simulation.tif']);
                  if f
                    saveastiff(saveim,[path f])
                    finfo=strrep(f,'.tif','_gt.mat');
                    save([path finfo], locgt);
                  end
              end
              
        out=[];
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        
        function prerun(obj,p)   
            if  obj.getPar('loc_preview')
                f='';path='';
            else
                if p.savetiffs
                    lf=obj.getPar('last_SMLFile');
                    if ~isempty(lf)
                        pin=fileparts(lf);
                    else
                        lf=obj.getPar('loc_fileinfo');
                        if ~isempty(lf.basefile)
                            pin=fileparts(lf.basefile);
                        else
                            pin=pwd;
                        end
                    end
                    [f , path]=uiputfile([pin filesep 'simulationfit_sml.mat']);
                else
                    path='';
                    f='simulation_nosave_sml.mat';
                end
                
            end
            obj.file=[path f];
%             obj.par.basefile
           obj.setcampar;
           %get PSF
           switch p.psfmodel.selection
               case {'Symm Gauss','Astig Gauss' }
                   obj.PSF=GaussPSF;
                   
                   if strcmp(p.psfmodel.selection,'Symm Gauss')
                       p.gausspar(3)=0;
                   end
                   obj.PSF.gausspar=p.gausspar;
%                    warndlg('not implemented')
                   
               case {'Spline'}
                   f=p.cal_3Dfile;
                   if isempty(f)
                       errordlg('please select first a 3D calibration file with a PSF model')
                       return
                   end
                   psf=splinePSF;
                   psf.loadmodel(f);
                   obj.PSF=psf;
                   roisize=obj.getPar('loc_ROIsize');
                   scoeff=size(psf.modelpar.coeff);
                   obj.PSF.roisize=min(scoeff(1)-1,roisize+4);
                   
           end
           [~,par]=simulatecamera(obj.locs,obj.par,1,obj.PSF);
           obj.par=par;
        end
        function setcampar(obj,a,b)
            cs=obj.getPar('loc_cameraSettings');
           
            p=obj.getAllParameters;

            if p.autorange
                     if isempty(obj.locs)
                         p.xrange=[0 1]; p.yrange=[0 1];
                     else
                      rim=1000;
                      p.xrange=[min(obj.locs.x)-rim max(obj.locs.x)+rim]; 
                      p.yrange=[min(obj.locs.y)-rim max(obj.locs.y)+rim]; 
                     end
            end
            obj.par=copyfields(p,cs);
            par=obj.par;
            xrp=(par.xrange/par.pixelsize);
            yrp=(par.yrange/par.pixelsize);        
            info=interfaces.metadataSMAP;
%             info.roi=round([xrp(1) yrp(1) xrp(2)-xrp(1) yrp(2)-yrp(1)]);
            info.roi=([xrp(1)+.5 yrp(1)+.5 xrp(2)-xrp(1) yrp(2)-yrp(1)]);
            info.cam_pixelsize_um=p.pixelsize/1000;
           
            if p.usecam
               info.offset=par.offset;
               if (par.cameramodel.Value==3)
                   info.emgain=par.emgain;
               else
                   info.emgain=1;
               end
               info.conversion=par.conversion;     
            else
               info.offset=0;
               info.emgain=1;
               info.conversion=1;
            end
           info.EMon=(par.emgain>0)&(par.cameramodel.Value==3);
           info.basefile=obj.file;
           if ~isempty(obj.locs)
           info.numberOfFrames=max(obj.locs.frame);
           end
           info.assigned.roi=true;
           obj.setPar('loc_fileinfo',info);
           obj.setPar('loc_fileinfo_set',info);
           obj.par=copyfields(obj.par,info,{'offset','emgain','conversion','EMon'});
            
        end
    end
end

function [locs,p]=storelocs(obj,p)
  switch p.simulationsource.Value
      case 1 %current locs
          locs=obj.locData.getloc({'xnm','ynm','znm','xnm_gt','ynm_gt','frame','phot'},'position','roi','layer',1);

      case 2 %load localizations
          path=fileparts(obj.getPar('lastSMLFile'));
          [f,path]=uigetfile([path filesep '*.mat']);
          if ~f
              return
          end
          l=load([path f]);
          locs=l.saveloc.loc;
      case 3
            locs=obj.simulator.locData.loc;
            obj.simulator.locData.clear;
      case 4
          locs=randomlocs(p);
  end
  
  if isfield(locs,'xnm_gt') && ~isempty(locs.xnm_gt)
        locs.x=locs.xnm_gt;
        locs.y=locs.ynm_gt;
        if isfield(locs,'znm_gt')
            locs.znm=locs.znm_gt;
        end
  else
      locs.x=locs.xnm;
      locs.y=locs.ynm;
  end
          
    if p.autorange
      p.xrange=[min(locs.x)-1000 max(locs.x)+1000];
      p.yrange=[min(locs.y)-1000 max(locs.y)+1000];
      p.frames=[min(locs.frame) max(locs.frame)];
    end
     if ~p.usecam
          p.EMon=false;
          p.conversion=1; 
          p.offset=0;
          p.emgain=1;

      end
      if p.emgain==0
              p.EMon=false;
      else
              p.EMon=true;
      end

end

function locs=randomlocs(p)
[settings, button] = settingsdlg(...
    'Description', 'Define simulation parameters',...
    {'Density µm^-2';'density'}, 0.02,...
    {'Photons';'photons'}, 5000,...
    {'Frames';'frames'}, 5000,...
    {'Lifetime (frames)';'lifetime'}, 1,...
    {'Zrnage (nm))';'zrange'}, 1000,...
    {'Z offset (nm))';'zoffset'}, 0);
if ~strcmpi(button,'ok')
    locs=[];
    return
end
if p.autorange
    p.xrange=[0 10000]; p.yrange=[0 10000];
end

area=(p.xrange(2)-p.xrange(1))*(p.yrange(2)-p.yrange(1))/1e6; %in um2
locsperframe=area*settings.density;
totallocs=locsperframe*settings.frames;
startf=rand(totallocs,1)*settings.frames;
lifet=exprnd(settings.lifetime,totallocs,1);
peakint=settings.photons/settings.lifetime;

ind=1;
if settings.lifetime==0
    lttot=1;
else
    lttot=settings.lifetime*1.5;
end
nloc=ceil(totallocs*lttot);
frames=zeros(nloc,1);
photons=zeros(nloc,1);
coord=zeros(nloc,3);
for k=1:totallocs
    onf=ceil(startf(k));
    offf=ceil(startf(k)+lifet(k));
    ind0=ind;
    frames(ind)=onf;
    if settings.lifetime==0
        photons(ind)=settings.photons;
    else
        photons(ind)=(onf-startf(k))*peakint;
    end
    ind=ind+1;
    while offf-onf>1
        onf=onf+1;
        frames(ind)=onf;
        photons(ind)=peakint;
        ind=ind+1;
    end
    if offf-onf>0 %last frame    
        frames(ind)=onf+1;
        photons(ind)=(startf(k)+lifet(k)-offf+1)*peakint;
        ind=ind+1;
    end 
    coord(1,ind0:ind-1)=rand(ind-ind0,1)*(p.xrange(2)-p.xrange(1))+p.xrange(1);
    coord(2,ind0:ind-1)=rand(ind-ind0,1)*(p.yrange(2)-p.yrange(1))+p.yrange(1);
    coord(3,ind0:ind-1)=rand(ind-ind0,1)*settings.zrange(1)-settings.zrange(1)/2+settings.zoffset;
    
end
frames(ind:end)=[];
photons(ind:end)=[];
coord(:,ind:end)=[];
[~,indsort]=sort(frames);
locs.xnm_gt=coord(1,indsort)';
locs.ynm_gt=coord(2,indsort)';
locs.znm_gt=coord(3,indsort)';
locs.frame=frames(indsort);
locs.phot=photons(indsort);

end

function getlocs_callback(a,b,obj)
p=obj.getAllParameters;
[obj.locs,obj.par]=storelocs(obj,p);
filestruct=initfile('');
filestruct.info=interfaces.metadataSMAP;
filestruct.info.numberOfFrames=max(obj.locs.frame);
obj.setPar('loc_fileinfo',filestruct.info);
obj.setPar('status','localizations transferred');
end

function simulation_callback(a,b,obj)
if a.Value==3 %make with simulation plugin
    if isempty(obj.simulator) || isempty(obj.simulator.handle)||~isvalid(obj.simulator.handle)
        obj.simulator=ROIManager.Segment.SimulateSites;
%         obj.simulator.attchLocData(obj.locData);
        obj.simulator.attachPar(obj.P);
        obj.simulator.locData.attachPar(obj.P);

%     end
    obj.simulator.locData.clear;
%     if isempty(obj.simulator.handle)||~isvalid(obj.simulator.handle)
        p.Xrim=10;
        p.Vrim=100;
        obj.simulator.setGuiAppearence(p)
        obj.simulator.handle=figure('MenuBar','none','Toolbar','none','Name','simulate locs');
        obj.simulator.makeGui;
        disp('press get localizations after running the simulator')
        obj.children.simulator=obj.simulator;
    end
    figure(obj.simulator.handle)
    
end
end

function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
if isempty(p.cal_3Dfile)
    path=obj.getGlobalSetting('DataDirectory');
    fh=obj.getPar('loc_fileinfo');
    if ~isempty(fh) && ~isempty(fh.imagefile)
        path=fileparts(fh.imagefile);
    end  
    p.cal_3Dfile=[path filesep '*3dcal.mat'];
end
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY') && ~isfield(l,'cspline')
        msgbox('no 3D data recognized. Select other file.');
    end
    obj.setGuiParameters(struct('cal_3Dfile',[p f]));
        obj.setPar('cal3Dfile',[p f]);
    
end
end


function pard=guidef(obj)

% 
pard.simulationsource.object=struct('String',{{'Use current localizations','Load localizations','Make with simulation plugin','Random'}},'Style','popupmenu',...
    'Callback',{{@simulation_callback,obj}});
pard.simulationsource.position=[1,1];
pard.simulationsource.Width=1.8;

pard.getlocalizations.object=struct('String','get Localizations','Style','pushbutton','Callback',{{@getlocs_callback,obj}},'FontWeight','bold');
pard.getlocalizations.position=[1,2.8];
pard.getlocalizations.Width=1.2;
p(1).on={'gausspar','gausspart'};p(1).value=1;p(1).off={'cal_3Dfile','loadcal'};
p(2)=p(1);p(2).value=2;
p(3).value=3;p(3).on={'cal_3Dfile','loadcal'};p(3).off={'gausspar','gausspart'};
pard.psfmodel.object=struct('String',{{'Symm Gauss','Astig Gauss','Spline'}},'Style','popupmenu','Callback',{{@obj.switchvisible,p}},'Value',3);
pard.psfmodel.position=[2,1];
pard.psfmodel.Width=1.5;

pard.gausspar.object=struct('String','1 500 200','Style','edit','Visible','off');
pard.gausspar.position=[2,2.5];
pard.gausspar.Width=1.;

pard.gausspart.object=struct('String','sig (pix), zR (nm), d(nm)','Style','text','Visible','off');
pard.gausspart.position=[2,3.5];
pard.gausspart.Width=1.5;

pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[2,2.5];
pard.cal_3Dfile.Width=2;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');

pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal','Callback',{{@loadcall_callback,obj}});
pard.loadcal.position=[2,4.5];
pard.loadcal.Width=0.5;

% pard.getpsfpar.object=struct('String','PSF model','Style','pushbutton','Callback',{{@psfpar_callback,obj}});
% pard.getpsfpar.position=[2,2];
% pard.getpsfpar.Width=2;

lp=3;

pard.savetiffs.object=struct('String','Save Tiff files','Style','checkbox','Value',0);
pard.savetiffs.position=[1,4];
pard.savetiffs.Width=1;

%camera model
p(1).on={};p(1).value=1;p(1).off={'scmosnoiset','scmosnoise'};
p(2).value=2;p(2).on={'scmosnoiset','scmosnoise'};p(2).off={};
p(3)=p(2);p(3).value=3;
pard.cameramodel.object=struct('String',{{'perfect cam','sCMOS','EMCCD'}},'Style','popupmenu','Callback',{{@obj.switchvisible,p}},'Value',1);
pard.cameramodel.position=[lp,1];
pard.cameramodel.Width=1.15;

pard.scmosnoiset.object=struct('String','noise e-','Style','text','Visible','off');
pard.scmosnoiset.position=[lp,2.15];
pard.scmosnoiset.Width=0.5;
pard.scmosnoise.object=struct('String','1.5','Style','edit','Visible','off');
pard.scmosnoise.position=[lp,2.65];
pard.scmosnoise.Width=0.25;

pard.t1.object=struct('String','Pixelsize (nm)','Style','text');
pard.t1.position=[lp,2.9];
pard.t1.Width=1;

pard.pixelsize.object=struct('String','100','Style','edit');
pard.pixelsize.position=[lp,3.7];
pard.pixelsize.Width=0.3;

pard.t4.object=struct('String','Background','Style','text');
pard.t4.position=[lp,4];
pard.t4.Width=0.7;
pard.background.object=struct('String','20','Style','edit');
pard.background.Width=.3;
pard.background.position=[lp,4.7];


la=4;
p=[];
p(1).value=0;p(1).on={'t2','xrange','t3','yrange','tf','frames'};p(1).off={};
p(2).value=1;p(2).off={'t2','xrange','t3','yrange','tf','frames'};p(2).on={};
pard.autorange.object=struct('String','auto','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p,{@obj.setcampar}}});
pard.autorange.position=[la,1];
pard.autorange.Width=0.5;

pard.t2.object=struct('String','X range','Style','text','Visible','off');
pard.t2.position=[la,1.5];
pard.t2.Width=0.5;

pard.xrange.object=struct('String','0 50000','Style','edit','Visible','off','Callback',@obj.setcampar);
pard.xrange.position=[la,2];
pard.xrange.Width=0.5;

pard.t3.object=struct('String','Y range','Style','text','Visible','off','Callback',@obj.setcampar);
pard.t3.position=[la,2.5];
pard.t3.Width=0.5;

pard.yrange.object=struct('String','0 50000','Style','edit','Visible','off','Callback',@obj.setcampar);
pard.yrange.position=[la,3];
pard.yrange.Width=0.5;

pard.tf.object=struct('String','Frames','Style','text','Visible','off');
pard.tf.position=[la,3.5];
pard.tf.Width=0.5;

pard.frames.object=struct('String','1 inf','Style','edit','Visible','off');
pard.frames.Width=.5;
pard.frames.position=[la,4];

lu=5;
p(1).value=0;p(1).on={};p(1).off={'t5','emgain','t6','conversion','t7','offset'};
p(2).value=1;p(2).off={};p(2).on={'t5','emgain','t6','conversion','t7','offset'};

pard.usecam.object=struct('String','Use camera units','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p,{@obj.setcampar}}});
pard.usecam.position=[lu,1];
pard.usecam.Width=1.3;

pard.t5.object=struct('String','EMgain','Style','text');
pard.t5.position=[lu,2.5];
pard.t5.Width=.8;
pard.emgain.object=struct('String','100','Style','edit','Callback',@obj.setcampar);
pard.emgain.Width=.3;
pard.emgain.position=[lu,3];

pard.t6.object=struct('String','Conversion','Style','text');
pard.t6.position=[lu,3.3];
pard.t6.Width=.7;
pard.conversion.object=struct('String','5','Style','edit','Callback',@obj.setcampar);
pard.conversion.Width=.3;
pard.conversion.position=[lu,4];

pard.t7.object=struct('String','Offset','Style','text');
pard.t7.position=[lu,4.3];
pard.t7.Width=0.4;
pard.offset.object=struct('String','100','Style','edit','Callback',@obj.setcampar);
pard.offset.Width=.3;
pard.offset.position=[lu,4.7];

% pard.save.object=struct('String','Save','Style','checkbox');
% pard.save.Width=1;
% pard.save.position=[1,4];

pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};

pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='Plugin to simulate camera raw frames using a realistic fluorophore blinking model, PSF model  and camera noise model. Useful to test fitting workflows.';
end
