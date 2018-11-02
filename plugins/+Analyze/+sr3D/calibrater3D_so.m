classdef calibrater3D_so<interfaces.DialogProcessor
    properties
        SXY
    end
    methods
        function obj=calibrater3D_so(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
            
            
%              obj.showresults=true;
%              obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-1;obj.guiPar.Vrim=obj.guiPar.Vrim-20;
        end
        
        function initGui(obj)
%             obj.run;
%             close(obj.handle)
        end
        function out=run(obj,p)
            
            fit3ddir=strrep(pwd,'SMAP','fit3D');
            if exist(fit3ddir,'file')
            addpath(fit3ddir);
            end
            
            fit3ddir=strrep(pwd,'SMAP','fit3Dcspline');
            if exist(fit3ddir,'file')
            addpath(fit3ddir);
            end
            
            smapdir=pwd;
           
                       
             %get bead positions
%              minframes=800/p.dz*1.5;
%              
%              if isempty(obj.locData.loc) %no files loaded
%                  posall={};
%                  posallnm={};
%                  rois={[0 0 512 512]};
%                  roi=rois{1};
%                  pixelsizes={[0.1 0.1]};
%                  files={};
%                  file='';
%              else
%                 locDatacopy=obj.locData.copy;
%                 locDatacopy.regroup(150,round(minframes/2));
%                 locg=locDatacopy.getloc({'xnm','ynm','numberInGroup','filenumber'},'layer',1,'Position','roi','removeFilter','filenumber','grouping','grouped');
%                 beadind=find(locg.numberInGroup>minframes);
%                 xb=locg.xnm(beadind);
%                 yb=locg.ynm(beadind);
%                 filenumber=locg.filenumber(beadind);
%                 files=[];
% 
%                  locfiles=obj.locData.files.file;
%                  minxy=[inf inf];maxxy=[0 0];
% 
%                  for k=1:length(locfiles) 
%                     if isfield(locfiles(k).info,'imagefile')&&~isempty(locfiles(k).info.imagefile)
%                         filename=locfiles(k).info.imagefile;
%                     else
%                         filename=locfiles(k).info.basefile;
%                     end
%                     filename=strrep(strrep(filename,'\',filesep),'/',filesep);
%                     file{k}.name=filename;
%                     file{k}.metadata=locfiles(k).info;
% %                     il=getimageloader(obj,filename);
% %                     file{k}.name=il.file;
% %                     file{k}.metadata=il.metadata;
% %                     infile=filenumber==k;
% %                     roi=il.metadata.roi;
% %                     pixelsize=il.metadata.cam_pixelsize_um;
% %                     file{k}.posx=xb(infile)/pixelsize(1)/1000-roi(1);
% %                     file{k}.posy=yb(infile)/pixelsize(end)/1000-roi(2);
% %                     il.close;
% %                     posall{k}=horzcat(file{k}.posx,file{k}.posy);
% %                     minxy=min(minxy,min(posall{k}));
% %                     maxxy=max(maxxy,max(posall{k}));
% %                     posallnm{k}=horzcat(xb,yb);
% %                     pixelsizes{k}=pixelsize;
% %                     rois{k}=roi;
%                  end
%                  files=obj.locData.files.file;
%              end
%              
%              smappos.positions=posall;
%              smappos.positionsnm=posallnm;
%              smappos.imageROI=roi;
%              
%              smappos.pixelsize=pixelsizes;
%              smappos.roi=rois;
%              smappos.files=files;
             smappos.P=obj.P;
             
             
%              if p.setrange
%                  if isempty(p.xrange)
%                      smappos.xrange=[minxy(1)-1 maxxy(1)+1];
%                  else
%                      smappos.xrange=p.xrange;
%                  end
%                  if isempty(p.yrange)
%                      smappos.yrange=[minxy(2)-1 maxxy(2)+1];
%                  else
%                      smappos.yrange=p.yrange;
%                  end                 
%              end
             
%              if p.setframerange
%                  smappos.framerangeuse=p.framerangeuse;
%              end
%              if exist('settings_3D.txt','file')
%                  smappos.settings=readstruct('settings_3D.txt');
%              else
%                  smappos.settings=[];
%              end
             smappos.fit3ddir=fit3ddir;
             smappos.smapdir=smapdir;
               cg=calibrate3D_GUI_g(smappos);    
%              switch p.global.selection
%                  case 'global'
%                      cg=calibrate3D_GUI_g(smappos);
%                  case '4pi'
%                      cg=calibrate3D_GUI_g(smappos);
%                  otherwise
%                     cg=calibrate3D_GUI(smappos); 
%              end
%              cg.guihandles.dz.String=num2str(p.dz);
%              if ~isempty(file)
%              cg.guihandles.filelist.String=getFieldAsVector(file,'name');
%              [p, f]=fileparts(file{1}.name);
%             
%              cg.guihandles.outputfile.String=[p  filesep f '_3dcal.mat'];
%              
%              if length(xb)>0
%                  cg.guihandles.posfromsmap.Value=true;
%              end
%              cg.guihandles.emgain.Value=file{1}.metadata.EMon;
%              end
             out=[];
        end
        
 
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function pard=guidef(obj)

% pard.dzt.object=struct('String','dz (nm)','Style','text');
% pard.dzt.position=[1,1];
% pard.dz.object=struct('String','10','Style','edit');
% pard.dz.position=[1,2];
% pard.global.object=struct('String',{{'single channel', 'global','4pi'}},'Style','popupmenu');
% pard.global.position=[1,3];
% pard.global.Width=2;
% 
% pard.setrange.object=struct('Style','checkbox','String','set range (empty = auto)');
% pard.setrange.position=[2,1];
% pard.setrange.Width=1.5;
% 
% pard.xrange.object=struct('String','','Style','edit');
% pard.xrange.position=[2,3];
% pard.yrange.object=struct('String','','Style','edit');
% pard.yrange.position=[2,4];
% 
% pard.setframerange.object=struct('Style','checkbox','String','set frame range ');
% pard.setframerange.position=[3,1];
% pard.setframerange.Width=1.5;
% 
% 
% pard.framerangeuse.object=struct('String','','Style','edit');
% pard.framerangeuse.position=[3,3];

pard.inputParameters={'cam_pixelsize_um'};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.name='calibrate3DsplinePSF';
pard.plugininfo.description=sprintf('Plugin to calibrate 3D PSF. \n According to Li et al, 2017');
end
