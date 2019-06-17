classdef SimulateCameraImagesP<interfaces.DialogProcessor
    properties
        locs
        par
    end
    methods
        function obj=SimulateCameraImagesP(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'se_sitefov','se_cellpixelsize','se_siteroi'};
            obj.history=true;
            obj.showresults=true;
        end
        
        function out=run(obj,p)  
          if obj.processorgui  
              [locs,p]=storelocs(obj,p);

              allframes=max(1,p.frames(1)):min(p.frames(end),max(locs.frame));
              lastsml=obj.getPar('lastSMLFile');
              if ~isempty(lastsml)
                path=fileparts(lastsml);
              else
                  path='';
              end
         
              [f,path]=uiputfile([path filesep '*.tif']);
              if f
                  p.plotaxis=obj.initaxis('camera image');
                [img,simulpar]=simulatecamera(locs,p,allframes);
                imout=uint16(img);
                 metadata = createMinimalOMEXMLMetadata(imout);
                 pixelSize = ome.units.quantity.Length(java.lang.Double(p.pixelsize/1000), ome.units.UNITS.MICROM);
                 metadata.setPixelsPhysicalSizeX(pixelSize, 0);
                 metadata.setPixelsPhysicalSizeY(pixelSize, 0);
%                  
                 metadata.setDetectorAmplificationGain(java.lang.Double(p.emgain),0,0);
                 metadata.setDetectorOffset(java.lang.Double(p.offset),0,0);
                 metadata.setDetectorGain(java.lang.Double(p.conversion),0,0);
%                  bfsave(imout,[path f],'XYTCZ','metadata',metadata);
                 if exist([path f],'file')
                     delete([path f]);
                 end
                 bfsave(imout,[path  f], 'metadata', metadata);
%                 saveastiff(imout,[path f])
              end
             
             
          elseif ~isempty(obj.parent)
                p=obj.par;   
                 if obj.getPar('loc_preview')
                    allframes=max(1,obj.getPar('loc_previewframe'));    
                 else
                    allframes=max(1,p.frames(1)):min(p.frames(end),max(obj.locs.frame));
                 end
              %inti
%               filestruct=initfile('simulation');
%               obj.setPar('loc_fileinfo',filestruct);
              
              for k=1:length(allframes)
                  data=interfaces.WorkflowData;
                  data.frame=allframes(k);
                  data.ID=k;
                  [img,simulpar]=simulatecamera(obj.locs,p,allframes(k));
                  data.data=img;
                  obj.parent.output(data);
              end
              data=interfaces.WorkflowData;
              data.eof=true;
              data.ID=k+1;
              data.frame=allframes(end)+1;
              obj.parent.output(data);
          end
        out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        
        function prerun(obj)
            if  obj.getPar('loc_preview')
                f='';path='';
            else
            [f , path]=uiputfile('simulation');
            end
            filestruct=initfile([path f]);
            p=obj.getAllParameters;
            [obj.locs,p]=storelocs(obj,p);
            [~,par]=simulatecamera(obj.locs,p,1);
            filestruct.info=interfaces.metadataSMAP;
            filestruct.info.cam_pixelsize_um=p.pixelsize([1 1])/1000;
           filestruct.info.Width=par.sizex;
           filestruct.info.Height=par.sizey;
           filestruct.info.offset=par.offset;
           filestruct.info.emgain=par.emgain;
           filestruct.info.conversion=par.conversion;
           filestruct.info.basefile=[path f];
           filestruct.info.numberOfFrames=max(obj.locs.frames);
           obj.par=par;
           obj.setPar('loc_fileinfo',filestruct.info);
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
          sim=ROIManager.Segment.SimulateSites;
          sim.attachPar(obj.P);
             sim.handle=figure('MenuBar','none','Toolbar','none','Name','simulate locs');
    p.Xrim=10;
    p.Vrim=100;
    sim.setGuiAppearence(p)
    sim.makeGui;
    disp('close simulate localization gui after caluclating localizations')
    waitfor(sim.handle)
    locs=sim.locData.loc;
%           sim.attachLocData(obj.locData);
          
          
  end
  
  if isfield(locs,'xnm_gt') && ~isempty(locs.xnm_gt)
        locs.x=locs.xnm_gt;
        locs.y=locs.ynm_gt;
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
          p.emgain=0;

      end
      if p.emgain==0
              p.EMon=false;
          else
              p.EMon=true;
      end

end

function pard=guidef(obj)

% 
pard.simulationsource.object=struct('String',{{'Use current localizations','Load localizations','Make with simulation plugin'}},'Style','popupmenu');
pard.simulationsource.position=[1,1];
pard.simulationsource.Width=2;

pard.t1.object=struct('String','Pixelsize (nm)','Style','text');
pard.t1.position=[2,1];
pard.t1.Width=0.75;

pard.pixelsize.object=struct('String','100','Style','edit');
pard.pixelsize.position=[2,1.75];
pard.pixelsize.Width=0.5;

pard.autorange.object=struct('String','auto','Style','checkbox','Value',1);
pard.autorange.position=[4,1];
pard.autorange.Width=0.5;

pard.t2.object=struct('String','X range','Style','text');
pard.t2.position=[4,1.5];
pard.t2.Width=0.5;

pard.xrange.object=struct('String','0 50000','Style','edit');
pard.xrange.position=[4,2];
pard.xrange.Width=0.5;

pard.t3.object=struct('String','Y range','Style','text');
pard.t3.position=[4,2.5];
pard.t3.Width=0.5;

pard.yrange.object=struct('String','0 50000','Style','edit');
pard.yrange.position=[4,3];
pard.yrange.Width=0.5;

pard.t4.object=struct('String','Background','Style','text');
pard.t4.position=[2,3];
pard.t4.Width=0.75;
pard.background.object=struct('String','20','Style','edit');
pard.background.Width=.5;
pard.background.position=[2,3.75];

pard.tf.object=struct('String','Frames','Style','text');
pard.tf.position=[4,3.5];
pard.tf.Width=0.5;
pard.frames.object=struct('String','1 inf','Style','edit');
pard.frames.Width=.5;
pard.frames.position=[4,4];

pard.usecam.object=struct('String','Convert to camera values','Style','checkbox','Value',1);
pard.usecam.position=[6,1];
pard.usecam.Width=2;

pard.t5.object=struct('String','EM gain (0=off)','Style','text');
pard.t5.position=[7,1];
pard.t5.Width=1;
pard.emgain.object=struct('String','100','Style','edit');
pard.emgain.Width=.5;
pard.emgain.position=[7,2];

pard.t6.object=struct('String','Conversion','Style','text');
pard.t6.position=[7,2.5];
pard.t6.Width=1;
pard.conversion.object=struct('String','5','Style','edit');
pard.conversion.Width=.5;
pard.conversion.position=[7,3.5];

pard.t7.object=struct('String','Offset','Style','text');
pard.t7.position=[7,4];
pard.t7.Width=0.5;
pard.offset.object=struct('String','100','Style','edit');
pard.offset.Width=.5;
pard.offset.position=[7,4.5];

% pard.save.object=struct('String','Save','Style','checkbox');
% pard.save.Width=1;
% pard.save.position=[1,4];


pard.plugininfo.type='ROI_Analyze';
end
