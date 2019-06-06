classdef TifLoader<interfaces.WorkflowModule
    % TifLoader
    %     Loads raw camera images (single file or direcotry with single
    %     images). These can be tiff, files, but also any OME compatible image
    %     files. This plugin can load and process images during acquisition for
    %     online fitting. With the help of the CameraManager metadata is parsed
    %     and passed on to the CameraConverter.';

    properties   
        mirrorem
        loaders
        imloader
        framestop
        framestart
        numberOfFrames
        timerfitstart
        edgesize
    end
    methods
        function obj=TifLoader(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1;
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
             obj.loaders={'auto',@imageloaderAll;'MMStack',@imageloaderMM;'MMSingle',@imageloaderMMsingle;'OME',@imageloaderOME;'simple Tif',@imageloaderTifSimple};
            initGui@interfaces.WorkflowModule(obj);
            obj.inputParameters={'loc_subtractbg','loc_blocksize_frames'};            
            obj.guihandles.loadtifbutton.Callback={@loadtif_callback,obj};
            obj.addSynchronization('filelist_localize',obj.guihandles.tiffile,'String',{@loadtif_ext,obj});
            obj.guihandles.loaderclass.String=obj.loaders(:,1);
        end
        function prerun(obj,p)
            tf=p.tiffile;
            if iscell(tf)
                tf=tf{1};
            end
            if ~exist(tf,'file')
                obj.status('TifLoader: localization file not found')
                 error('TifLoader: localization file not found')
            end
%      
            if isempty(obj.imloader)
                obj.addFile(p.tiffile)
            else
                try 
                    img=(obj.imloader.getimage(1));
                    if isempty(img)
                        obj.addFile(p.tiffile)
                    end
                catch
                    disp('reload file in image loader')
                    obj.addFile(p.tiffile);
                end
               
            end

            obj.imloader.onlineAnalysis=p.onlineanalysis;
            obj.imloader.waittime=p.onlineanalysiswaittime;
            obj.imloader.setImageNumber(p.framestart-1);
            obj.numberOfFrames=obj.imloader.metadata.numberOfFrames;

            if p.onlineanalysis
                obj.framestop=inf;
            else
                obj.framestop=p.framestop;
            end
            obj.framestart=p.framestart;
            if obj.getPar('loc_preview')
                previewframe=obj.getPar('loc_previewframe');               
                if p.loc_subtractbg
                dt=p.loc_blocksize_frames;
                frameload=max(1,previewframe-floor(dt/2));
                obj.framestart=frameload;
                obj.imloader.setImageNumber(frameload-1);
                obj.framestop=frameload+dt-1;
                else
                    obj.imloader.setImageNumber(previewframe-1);
                    obj.framestop=previewframe;
                end               
            end
            obj.timerfitstart=tic;   
            obj.setPar('savefit',0) %delete, reset
        end
        function run(obj,data,p)
            global SMAP_stopnow
            if SMAP_stopnow
                disp('STOP button pressed. To localize, unpress')
            end
            
            if nargin>1&& ~isempty(data)&&~isempty(data.data) %optional input channel
                file=data.data;
                obj.addFile(file)
            end
            id=1;

            imloader=obj.imloader;
            image=imloader.readNext ; 

            tall=0;
            tfitall=0;
            
            % parallel
            if p.parallelload
                parp=gcp;
            end
            allimages=0;
            imcounter=0;
            while ~isempty(image)&&imloader.currentImageNumber<=obj.framestop&&~SMAP_stopnow                
                datout=interfaces.WorkflowData;

                if obj.mirrorem
                    image=image(:,end:-1:1);
                end
                allimages=double(image)+allimages;
                imcounter=imcounter+1;
                datout.data=image;

                datout.frame=imloader.currentImageNumber;
                datout.ID=id;
                id=id+1;
                th=tic;
                obj.output(datout)
                tfitall=tfitall+toc(th);

      

                
                %display
                if mod(datout.frame,10)==0
                    obj.setPar('loc_currentframe',struct('frame',datout.frame,'image',image));
                    if p.onlineanalysis
                        
                        totalf=imloader.metadata.numberOfFrames-obj.framestart;
                        elapsed=toc(obj.timerfitstart);
                        totaltime=elapsed/(datout.frame-obj.framestart+1)*totalf;
                    else
                        totalf=min(obj.numberOfFrames,obj.framestop)-obj.framestart;
                        elapsed=toc(obj.timerfitstart);
                        totaltime=elapsed/(datout.frame-obj.framestart+1)*totalf;
                    end
                    
                    numlocs=obj.getPar('fittedLocs');
                    sep='''';
                    statuss=['frame ' separatethousands(datout.frame-obj.framestart,sep) ' of ' separatethousands(totalf,sep) ...
                        '. Time: ' separatethousands(elapsed,sep,0) ' of ' separatethousands(totaltime,sep,0) sprintf('s.')...
                        separatethousands(numlocs,sep) ' locs, ' separatethousands(numlocs/elapsed,sep,0) ' locs/s.'];
                    obj.status(statuss);
                end
                th=tic;
                image=imloader.readNext;
                tall=tall+toc(th);
            end
            
            obj.setPar('tiffloader_loadingtime',tall);
            obj.setPar('tiffloader_fittime',tfitall);
            obj.setPar('tiffloader_averagetiff',cast(allimages/imcounter,'like',image));
            dateof=interfaces.WorkflowData;
            dateof.frame=imloader.currentImageNumber+1;
            dateof.ID=id;
            dateof.eof=true;
            obj.output(dateof)
            if ~obj.getPar('loc_preview') %if real fitting: close
                obj.imloader.close;
                disp('close imloader')
            end
            disp('fitting done')
        end
        function addFile(obj,file,setinfo)
            if nargin<3
                setinfo=false;
            end
            try
                obj.imloader.close;
            catch err
                
            end
             lc=obj.getSingleGuiParameter('loaderclass');
            loader=obj.loaders{lc.Value,2};
            obj.imloader=loader(file,obj.getPar('loc_fileinfo'),obj.P);
            
            if ~isempty(obj.imloader.metadata.allmetadata)&&isfield(obj.imloader.metadata.allmetadata,'metafile')
             obj.setPar('loc_metadatafile',obj.imloader.metadata.allmetadata.metafile);
            end
            allmd=obj.imloader.metadata;
            warnmissingmeta(allmd);

            p=obj.getGuiParameters;
            fileinf=obj.imloader.metadata;
            if fileinf.EMon && p.mirrorem  %if em gain on and mirrorem on: switch roi
                %It seems that on the Andor the roi is independent on the
                %mode, 
%                 if any(fileinf.roi(1:2)>0) %if roi(1:2)=[0 0] it is likely that roi was not read out and set to default.
%                     fileinf.roi(1)=512-fileinf.roi(1)-fileinf.roi(3);
%                 end
                fileinf.EMmirror=true;
            else 
                fileinf.EMmirror=false;
            end
            obj.mirrorem=fileinf.EMmirror;

            obj.setPar('loc_fileinfo',fileinf);
            obj.setPar('loc_filename',file);
            
            if iscell(obj.imloader.file)
                obj.guihandles.tiffile.Max=100;
            else
                obj.guihandles.tiffile.Max=1;
            end
            obj.guihandles.tiffile.String=obj.imloader.file;

             p=obj.getAllParameters;
             if p.onlineanalysis
                 obj.guihandles.framestop.String='inf';
             else
                numf=obj.imloader.metadata.numberOfFrames;
                if isnan(numf)
                    numf=inf;
                end
                obj.guihandles.framestop.String=int2str(numf);
             end
             if setinfo
                 obj.setPar('loc_fileinfo_set',obj.getPar('loc_fileinfo'));
             end

        end
        function loadtif_ext(obj)
            p=obj.getAllParameters;
            obj.addFile(p.tiffile,true)
            
        end
        function loadedges(obj,a,b)
            tiffile=obj.getSingleGuiParameter('tiffile');
            if exist(tiffile,'file')
                obj.addFile(tiffile);
            end
        end
        function setoutputfilename(obj)
            outfile=[obj.getPar('loc_fileinfo').basefile];
            [path, file]=fileparts(outfile);
            obj.setPar('loc_outputfilename',[path filesep file  '_sml.mat'])
        end
    end
end

function loadtif_callback(a,b,obj)
p=obj.getGuiParameters;

if p.ismultifile %later: check filename (e.g. _q1 _q2 etc). Also make sure quadrants are not mixed /rearranged
    if iscell(p.tiffile)
        p.tiffile=p.tiffile{1};
    end
    sf=selectManyFiles(fileparts(p.tiffile));
    filelisth=obj.guihandles.tiffile.String;
    if ~iscell(filelisth)
        filelisth={filelisth};
    end
    sf.guihandles.filelist.String=filelisth;
    waitfor(sf.handle);
    f=sf.filelist;
else
    try
        fe=bfGetFileExtensions;
    catch
        fe='*.*';
    end
    filelisth=p.tiffile;
    if iscell(filelisth)
        filelisth=filelisth{1};
    end
    [f,path]=uigetfile(fe,'select camera images',[fileparts(filelisth) filesep '*.tif']);
    if f
        f=[path f];
    else 
        f=[];
    end
end
if ~isempty(f)
    obj.addFile(f);
end  
%set output filename
outfile=[obj.getPar('loc_fileinfo').basefile];
[path, file]=fileparts(outfile);
obj.setPar('loc_outputfilename',[path filesep file  '_sml.mat'])
obj.setPar('loc_fileinfo_set',obj.getPar('loc_fileinfo'));
end

function warnmissingmeta(md)
expected={'emgain','conversion','offset','EMon','cam_pixelsize_um'};
missing=setdiff(expected,fieldnames(md));
if isempty(missing)
    return
end
str={'essential metadata missing (set manually in camera settings): ', missing{:}};
warndlg(str)
end

function changeloader(a,b,obj)
file=obj.getSingleGuiParameter('tiffile');
try
addFile(obj,file)
catch err
    obj.status('error in loading file. Choose different tiff loader')
    err
end  
end

function mirrorem_callback(a,b,obj)
if ~isempty(obj.imloader)
     fileinf=obj.imloader.metadata;
            if fileinf.EMon && obj.getSingleGuiParameter('mirrorem')  %if em gain on and mirrorem on: switch roi
                fileinf.EMmirror=true;
            else 
                fileinf.EMmirror=false;
            end
            obj.mirrorem=fileinf.EMmirror;
            if obj.getSingleGuiParameter('padedges') 
                dr=obj.getSingleGuiParameter('padedgesdr');
                fileinf.roi(1:2)=fileinf.roi(1:2)-dr;
                fileinf.roi(3:4)=fileinf.roi(3:4)+2*dr;
                obj.edgesize=dr;
            end
            obj.setPar('loc_fileinfo',fileinf);
end
end

function pard=guidef(obj)
pard.text.object=struct('Style','text','String','Image source file:');
pard.text.position=[1,1];
pard.text.Width=1.5;
pard.text.Optional=true;

pard.tiffile.object=struct('Style','edit','String',' ','HorizontalAlignment','right','Max',100);
pard.tiffile.position=[2,1];
pard.tiffile.Width=4;

pard.loadtifbutton.object=struct('Style','pushbutton','String','load images','Visible','on');
pard.loadtifbutton.position=[3,1];
pard.loadtifbutton.TooltipString=sprintf('Open raw camera image tif files. \n Either single images in directory. \n Or multi-image Tiff stacks');

pard.loaderclass.object=struct('Style','popupmenu','String','auto','Callback',{{@changeloader,obj}});
pard.loaderclass.position=[3,4];
pard.loaderclass.Width=1;
pard.loaderclass.Optional=true;

pard.ismultifile.object=struct('Style','checkbox','String','channels in different files','Value',0);
pard.ismultifile.position=[3,2];
pard.ismultifile.Width=2;
pard.ismultifile.Optional=true;


pard.onlineanalysis.object=struct('Style','checkbox','String','Online analysis. Waittime (s):','Value',0);
pard.onlineanalysis.position=[4,2];
pard.onlineanalysis.Width=1.75;
pard.onlineanalysis.TooltipString='Fit during acquisition. If checked, max frames is ignored. Waits until no more images are written to file.';
pard.onlineanalysis.Optional=true;
pard.onlineanalysiswaittime.object=struct('Style','edit','String','5');
pard.onlineanalysiswaittime.position=[4,3.75];
pard.onlineanalysiswaittime.Width=0.5;
pard.onlineanalysiswaittime.TooltipString='Waiting time for online analysis.';
pard.onlineanalysiswaittime.Optional=true;

pard.parallelload.object=struct('Style','checkbox','String','parallel','Value',0);
pard.parallelload.position=[4,4.25];
pard.parallelload.Width=0.75;
pard.parallelload.TooltipString='Parallel load and process. Makes sense only for very slow loading process.';
pard.parallelload.Optional=true;


pard.textf.object=struct('Style','text','String','Frame range');
pard.textf.position=[5.2,1.25];
pard.textf.Width=0.75;
pard.textf.Optional=true;
pard.framestart.object=struct('Style','edit','String','1');
pard.framestart.position=[5.2,2];
pard.framestart.Width=0.5;
pard.framestart.Optional=true;
pard.framestop.object=struct('Style','edit','String','1000000');
pard.framestop.position=[5.2,2.5];
pard.framestop.Width=0.5;
pard.framestop.Optional=true;

pard.mirrorem.object=struct('Style','checkbox','String','EM mirror','Value',1,'Callback',{{@mirrorem_callback,obj}});
pard.mirrorem.position=[5.2,4.2];
pard.mirrorem.TooltipString=sprintf('calibrate gain and offset from images');
pard.mirrorem.Optional=true;
pard.mirrorem.Width=0.8;

pard.plugininfo.type='WorkflowModule'; 
t1='Loads raw camera images (single file or direcotry with single images). These can be tiff, files, but also any OME compatible image files. This plugin can load and process images during acquisition for online fitting. With the help of the CameraManager metadata is parsed and passed on to the CameraConverter.';
pard.plugininfo.description=t1;
end