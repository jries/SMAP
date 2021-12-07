classdef DECODE_fitting<interfaces.WorkflowModule

    properties
        decodepid
        workingdir
        yamldefault='DECODE_local.yaml';
        imagefile
    end
    methods
        function obj=DECODE_fitting(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function prerun(obj,p)
            if p.overwritePixelsize
                obj.setPar('overwrite_pixelsize',[p.pixelsizex p.pixelsizey])
                cs=obj.getPar('loc_cameraSettings');
                cs.cam_pixelsize_um=[p.pixelsizex p.pixelsizey];
                obj.setPar('loc_cameraSettings',cs);
            else
                obj.setPar('overwrite_pixelsize',[])
            end
        end
        function outreturn=run(obj,data,p)
            outreturn=[];
            cam_settings=obj.getPar('loc_cameraSettings');
            [workingdirlocal, outname]=fileparts(p.outputpath);
            frameshere=cam_settings.imagefile;
            if strcmpi(p.runwhere.selection,'local')
                if ismac
                    gpurec='cpu';
                else
                    gpurec='cuda:0';
                end
                frame_path=cam_settings.imagefile;
                model_path=p.model_path;
                emitter_path=p.outputpath;
            else %server
                server=obj.getGlobalSetting('DECODE_server');
                gpustat=webread([server '/status_gpus']);
                [gpus,gpurec]=parsegpustathttp(gpustat);  
                            % frames for server 
                decodenetwork=obj.getGlobalSetting('DECODE_network_data');
                pstart=length(decodenetwork);
% %                 frameshere=cam_settings.imagefile;
                framesdir=frameshere(1:pstart);
                if ~strcmp(framesdir,decodenetwork) %not the same: later copy here
                    disp('images need to be on the decode network storage for data')
                    return
                end
                frame_path=[frameshere(pstart+2:end)];

                modeldir=p.model_path(1:pstart);
                if ~strcmp(modeldir,decodenetwork) %not the same: later copy here
                    disp('model need to be on the decode network storage for data')
                    return
                end            
                model_path=[p.model_path(pstart+2:end)];   

                outdir=p.outputpath(1:pstart);
                if ~strcmp(outdir,decodenetwork) %not the same: later copy here
                    disp('output need to be on the decode network storage for data')
                    return
                end  
                emitter_path=[p.outputpath(pstart+2:end)];
                yamlwrappathremote=[fileparts(emitter_path) '/' outname '_fitwrap.yaml'];
            end

            %for both
            wrapyaml.Hardware.device=gpurec;
            wrapyaml.Hardware.worker=(uint16(4));
            wrapyaml.Frames.path=frame_path;
            wrapyaml.Model.path=model_path;
            wrapyaml.Model.param_path=[fileparts(model_path) '/param_run.yaml'];
            wrapyaml.Output.path=emitter_path;

            % overwrite Camera parameters in this file: offset mirror,
            % gain,...
            if cam_settings.EMmirror
                  mirrorax=int16(-1);
            else
                  mirrorax=[];
            end
            wrapyaml.Camera.mirror_dim=mirrorax;
            wrapyaml.Camera.baseline=cam_settings.offset;
            wrapyaml.Camera.e_per_adu=cam_settings.conversion;
            wrapyaml.Camera.em_gain=cam_settings.emgain*cam_settings.EMon;
            if p.overwritePixelsize
                wrapyaml.Camera.px_size=[p.pixelsizex p.pixelsizey]*1000; %check XXXXX! could be swapped
            else
                wrapyaml.Camera.px_size=cam_settings.cam_pixelsize_um*1000;
            end

            %frames
            if obj.getPar('loc_preview')
                f1=max(0,obj.getPar('loc_previewframe')-2);
                f2=f1+2; 
                framerange=[f1 f2];
                starttext='Preview started, this can take a while (10s of seconds)';
            else
                frames=obj.getPar('loc_frames_fit');
                if frames(2)>=cam_settings.numberOfFrames 
                    if frames(1)==1
                        framerange=[];
                    else
                        framerange=[0 cam_settings.numberOfFrames];
                    end
                else
                    framerange=frames-1;
                end
                starttext='DECODE fitting started';
            end
            wrapyaml.Frames.range=uint32(framerange);

            yamlwrappathlocal=[workingdirlocal '/' outname '_fitwrap.yaml'];
            % make wrapper yaml
            WriteYamlSimple(yamlwrappathlocal, wrapyaml);
            fileh5=strrep(frameshere,'.tif','.h5'); % read h5
            %start fitting
            if strcmpi(p.runwhere.selection,'local')
                pdecode=obj.getGlobalSetting('DECODE_path');
%                 pcall='conda activate decode_env';
%                 pcall{2}=[pdecode 'python decode.neuralfitter.inference.infer --fit_meta_path ' yamlwrappathlocal];
                 pcall=[pdecode '/bin/python -m decode.neuralfitter.inference.infer --fit_meta_path ' yamlwrappathlocal];
                gitdecodepath=[fileparts(pwd) filesep 'DECODE'];
                logfile=[workingdirlocal '/' outname '_log.txt'];

                cpath=pwd;
                cd('../DECODE');
                [status,cmdout]=system([pcall '&>' logfile ' &']);
                cd(cpath);
%                 fi=dir(logfile);
%                 changed=fi.datenum;
                starttime=now;
                line="";
                while 1
                    pause(2)
                    if exist(logfile,'file')
                        alllines=readlines(logfile,'WhitespaceRule','trim','EmptyLineRule','skip');
                        if isempty(alllines)
                            continue
                        end
                        line=alllines(end);
                        if ~isempty(line)
                            obj.status(line);
                            drawnow
                        end
                        
                    end

                    %determine when to stop
                    if contains(line,"Fit done and emitters saved") 
                        disp('logfile contains line: fitting done')
                        break
                    end
                    
                    if exist(fileh5,'file') 
                        fi=dir(fileh5);
                        if fi.datenum>starttime
                            disp('h5 written')
                            break
                        end
                    end

%                     fi=dir(logfile);
%                     changednew=fi.datenum;
%                     if changednew == changed
% %                         pause(20)
%                         fi=dir(logfile);
%                         changednew=fi.datenum;    
%                         if changednew == changed
%                             disp('logfile did not change')
% %                             break
%                         end
%                     end
%                     changed=changednew;
                end   
            else %server
                % call decode fitter
                obj.status(starttext); drawnow;
                url = [obj.getGlobalSetting('DECODE_server') '/submit_fit'];
                options = weboptions('RequestMethod', 'post', 'ArrayFormat','json');
                pid = webread(url, 'path_fit_meta', yamlwrappathremote, options);
                obj.decodepid=pid;
                
                %update staus
                logfile=[workingdirlocal '/out.log' ];
                fittingstat=webread([server '/status_processes']);
                while strcmp(fittingstat.fit.(['x' num2str(pid)]),'running') || contains(fittingstat.fit.(['x' num2str(pid)]),'sleep')
                    pause(2)
                    if exist(logfile,'file')
                        alllines=readlines(logfile,'WhitespaceRule','trim','EmptyLineRule','skip');
                        if isempty(alllines)
                            continue
                        end
                        line=alllines(end);
                        if ~isempty(line)
                            obj.status(line);
                            drawnow
                        end
                    end
                    fittingstat=webread([server '/status_processes']);
                end
            end
           
            [locs,info]=decodeh5ToLoc(fileh5);
            locs.xpix=locs.xnm/info.pix2nm(1)+1;
            locs.ypix=locs.ynm/info.pix2nm(2)+2; %check ROI XXXXX
            locs.xpixerr=locs.xnmerr/info.pix2nm(1);
            locs.ypixerr=locs.ynmerr/info.pix2nm(2);

            if p.userefractive_index_mismatch %RI mismatch,
                locs.znm=locs.znm *p.refractive_index_mismatch;
                locs.locprecznm=locs.locprecznm *p.refractive_index_mismatch;
            end
            %check sign of z and pixel size. Look for offset compared to
            %fitting.

            obj.setPar('loc_fitinfo',wrapyaml) % setPer('fitinfo') all training and fitting yaml parameters,
            output=interfaces.WorkflowData;
            output.eof=true;
            output.data=locs;
            output.ID=1;
            output.frame=1;
            if obj.getPar('loc_preview')
                output.frame=framerange(1);
                obj.setPar('preview_locs',locs);
                obj.setPar('preview_peakfind',locs);
                obj.output(output,2);
                obj.status('Preview done.')             
            else
                obj.output(output);
                obj.status('DECODE fitting done.')
            end
        end
        function addFile(obj,file,setinfo)   %for batch processing?
            if isempty(file)
                fileinfo=obj.getPar('loc_fileinfo_set');
                file=fileinfo.imagefile;
            end
            if setinfo
                if obj.getSingleGuiParameter('runwhere').Value==2 || contains(file,'/decode/fits')
                    outfile=strrep(file,'.tif','.h5');
                else
                    [~,fn]=fileparts(file);
                    [~,dir]=fileparts(fileparts(file));
                    decodenetwork=obj.getGlobalSetting('DECODE_network_data');
                    outfile=[decodenetwork filesep 'fits' filesep dir filesep fn '.h5'];
                end
                obj.setGuiParameters(struct('outputpath',outfile))
                obj.imagefile=file;
                %later: selection if h5 or csv
            end
            %update output file
        %     obj.workingdir=fileparts(fileparts(p));
        %     if isempty(obj.getSingleGuiParameter('outputpath'))
        %         
        %     end
            %output on channel 2 data.data=filename
        end
        function setoutputfilename(obj)
            %can be dummy but called from batch processor
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.addSynchronization('loc_fileinfo_set',[],[],{@obj.addFile,[],1});
            % create global parameter: DECODE_network_data
            obj.createGlobalSetting('DECODE_path','DECODE','The anaconda environmet path of decode (eg. /decode_env/):',struct('Style','dir','String','decode')) 
            obj.createGlobalSetting('DECODE_network_data','DECODE','network directory for DECODE training and fitting',struct('Style','dir','String',' '))
            obj.createGlobalSetting('DECODE_server','DECODE','network directory for DECODE training and fitting',struct('Style','edit','String','http://pc-ries25:8000'))
        end            
    end
end

function loadmodel(a,b,obj)
if isempty(obj.workingdir)
    defaultp=obj.getGlobalSetting('DECODE_network_data');
else
    defaultp=fileparts(obj.getSingleGuiParameter('model_path'));
end
[f,p]=uigetfile([defaultp filesep 'model_*.pt']);
if f
    model_path=[p f];
    obj.setGuiParameters(struct('model_path',model_path))
    param_path=[p 'param_run.yaml'];
    if ~exist(param_path,'file')
        warndlg('param_run.yaml expected in model_0.pt path')
    end 
    obj.workingdir=fileparts(fileparts(p));
end
end


function selectoutput(a,b,obj)
p=obj.getAllParameters;
outputp=p.outputpath;
if ~isempty(obj.imagefile)
    [plocal,flocal]=fileparts(obj.imagefile);
    flocal=strrep(flocal,'.tif','.h5');
else
    flocal='fit.h5';
end

if isempty(outputp) || ~exist(fileparts(outputp),'dir')
    if p.runwhere.Value==2 %local
        
        outputp=plocal;
    else
        outputp=obj.getGlobalSetting('DECODE_network_data');
    end
    outputp=[outputp filesep flocal];
end

[file,pfad]=uiputfile(outputp);
if ~isempty(dir)
    obj.setGuiParameters(struct('outputpath',[pfad file]))
end
end

function training_callback(a,b,obj)
    name='DECODE training';
    module=plugin('Analyze','calibrate','DECODE_training_estimates');
    p.Vrim=100;
    module.handle=figure('MenuBar','none','Toolbar','none','Name',name);
    module.attachPar(obj.P);
    module.attachLocData(obj.locData);
    p.Xrim=10;
    module.setGuiAppearence(p)
    module.makeGui;
end

function pard=guidef(obj)

pard.loadcal.object=struct('Style','pushbutton','String','Load model','Callback',{{@loadmodel,obj}});
pard.loadcal.position=[1,1];
pard.loadcal.Width=1;
pard.model_path.object=struct('Style','edit','String','');
pard.model_path.position=[1,2];
pard.model_path.Width=3;
pard.model_path.TooltipString=sprintf('DECODE model');
% pard.modelyaml.object=struct('Style','edit','String','');
% pard.modelyaml.position=[2,1.75];
% pard.modelyaml.Width=1.5;

pard.select_output.object=struct('Style','pushbutton','String','Select output','Callback',{{@selectoutput,obj}});
pard.select_output.position=[2,1];
pard.select_output.Width=1;
pard.outputpath.object=struct('Style','edit','String','');
pard.outputpath.position=[2,2];
pard.outputpath.Width=3;
pard.outputpath.TooltipString=sprintf('DECODE output');


p(1).value=1; p(1).on={}; p(1).on={'server'};
p(2).value=2; p(2).off={'server'}; p(2).on={};

pard.runwhere.object=struct('Style','popupmenu','String',{{'Server','Local'}},'Callback',{{@obj.switchvisible,p}});
pard.runwhere.position=[3,1];
pard.runwhere.Width=1;
% pard.server.object=struct('Style','edit','String','riesgroup@pc-ries25:8000');
% pard.server.position=[3,2];
% pard.server.Width=3;


pard.starttraining.object=struct('Style','pushbutton','String','DECODE training','Callback',{{@training_callback,obj}});
pard.starttraining.position=[3,3];

p(1).value=0; p(1).on={}; p(1).off={'refractive_index_mismatch'};
p(2).value=1; p(2).on={'refractive_index_mismatch'}; p(2).off={};
pard.userefractive_index_mismatch.object=struct('Style','checkbox','String','RI mismatch:','Callback',{{@obj.switchvisible,p}});
pard.userefractive_index_mismatch.position=[5,3.5];
pard.userefractive_index_mismatch.Width=1.5;
pard.userefractive_index_mismatch.Optional=true;

pard.refractive_index_mismatch.object=struct('Style','edit','String','.8','Visible','off');
pard.refractive_index_mismatch.position=[5,4.5];
pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
pard.refractive_index_mismatch.Optional=true;
pard.refractive_index_mismatch.Width=0.5;


p(1).value=0; p(1).on={}; p(1).off={'pixelsizex','pixelsizey'};
p(2).value=1; p(2).on={'pixelsizex','pixelsizey'}; p(2).off={};
pard.overwritePixelsize.object=struct('Style','checkbox','String','New pixelsize X,Y (um):','Callback',{{@obj.switchvisible,p}});
pard.overwritePixelsize.position=[5,1];
pard.overwritePixelsize.Width=1.5;
pard.overwritePixelsize.Optional=true;

pard.pixelsizex.object=struct('Style','edit','String','.1','Visible','off');
pard.pixelsizex.position=[5,2.5];
pard.pixelsizex.Width=0.5;
pard.pixelsizex.Optional=true;

pard.pixelsizey.object=struct('Style','edit','String','.1','Visible','off');
pard.pixelsizey.position=[5,3];
pard.pixelsizey.Width=0.5;
pard.pixelsizey.Optional=true;



% pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};

pard.plugininfo.type='WorkflowModule';
pard.plugininfo.description='';
end