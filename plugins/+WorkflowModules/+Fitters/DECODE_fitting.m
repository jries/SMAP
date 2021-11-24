classdef DECODE_fitting<interfaces.WorkflowModule

    properties
        decodepid
        workingdir
        yamldefault='DECODE_local.yaml';
    end
    methods
        function obj=DECODE_fitting(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1; 
            obj.isstartmodule=true;
%              
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function prerun(obj,p)
           
        end
        function run(obj,data,p)
            cam_settings=obj.getPar('loc_cameraSettings');
            server=p.server;
            gpustat=webread([server '/status_gpus']);
            [gpus,gpurec]=parsegpustathttp(gpustat);   
            wrapyaml.Hardware.device=gpurec;
            wrapyaml.Hardware.worker=(uint16(4));

            % frames for server 
            decodenetwork=obj.getGlobalSetting('DECODE_network_data');
            pstart=length(decodenetwork);
            frameshere=cam_settings.imagefile;
            framesdir=frameshere(1:pstart);
            if ~strcmp(framesdir,decodenetwork) %not the same: later copy here
                disp('images need to be on the decode network storage for data')
                return
            end
            frame_path=[frameshere(pstart+2:end)];
            wrapyaml.Frames.path=frame_path;

            modeldir=p.model_path(1:pstart);
            if ~strcmp(modeldir,decodenetwork) %not the same: later copy here
                disp('model need to be on the decode network storage for data')
                return
            end            
            model_path=[p.model_path(pstart+2:end)];
            wrapyaml.Model.path=model_path;
            wrapyaml.Model.param_path=[fileparts(model_path) '/param_run.yaml'];
            
            outdir=p.outputpath(1:pstart);
            if ~strcmp(outdir,decodenetwork) %not the same: later copy here
                disp('output need to be on the decode network storage for data')
                return
            end  
            emitter_path=[p.outputpath(pstart+2:end)];
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
            end
            [workingdirlocal, outname]=fileparts(p.outputpath);
            yamlwrappathlocal=[workingdirlocal '/' outname '_fitwrap.yaml'];
            yamlwrappathremote=[fileparts(emitter_path) '/' outname '_fitwrap.yaml'];
            
            % make wrapper yaml
            WriteYaml(yamlwrappathlocal, wrapyaml);
            % call decode fitter
            url = [p.server '/submit_fit'];
            options = weboptions('RequestMethod', 'post', 'ArrayFormat','json');
            pid = webread(url, 'path_fit_meta', yamlwrappathremote, options);
            obj.decodepid=pid;
            
            %update staus
            logfile=[workingdirlocal '/out.log' ];
            fittingstat=webread([server '/status_processes']);
            while strcmp(fittingstat.fit.(['x' num2str(pid)]),'running')
                pause(2)
                alllines=readlines(logfile);
                line=alllines(end);
                obj.status(line)
                drawnow
                fittingstat=webread([server '/status_processes']);
            end
            % read h5, 
            % %add ROI, 
            % setPer('fitinfo') all training and fitting yaml parameters,
            % adjust LocSaver to save those as well
            % %change pixelsize, 
            % %RI mismatch, 
            % %save as _sml.mat. Use loc saver plugin?
        end
        function addFile(obj,file,setinfo)   %for batch processing?
            if isempty(file)
                fileinfo=obj.getPar('loc_fileinfo_set');
                file=fileinfo.imagefile;
            end
            if setinfo
                if contains(file,'/decode/fits')
                    outfile=strrep(file,'.tif','.h5');
                else
                    [~,fn]=fileparts(file);
                    [~,dir]=fileparts(fileparts(file));
                    decodenetwork=obj.getGlobalSetting('DECODE_network_data');
                    outfile=[decodenetwork filesep 'fits' filesep dir filesep fn '.h5'];
                end
                obj.setGuiParameters(struct('outputpath',outfile))
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
            obj.createGlobalSetting('DECODE_network_data','Plugins','network directory for DECODE training and fitting',struct('Style','dir','String',' '))
            yamldefault=[obj.getPar('SettingsDirectory') filesep 'temp' filesep obj.yamldefault];
            if ~exist(yamldefault,'file')
                yamlold=[obj.getPar('SettingsDirectory') filesep 'cameras' filesep 'DECODE_default.yaml'];
                copyfile(yamlold, yamldefault)
            end
            ypar=ReadYaml(yamldefault);
            obj.setGuiParameters(struct('server',ypar.Connect.remote_workstation));
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
%     if isempty(obj.getSingleGuiParameter('outputpath'))
%         obj.setGuiParameters(struct('outputpath',obj.workingdir))
%     end
end
end


function selectoutput(a,b,obj)
outputp=obj.getSingleGuiParameter('outputpath');
if ~exist(outputp,"dir")
    outputp=obj.getGlobalSetting('DECODE_network_data');
end
dir=uigetdir([outputp filesep]);
if ~isempty(dir)
    obj.setGuiParameters(struct('outputpath',dir))
end
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


pard.servert.object=struct('Style','text','String','Server');
pard.servert.position=[3,1];
pard.servert.Width=1;
pard.server.object=struct('Style','edit','String','riesgroup@pc-ries25:8000');
pard.server.position=[3,2];
pard.server.Width=3;

p(1).value=0; p(1).on={}; p(1).off={'refractive_index_mismatch'};
p(2).value=1; p(2).on={'refractive_index_mismatch'}; p(2).off={};
pard.userefractive_index_mismatch.object=struct('Style','checkbox','String','RI mismatch:','Callback',{{@obj.switchvisible,p}});
pard.userefractive_index_mismatch.position=[4,3.5];
pard.userefractive_index_mismatch.Width=1.5;
pard.userefractive_index_mismatch.Optional=true;

pard.refractive_index_mismatch.object=struct('Style','edit','String','.8','Visible','off');
pard.refractive_index_mismatch.position=[4,4.5];
pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
pard.refractive_index_mismatch.Optional=true;
pard.refractive_index_mismatch.Width=0.5;


p(1).value=0; p(1).on={}; p(1).off={'pixelsizex','pixelsizey'};
p(2).value=1; p(2).on={'pixelsizex','pixelsizey'}; p(2).off={};
pard.overwritePixelsize.object=struct('Style','checkbox','String','New pixelsize X,Y (um):','Callback',{{@obj.switchvisible,p}});
pard.overwritePixelsize.position=[4,1];
pard.overwritePixelsize.Width=1.5;
pard.overwritePixelsize.Optional=true;

pard.pixelsizex.object=struct('Style','edit','String','.1','Visible','off');
pard.pixelsizex.position=[4,2.5];
pard.pixelsizex.Width=0.5;
pard.pixelsizex.Optional=true;

pard.pixelsizey.object=struct('Style','edit','String','.1','Visible','off');
pard.pixelsizey.position=[4,3];
pard.pixelsizey.Width=0.5;
pard.pixelsizey.Optional=true;



% pard.syncParameters={{'cal_3Dfile','cal_3Dfile',{'String'}}};

pard.plugininfo.type='WorkflowModule';
pard.plugininfo.description='';
end