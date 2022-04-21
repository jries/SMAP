classdef learnPSF_invmodeling<interfaces.DialogProcessor
    properties        
        outputfile
        cameraSettings=struct('EMon',1,'emgain',1,'conversion',1,'offset',400,'cam_pixelsize_um',0.1);        
        selectedROI=[0 0 0];
        %define class properties if needed
    end
    methods
        function obj=learnPSF_invmodeling(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=true; %set true, if results are shown by default
            obj.guiselector.show=false; %if true, the selector for simple vs complex GUI is shown.
        end       
        function initGui(obj)
            obj.createGlobalSetting('PSFlearning_env','Python','The anaconda environmet path for the PSF learning (eg. /psf_env/):',struct('Style','dir','String','')) 
        end
        function out=run(obj,p)
            envpath=[obj.getGlobalSetting('PSFlearning_env') ];
            if ~exist(envpath,"dir")
                warndlg('please select the PSF learning conda environment in the parameters menu.')
                return
            end
            
              
%             envpath = '/Users/jonasries/opt/anaconda3/envs/psfenv'; %to preferences. Also use different name for env.
            runpath = [fileparts(pwd) '/psfmodelling/examples'];
            paramtemplate=[runpath filesep 'params_default.json'];
            fid = fopen(paramtemplate); 
            raw = fread(fid,inf); 
            str = char(raw'); 
            fclose(fid); 
            pf = jsondecode(str);

            %overwrite with updated parameters:
            pf.plotall=false;
      
            pf.filelist=p.filelist;
            fn1=p.filelist{1};
            [~,~,pf.format]=fileparts(fn1);
            %loss
            loss.mse1=p.loss1(1);
            loss.mse2=p.loss1(2);
            loss.smooth=p.loss1(3);
            loss.edge=p.loss1(4)* 0.01;
            loss.psf_min=p.loss2(1);
            loss.bg_min=p.loss2(2);
            loss.photon_min=p.loss2(3)*1e-6;
            loss.Inorm=p.loss2(4);

            % modality
            switch p.representation.selection
                case 'Voxels'
                    PSFtype='voxel';
                case 'Pupil'
                    PSFtype='pupil';
                    loss.smooth=0;
                case 'Zernike'
                    PSFtype = 'zernike_vector'; 
                    loss.smooth=0;
            end
            roisize=[1 1]*p.roisize;
            pf.gaus_sigma = [2,2];
            pf.max_kernel = [3,3];

            switch p.modality.selection
                case '1 Ch'
                    pf.PSFtype=PSFtype;
                    pf.channeltype='single';
                case '2 Ch'
%                     pf.PSFtype='voxel';
                    pf.channeltype='multi';
                case '4 Pi'
                    pf.PSFtype='voxel';
                    pf.channeltype='4pi';  
                case 'LLS'
                    pf.PSFtype='voxel';
                    pf.channeltype='single'; 
                    rz=round(2700/p.dz); %XXX change?
                    roisize=[rz,roisize];
                    pf.gaus_sigma = [6,2,2];
                    pf.max_kernel = [9,3,3];
            end

            % read meta data
            r=imageloaderAll(fn1,[],obj.P);
            md=r.metadata;
            pf.gain = md.conversion;
            pf.ccd_offset=md.offset;
            pf.pixelsize_x=md.cam_pixelsize_um(1);
            pf.pixelsize_y=md.cam_pixelsize_um(end);

            pf.pixelsize_z=p.dz/1000;
            pf.datapath=[fileparts(fn1) filesep];
            pf.bead_radius=p.beadsize/1000;
            pf.estimate_drift=p.estimate_drift==1;
            pf.vary_photon=p.vary_photon==1;
            pf.usecuda=p.usecuda==1;
            pf.iteration=p.iteration;

            %ROI goes here
            FOV = struct('x_center',obj.selectedROI(1),'y_center',obj.selectedROI(2),'radius',obj.selectedROI(3),'z_start',p.skipframes(1),'z_end',-p.skipframes(end),'z_step',1); % if width and height are zero, use the full FOV, 'z_start' is counting as 0,1,2..., 'z_end' is counting as 0,-1,-2...
            pf.FOV = FOV;

            pf.roi_size=roisize; 
            pf.peak_height=p.segcutoff;


            pf.loss_weight=loss;
            [pfad,fnh]=fileparts(fn1);
            paramfile=fullfile(pfad,[fnh '_par.json']);
            pf.savename=[pfad filesep 'psfmodel_' fnh];

            encode_str = jsonencode(pf,'PrettyPrint',true);
            fid = fopen(paramfile,'w'); 
            fwrite(fid,encode_str); 
            fclose(fid);

            %% run python script
            [p1,env]=fileparts(envpath);
%             condapath=fileparts(p1);
            pythonfile = 'learn_psf.py';
            command = ['python ' pythonfile ' "' paramfile '"'];
            currentpath=pwd;

            logfile=strrep(paramfile,'.json','_log.txt');

            [pid,status, results]=systemcallpython(envpath,command,runpath,logfile);
             disp(pid)
            cd(currentpath);
            out=[]; %no output
            % use PID to see if still running

            t=timer('StartDelay',1,'Period',1,'TasksToExecute',100,'ExecutionMode','fixedDelay');
            pt.statusHandle=obj.P.par.mainGui.content.guihandles.status;
            pt.PID=pid;
            pt.savefile=pf.savename;
            t.TimerFcn={@displayprogress_timer,logfile,pt};
            t.StartFcn={@displayprogress_timer,'',pt};
            t.StopFcn={@fitting_done,obj,pt};
            t.start;       

        end

        
        function pard=guidef(obj)
            pard.load_filest.object=struct('String','beads:','Style','text');
            pard.load_filest.position=[1,1];
            pard.load_filest.Width=0.5;
            pard.load_files.object=struct('String','load','Style','pushbutton','Callback',{{@load_files_callback,obj}});
            pard.load_files.position=[1,1.5];
            pard.load_files.Width=0.5;
            pard.filelist.object=struct('String','','Style','edit','Max',10);
            pard.filelist.position=[1,2.];            
            pard.filelist.Width=2;

            pard.camparbutton.object=struct('String','Cam par','Style','pushbutton','Callback',{{@campar_callback,obj}});
            pard.camparbutton.position=[1,4];
            pard.camparbutton.Width=0.5;
            pard.lockcampar.object=struct('String','lock','Style','checkbox');
            pard.lockcampar.position=[1,4.5];
            pard.lockcampar.Width=0.5;

            pard.dzt.object=struct('String','dz (nm)','Style','text');
            pard.dzt.position=[2,3];
            pard.dzt.Width=0.5;
            pard.dz.object=struct('String','50','Style','edit');
            pard.dz.position=[2,3.5];
            pard.dz.Width=0.5;

            pard.modalityt.object=struct('String','Modality:','Style','text');
            pard.modalityt.position=[2,1];            
            pard.modality.object=struct('String',{{'1 Ch','2 Ch','4 Pi','LLS'}},'Style','popupmenu','Tag','modality','Callback',{{@modechanged,obj}});
            pard.modality.position=[2,1.5];  
            pard.modality.Width=0.75;

            pard.representation.object=struct('String',{{'Voxels','Pupil','Zernike'}},'Style','popupmenu','Tag','representation','Callback',{{@modechanged,obj}});
            pard.representation.position=[2,2.25]; 
            pard.representation.Width=0.75;

            pard.roisizet.object=struct('String','roi (pix)','Style','text');
            pard.roisizet.position=[2,4];
            pard.roisizet.Width=0.5;
            pard.roisize.object=struct('String','21','Style','edit');
            pard.roisize.position=[2,4.5];
            pard.roisize.Width=0.5;
            
            lw=3;
            %2Ch
            pard.mirrortypet.object=struct('String','mirror','Style','text','Visible','off');
            pard.mirrortypet.position=[lw,3]; 
            pard.mirrortypet.Width=.5;
            pard.mirrortype.object=struct('String',{{'none','up-down','right-left'}},'Style','popupmenu','Visible','off');
            pard.mirrortype.position=[lw,3.5]; 
            pard.mirrortype.Width=1;
            pard.channelarranget.object=struct('String','channel','Style','text','Visible','off');
            pard.channelarranget.position=[lw,1]; 
            pard.channelarranget.Width=.5;
            pard.channelarrange.object=struct('String',{{'up-down','right-left'}},'Style','popupmenu','Visible','off');
            pard.channelarrange.position=[lw,1.5]; 
            pard.channelarrange.Width=1;
            %4Pi
            pard.zTt.object=struct('String','Period (µm)','Style','text','Visible','off');
            pard.zTt.position=[lw,1]; 
            pard.zTt.Width=1;
            pard.zT.object=struct('String','0.26','Style','edit','Visible','off');
            pard.zT.position=[lw,2]; 
            pard.zT.Width=0.5;

            lw=4;
            pard.segmentationt.object=struct('String','Segmenation: cutoff','Style','text');
            pard.segmentationt.position=[lw,1];  
            pard.segmentationt.Width=1.5;
            pard.segcutoff.object=struct('String','0.2','Style','edit');
            pard.segcutoff.position=[lw,2.]; 
            pard.segcutoff.Width=0.5;

            pard.skipframest.object=struct('String','skip frames','Style','text');
            pard.skipframest.position=[lw,3]; 
            pard.skipframest.Width=0.75;
            pard.skipframes.object=struct('String','0 0','Style','edit');
            pard.skipframes.position=[lw,3.75]; 
            pard.skipframes.Width=0.5;
            pard.selectroi.object=struct('String','Select ROI','Style','pushbutton','Callback',{{@selectroi_callback,obj}});
            pard.selectroi.position=[lw,4.25];     
            pard.selectroi.Width=0.75;
            
            lw=5;
            pard.estimate_drift.object=struct('String','est drift','Style','checkbox');
            pard.estimate_drift.position=[lw,1]; 
            pard.vary_photon.object=struct('String','vary N','Style','checkbox');
            pard.vary_photon.position=[lw,2]; 
            pard.usecuda.object=struct('String','use cuda','Style','checkbox');
            pard.usecuda.position=[lw,3]; 
            pard.iterationt.object=struct('String','iterations','Style','text');
            pard.iterationt.position=[lw,4]; 
            pard.iterationt.Width=0.5;
            pard.iteration.object=struct('String','100','Style','edit');
            pard.iteration.position=[lw,4.5]; 
            pard.iteration.Width=0.5;

            lw=6;
            pard.beadsizet.object=struct('String','Bead radius nm','Style','text');
            pard.beadsizet.position=[lw,1];  
            pard.beadsizet.Width=1.5;
            pard.beadsize.object=struct('String','0','Style','edit');
            pard.beadsize.position=[lw,2];  
            pard.beadsize.Width=0.5;

            lw=7;
            pard.lmse1t.object=struct('String','Loss: mse1, mse2, smooth, edge','Style','text');
            pard.lmse1t.position=[lw,1];  
            pard.lmse1t.Width=1.75;
            pard.loss1.object=struct('String','1 1 1 1','Style','edit');
            pard.loss1.position=[lw,2.75];  
            pard.loss1.Width=0.5;

            pard.loss2t.object=struct('String','min: psf, bg, phot, Inorm','Style','text');
            pard.loss2t.position=[lw,3.25];  
            pard.loss2t.Width=1.25;
            pard.loss2.object=struct('String','1 1 1 0','Style','edit');
            pard.loss2.position=[lw,4.5];  
            pard.loss2.Width=0.5;

            pard.plugininfo.type='ProcessorPlugin'; %type of plugin. Currently: ProcessorPlugin, WorkflowModule, WorkflowFitter, Renderer, LoaderPlugin, SaverPlugin, ROI_Analyze, ROI_Evaluate,WorkflowIntensity
        end
    end
end

function fitting_done(a,b,obj,pt)
fn=dir([pt.savefile '*.h5']);
if isempty(fn)
    disp('no output file written')
    return
end
[~,ind]=max([fn(:).datenum]);
% finf=h5info([fn(ind).folder filesep fn(ind).name]);
v=loadh5([fn(ind).folder filesep fn(ind).name]);

axb=obj.initaxis('beads');
xb=double(v.rois.cor(:,1));
yb=double(v.rois.cor(:,2));

plot(axb,xb,yb,'o')
axis(axb,'equal')
nbeads=size(v.rois.psf_data,1);
title([num2str(nbeads) ' beads found'])


zfit=reshape(v.locres.P(5,:),[],nbeads);
frames=(1:size(zfit,1))';
ax=obj.initaxis('zfit');
plot(ax,frames,frames,'k');
hold(ax,'on');
plot(ax,frames,zfit)
hold(ax,'off')
xlabel('frame')
ylabel('frame fit')
end

function modechanged(a,b,obj)
if strcmpi(a.Tag,'modality')
    p(1).value=1;p(1).off={'mirrortypet','mirrortype','channelarranget','channelarrange','zTt','zT'};p(1).on={};
    p(2).value=2;p(2).off={'zTt','zT'};p(2).on={'mirrortypet','mirrortype','channelarranget','channelarrange'};
    p(3).value=3;p(3).off={'mirrortypet','mirrortype','channelarranget','channelarrange'};p(3).on={'zTt','zT'};
    obj.switchvisible(a,b,p);
elseif strcmpi(a.Tag,'representation')
    l1=obj.getSingleGuiParameter('loss1');
    switch obj.getSingleGuiParameter('representation').selection
        case 'Voxels'
            l1(3)=1;
        case {'Zernike','Pupil'}
            l1(3)=0;
    end
    obj.setGuiParameters(struct('loss1',l1));

end
    %adjust loss etc and other specific settings
end

function load_files_callback(a,b,obj)
sf=selectManyFiles;
sf.guihandles.filelist.String=(obj.guihandles.filelist.String);
waitfor(sf.handle);
obj.guihandles.filelist.String=sf.filelist;
obj.guihandles.filelist.Value=1;
r=imageloaderAll(sf.filelist{1},[],obj.P);
md=r.metadata;
fn=fieldnames(obj.cameraSettings);
islock=obj.getSingleGuiParameter('lockcampar');
if ~islock
    obj.cameraSettings=copyfields(obj.cameraSettings,md,fn);
end
obj.selectedROI=[0 0 0];
end

function selectroi_callback(a,b,obj)
fl=obj.getSingleGuiParameter('filelist');
r=imageloaderAll(fl{1},[],obj.P);
img=r.getmanyimages([],'mat');
f=figure;
ax=gca;
imagesc(ax,mean(img,3))
title('Select ROI, double click when done')
axis equal
roi=drawcircle(ax);
obj.selectedROI=customWait(roi);
delete(f)
end
function pos = customWait(hROI)
l = addlistener(hROI,'ROIClicked',@clickCallback);
uiwait;
delete(l);
pos = [hROI.Position,hROI.Radius];
end
function clickCallback(~,evt)
if strcmp(evt.SelectionType,'double')
    uiresume;
end
end

function campar_callback(a,b,obj)
fn=fieldnames(obj.cameraSettings);
for k=length(fn):-1:1
    fields{k}=fn{k};
    defAns{k}=converttostring(obj.cameraSettings.(fn{k}));
end
answer=inputdlg(fields,'Acquisition settings',1,defAns);
if ~isempty(answer) && ~obj.getSingleGuiParameter('lockcampar')
    for k=1:length(fn)
        if isnumeric(obj.cameraSettings.(fn{k}))||islogical(obj.cameraSettings.(fn{k}))
            obj.cameraSettings.(fn{k})=str2num(answer{k});
        else
            obj.cameraSettings.(fn{k})=(answer{k});
        end
    end
end
if obj.getSingleGuiParameter('lockcampar')
    warning('cannot update camera paramters because they are locked')
end
end

function out=converttostring(in)
if iscell(in)
    out=join(in,',');
    out=out{1};
elseif ischar(in)
    out=in;
else
    out=num2str(in);
end
end

function displayprogress_timer(obj,event,logfile,pt)
if isempty(logfile)
    obj.UserData.starttime=now;
    obj.UserData.updatetime=datetime;
    obj.UserData.oldtextlength=0;
    pt.statusHandle.String='timer init';drawnow
    return
end

if exist(logfile,'file') && dir(logfile).datenum>obj.UserData.starttime
    alllines=readlines(logfile,'WhitespaceRule','trim','EmptyLineRule','skip');
    if isempty(alllines)
        return
    end
    line=alllines(end);
    if ~isempty(line) && length(alllines)>obj.UserData.oldtextlength
        pt.statusHandle.String=(line);
        obj.UserData.updatetime=datetime;
        obj.UserData.oldtextlength=length(alllines);
        drawnow
    end     
end
timeout=(datetime-obj.UserData.updatetime)>duration(0,0,60);
fn=dir([pt.savefile '*.h5']);
filesaved=false;
if ~isempty(fn)
    for k=1:length(fn)
        if fn(k).datenum>obj.UserData.starttime && datetime(fn(k).datenum,'ConvertFrom','datenum')+duration(0,0,10)<datetime('now')
            filesaved=true;
        end
    end
end

piddeleted=~processstatus(pt.PID);

if timeout || filesaved || piddeleted
    disp(['timer ended. timeout' num2str(timeout) 'filesaved' num2str(filesaved) 'piddeleted' num2str(piddeleted)])
    pt.statusHandle.String='timer done';drawnow
    obj.stop
end
end


