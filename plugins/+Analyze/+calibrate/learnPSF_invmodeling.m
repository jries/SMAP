classdef learnPSF_invmodeling<interfaces.DialogProcessor
    properties
        outputfile
        %define class properties if needed
    end
    methods
        function obj=learnPSF_invmodeling(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=false; %set true, if results are shown by default
            obj.guiselector.show=false; %if true, the selector for simple vs complex GUI is shown.
        end       
        function initGui(obj)
        end
        function out=run(obj,p)
            envpath = '/Users/jonasries/opt/anaconda3/envs/myenv'; %to preferences. Also use different name for env.
            runpath = [fileparts(pwd) '/psfmodelling/examples'];
            paramtemplate=[runpath filesep 'params.json'];
            fid = fopen(paramtemplate); 
            raw = fread(fid,inf); 
            str = char(raw'); 
            fclose(fid); 
            pf = jsondecode(str);

            %overwrite with updated parameters:
            fn1=p.filelist{1};

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
                    PSFtype='voxels';
                case 'Pupil'
                    PSFtype='pupil';
                    loss.smooth=0;
                case 'Zernike'
                    PSFtype = 'zernike_vector'; 
                    loss.smooth=0;
            end

            switch p.modality.selection
                case '1 Ch'
                    pf.PSFtype=PSFtype;
                    pf.channeltype='single';
                case '2 Ch'
                    pf.PSFtype='voxels';
                    pf.channeltype='multi';
                case '4 Pi'
                    pf.PSFtype='voxels';
                    pf.channeltype='4pi';          
            end


            
            % in future: read meta data
%             r=imageloaderAll(fn1);

            pf.pixelsize_z=p.dz;
            pf.datapath=fileparts(fn1);
            pf.bead_radius=p.beadsize;
            pf.estimate_drift=p.estimate_drift;
            pf.vary_photon=p.vary_photon;
            pf.usecuda=p.usecuda;
            pf.iteration=p.iteration;

            pf.loss_weight=loss;
            [pfad,fnh]=fileparts(fn1);
            paramfile=fullfile(pfad,[fnh '_par.json']);
            pf.savename=[pfad 'psfmodel_' fnh];

            encode_str = jsonencode(pf,'PrettyPrint',true);
            fid = fopen(paramfile,'w'); 
            fwrite(fid,encode_str); 
            fclose(fid);
            %% run python script


            [p1,env]=fileparts(envpath);
            condapath=fileparts(p1);
            pythonfile = 'learn_singleChannel.py';
            command = ['python ' pythonfile ' "' paramfile '"'];
            currentpath=pwd;

            logfile=strrep(paramfile,'.json','_log.txt');

            [pid,status, results]=systemcallpython(envpath,command,runpath,logfile);

            cd(currentpath);
            out=[]; %no output

            t=timer('StartDelay',1,'Period',1,'TasksToExecute',100,'ExecutionMode','fixedDelay');
            t.TimerFcn={@displayprogress_timer,logfile,obj.P.par.mainGui.content.guihandles.status};
            t.StartFcn={@displayprogress_timer,'',@obj.status};
            t.start;       
        end

        
        function pard=guidef(obj)
           %pard structure can be =[]; All fields are optional.

            pard.load_filest.object=struct('String','bead stacks:','Style','text');
            pard.load_filest.position=[1,1];
            pard.load_files.object=struct('String','load','Style','pushbutton','Callback',{{@load_files_callback,obj}});
            pard.load_files.position=[1,2];
            pard.filelist.object=struct('String','','Style','edit','Max',10);
            pard.filelist.position=[1,3];            
            pard.filelist.Width=1;
            pard.dzt.object=struct('String','dz (nm)','Style','text');
            pard.dzt.position=[1,4];
            pard.dzt.Width=0.5;
            pard.dz.object=struct('String','50','Style','edit');
            pard.dz.position=[1,4.5];
            pard.dz.Width=0.5;

%             pard.representationt.object=struct('String','Representation','Style','text');
%             pard.representationt.position=[2,3];            
              


            pard.modalityt.object=struct('String','Modality:','Style','text');
            pard.modalityt.position=[2,1];            
            pard.modality.object=struct('String',{{'1 Ch','2 Ch','4 Pi'}},'Style','popupmenu','Callback',{{@modechanged,obj}});
            pard.modality.position=[2,1.5];  
            pard.modality.Width=0.75;

            pard.representation.object=struct('String',{{'Voxels','Pupil','Zernike'}},'Style','popupmenu');
            pard.representation.position=[2,2.25]; 
            pard.representation.Width=0.75;
            
            lw=3;
            %2Ch
            pard.mirrortypet.object=struct('String','mirror','Style','text','Visible','off');
            pard.mirrortypet.position=[lw,3]; 
            pard.mirrortypet.Width=1;
            pard.mirrortype.object=struct('String',{{'none','up-down','right-left'}},'Style','popupmenu','Visible','off');
            pard.mirrortype.position=[lw,4]; 
            pard.mirrortype.Width=1;
            pard.channelarranget.object=struct('String','channel','Style','text','Visible','off');
            pard.channelarranget.position=[lw,1]; 
            pard.channelarranget.Width=1;
            pard.channelarrange.object=struct('String',{{'up-down','right-left'}},'Style','popupmenu','Visible','off');
            pard.channelarrange.position=[lw,2]; 
            pard.channelarrange.Width=1;
            %4Pi
            pard.zTt.object=struct('String','Period (µm)','Style','text','Visible','off');
            pard.zTt.position=[lw,1]; 
            pard.zTt.Width=1;
            pard.zT.object=struct('String','0.26','Style','edit','Visible','off');
            pard.zT.position=[lw,2]; 
            pard.zT.Width=1;

            lw=4;
            pard.segmentationt.object=struct('String','Segmenation:','Style','text');
            pard.segmentationt.position=[lw,1];  
            
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
            pard.selectroi.object=struct('String','Select ROI','Style','pushbutton','Callback',{{@selectroi_callback,obj}});
            pard.selectroi.position=[lw,4];              

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

%             pard.onofftoggle.object=struct('Style','checkbox','String','show','Callback',{{@obj.switchvisible,p}});
            pard.plugininfo.type='ProcessorPlugin'; %type of plugin. Currently: ProcessorPlugin, WorkflowModule, WorkflowFitter, Renderer, LoaderPlugin, SaverPlugin, ROI_Analyze, ROI_Evaluate,WorkflowIntensity
        end
    end
end

function modechanged(a,b,obj)
    p(1).value=1;p(1).off={'mirrortypet','mirrortype','channelarranget','channelarrange','zTt','zT'};p(1).on={};
    p(2).value=2;p(2).off={'zTt','zT'};p(2).on={'mirrortypet','mirrortype','channelarranget','channelarrange'};
    p(3).value=3;p(3).off={'mirrortypet','mirrortype','channelarranget','channelarrange'};p(3).on={'zTt','zT'};
   obj.switchvisible(a,b,p);
    
end

function load_files_callback(a,b,obj)
sf=selectManyFiles;
sf.guihandles.filelist.String=(obj.guihandles.filelist.String);
waitfor(sf.handle);
obj.guihandles.filelist.String=sf.filelist;
obj.guihandles.filelist.Value=1;

    ind=strfind(sf.filelist{1},';');
    if ~isempty(ind)
        fileh=sf.filelist{1}(1:ind-1);
    else
        fileh=sf.filelist{1};
    end
    [path,file]=fileparts(fileh);
    if length(sf.filelist)>1
        fileh2=sf.filelist{2};
        path2=fileparts(fileh2);
        if ~strcmp(path,path2) %not the same: look two hierarchies down
            if strcmp(fileparts(path),fileparts(path2))
                path=fileparts(path);
            elseif strcmp(fileparts(fileparts(path)),fileparts(fileparts(path2)))
                path=fileparts(fileparts(path));
            end
        end
    end
    obj.outputfile.String=[path filesep file '_3dcal.mat'];
    try
        r=imageloaderAll(fileh,[],obj.P);
        mirror=r.metadata.EMon;
%                     obj.guihandles.emgain.Value=mirror;
    catch err
        disp('EM mirror could not be defined automatically, set manually')
    end
    % if dz found: add to GUI

end
function selectroi_callback(a,b,obj)
end


function displayprogress_timer(obj,event,logfile,handle)
if isempty(logfile)
    obj.UserData.starttime=now;
    obj.UserData.updatetime=datetime;
    obj.UserData.oldtextlength=0;
    handle.String='timer init';drawnow
    return
end

if exist(logfile,'file') && dir(logfile).datenum>obj.UserData.starttime
    alllines=readlines(logfile,'WhitespaceRule','trim','EmptyLineRule','skip');
    if isempty(alllines)
        return
    end
    line=alllines(end);
    if ~isempty(line) && length(alllines)>obj.UserData.oldtextlength
        handle.String=(line);
        obj.UserData.updatetime=datetime;
        obj.UserData.oldtextlength=length(alllines);
        drawnow
    end     
end
if datetime-obj.UserData.updatetime>duration(0,0,20)
    disp('timer timeout')
    handle.String='timer done';drawnow
    obj.stop
end
end


