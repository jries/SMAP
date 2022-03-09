classdef learnPSF_invmodeling<interfaces.DialogProcessor
    % PLUGIN_TEMPLATE Summary of this plugin goes here
    % put a description of your plugin here.
        %replace Plugin_Template by filename   
    properties
        outputfile
        %define class properties if needed
    end
    methods
        function obj=learnPSF_invmodeling(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=true; %set true, if results are shown by default
            obj.guiselector.show=true; %if true, the selector for simple vs complex GUI is shown.
        end       
        function initGui(obj)
            %is called after the GUI, defined in guidef, is made. Here you
            %can add additional GUI components, e.g. some that are not
            %defined by simple uicontrols (e.g. uitable), set additional
            %synchronizations and callbacks.
        end
        function out=run(obj,p)
            paramtemplate=[fileparts(pwd) filesep 'psfmodelling' filesep 'examples' filesep 'param_single.json'];
            
            fid = fopen(paramtemplate); 
            raw = fread(fid,inf); 
            str = char(raw'); 
            fclose(fid); 
            pf = jsondecode(str);

            %overwrite with updated parameters:
            fn1=p.filelist{1};
            
            % in future: read meta data
%             r=imageloaderAll(fn1);

            pf.pixelsize_z=p.dz;
            pf.datapath=fileparts(fn1);
            pf.bead_radius=p.beadsize;
            pf.estimate_drift=p.estimate_drift;
            pf.vary_photon=p.vary_photon;
            pf.usecuda=p.usecuda;
            pf.iteration=p.iteration;
            loss.mse1=p.lmse1;
            loss.mse2=p.lmse2;
            loss.smooth=p.lsmooth;
            loss.edge=p.ledge * 0.01;
            loss.psf_min=p.lpsfmin;
            loss.bg_min=p.lbgmin;
            loss.photon_min=p.lphotonmin *1e-6;
            loss.Inorm=p.lInorm;
            pf.loss_weight=loss;
            [pfad,fnh]=fileparts(fn1);
            paramfile=fullfile(pfad,[fnh '_par.json']);

            encode_str = jsonencode(pf,'PrettyPrint',true);
            fid = fopen(paramfile,'w'); 
            fwrite(fid,encode_str); 
            fclose(fid);
            %% run python script
            envpath = '/Users/ries/opt/miniconda3/envs/myenv';
            runpath = '/Users/ries/git/psfmodelling/examples';
            %envpath = 'C:\Users\Ries Lab\anaconda3\envs\psfinv_tf2';
            %runpath = 'C:\Users\Ries Lab\git\PSFlearningTF2\psfmodelling\examples';
            [p1,env]=fileparts(envpath);
            condapath=fileparts(p1);
            pythonfile = 'learn_singleChannel.py';
            command = ['python ' pythonfile ' ' paramfile ' &'];
            currentpath=pwd;
            if ispc
                pcall=['call "' condapath '\Scripts\activate.bat" ' env ' & cd "' runpath '" & ' command ' & exit &'];
                
            else
                cd(runpath)
                pcall=[envpath '/bin/'  command ];
            end
            
            [obj.decodepid,status, results]=systemcallpython(pdecode,command,decodepath,logfile);
            [status, results]=system(pcall,'-echo');
            cd(currentpath);
           
            out=[]; %no output
           
        end
        function pard=guidef(obj)
           %pard structure can be =[]; All fields are optional.
            
            %define your GUI: for every GUI control define a structure as
            %follows, replace 'guiobject' by a name by which you want to
            %access the parameters
            % required fields: .object=struct(...). Defines a uicontrol
            % object. Just pass on all arguments you would otherwise pass
            % on to uicontrol. Careful: if you want to pass on a cell
            % array, you need to put it into double brackets: {{'a','b'}}.
            % Use i.e. for 'String' property of popupmenus or for
            % 'Callback' with options.
            %please note, that if you define your own callback function,
            %you need to take care of synchronizaiton with outputParameters
            %yourself, e.g. by calling obj.obj.updateGuiParameter(0,0,guihandle);
            %.position: relative position in a 4x11 grid. You can use
            %non-integer values
            %optional fields: Width, Height. In relative units.
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

            pard.representationt.object=struct('String','Representation','Style','text');
            pard.representationt.position=[2,1];            
            pard.representation.object=struct('String',{{'voxels','OTF','Zernike'}},'Style','popupmenu');
            pard.representation.position=[2,2];   
            pard.modalityt.object=struct('String','Modality:','Style','text');
            pard.modalityt.position=[2,3];            
            pard.modality.object=struct('String',{{'1 Ch','2 Ch','4 Pi'}},'Style','popupmenu');
            pard.modality.position=[2,4];  

            lw=3;
            pard.segmentationt.object=struct('String','Segmenation:','Style','text');
            pard.segmentationt.position=[lw,1];  
            
            lw=4;
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

            lw=5;
            pard.beadsizet.object=struct('String','Bead radius nm','Style','text');
            pard.beadsizet.position=[lw,1];  
            pard.beadsizet.Width=1.5;
            pard.beadsize.object=struct('String','0','Style','edit');
            pard.beadsize.position=[lw,2];  
            pard.beadsize.Width=0.5;
            pard.selectroi.object=struct('String','Select ROI','Style','pushbutton','Callback',{{@selectroi_callback,obj}});
            pard.selectroi.position=[lw,4];              

            lw=6;
            pard.lmse1t.object=struct('String','mse1','Style','text');
            pard.lmse1t.position=[lw,1];  
            pard.lmse1t.Width=0.5;
            pard.lmse2t.object=struct('String','mse2','Style','text');
            pard.lmse2t.position=[lw,1.5];  
            pard.lmse2t.Width=0.5;
            pard.lsmootht.object=struct('String','smooth','Style','text');
            pard.lsmootht.position=[lw,2];  
            pard.lsmootht.Width=0.5;
            pard.ledget.object=struct('String','edge','Style','text');
            pard.ledget.position=[lw,2.5];  
            pard.ledget.Width=0.5;
            pard.lpsfmint.object=struct('String','psf min','Style','text');
            pard.lpsfmint.position=[lw,3];  
            pard.lpsfmint.Width=0.5;
            pard.lbgmint.object=struct('String','bg min','Style','text');
            pard.lbgmint.position=[lw,3.5];  
            pard.lbgmint.Width=0.5;
            pard.lphotonmint.object=struct('String','phot min','Style','text');
            pard.lphotonmint.position=[lw,4.];  
            pard.lphotonmint.Width=0.5;
            pard.lInormt.object=struct('String','I norm','Style','text');
            pard.lInormt.position=[lw,4.5];  
            pard.lInormt.Width=0.5;
            lw=lw+1;
            pard.lmse1.object=struct('String','1','Style','edit');
            pard.lmse1.position=[lw,1];  
            pard.lmse1.Width=0.5;
            pard.lmse2.object=struct('String','1','Style','edit');
            pard.lmse2.position=[lw,1.5];  
            pard.lmse2.Width=0.5;
            pard.lsmooth.object=struct('String','1','Style','edit');
            pard.lsmooth.position=[lw,2];  
            pard.lsmooth.Width=0.5;
            pard.ledge.object=struct('String','1','Style','edit');
            pard.ledge.position=[lw,2.5];  
            pard.ledge.Width=0.5;
            pard.lpsfmin.object=struct('String','1','Style','edit');
            pard.lpsfmin.position=[lw,3];  
            pard.lpsfmin.Width=0.5;
            pard.lbgmin.object=struct('String','1','Style','edit');
            pard.lbgmin.position=[lw,3.5];  
            pard.lbgmin.Width=0.5;
            pard.lphotonmin.object=struct('String','1','Style','edit');
            pard.lphotonmin.position=[lw,4.];  
            pard.lphotonmin.Width=0.5;
            pard.lInorm.object=struct('String','1','Style','edit');
            pard.lInorm.position=[lw,4.5];  
            pard.lInorm.Width=0.5;
            
          
            %automatically hide and show other gui parameters based on the
            %value of a specific parameter
            %pass on structure array, each element with the fields
             %value choose which Value the description corresponds to
             % on: array of names of gui objects to switched on when this
             % value is selected
%              % off: array of names of gui objects to switched off
%             p(1).value=0; p(1).on={}; p(1).off={'guiobject2','guiobject'};
%             p(2).value=1; p(2).on={'guiobject2','guiobject'}; p(2).off={};
% 
%             pard.onofftoggle.object=struct('Style','checkbox','String','show','Callback',{{@obj.switchvisible,p}});
            pard.plugininfo.type='ProcessorPlugin'; %type of plugin. Currently: ProcessorPlugin, WorkflowModule, WorkflowFitter, Renderer, LoaderPlugin, SaverPlugin, ROI_Analyze, ROI_Evaluate,WorkflowIntensity
        end
    end
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

end
function selectroi_callback(a,b,obj)
end



