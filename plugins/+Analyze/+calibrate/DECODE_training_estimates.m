classdef DECODE_training_estimates<interfaces.DialogProcessor
%     Saves a training file for DECODE
    properties
        yamldefault='DECODE_local.yaml';
        yamlfile
        yamlpar
        jsontypes
        decodeprocess
        tensorboardprocess
        decodepid
        fitlocal
    end
    methods
        function obj=DECODE_training_estimates(varargin)  
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=false;
        end
        function out=run(obj,p)
           out=[];       
           if p.trainlocal
               obj.fitlocal=true;
               pdecode=strrep([obj.getGlobalSetting('DECODE_path') ],'\','/');
               if ~exist(pdecode,'dir')
                   warning('Decode not found, please specify in the SMAP/Preferences menu in the Plugin tab.');
                   return
               end

               [~,fname]=fileparts(obj.locData.files.file(1).name);
               yamlpath=[obj.yamlpar.InOut.experiment_out filesep 'DECODE_train_local_' fname '.yaml'];
               if ~exist(obj.yamlpar.InOut.experiment_out,'dir')
                   mkdir(obj.yamlpar.InOut.experiment_out);
               end
               finalizejson(obj);
               saveyaml(obj.yamlpar,yamlpath);
               logdir=[fileparts(obj.yamlpar.InOut.experiment_out) filesep 'log'];

               command=['python -m decode.neuralfitter.train.train -p ' yamlpath ' -l ' logdir];
               decodepath=[fileparts(pwd) filesep 'DECODE'];
               [pid,status, results]=systemcallpython(pdecode,command,decodepath);
               port=num2str(obj.yamlpar.SMAP.tensorboardport);
               tfcommand=['tensorboard --samples_per_plugin images=100 --port=' port ' --logdir=' logdir];
               getpids('tensorboard.exe','kill');
               [pidtb,status, results]=systemcallpython(pdecode,tfcommand,decodepath);
               obj.decodeprocess=pid;
               obj.tensorboardprocess.pid=pidtb;
               obj.tensorboardprocess.port=port;

           else %workstation via HTTP
               obj.fitlocal=false;
               %make output directory
%                outdir=[obj.yamlpar.Connect.local_network_storage  'training' filesep obj.yamlpar.InOut.experiment_out];
               local_network_storage=[obj.getGlobalSetting('DECODE_network_data') '/'];
               outdirlocal=[local_network_storage 'training' '/' obj.yamlpar.InOut.experiment_out ]; %later initialize with decode/experiment path. set when settign 3dcal.
               if ~exist(outdirlocal,'dir')
                   mkdir(outdirlocal);
               end
               %copy files
               finalizejson(obj);
               yamlp=obj.yamlpar;


               calfile=obj.yamlpar.InOut.calibration_file;
               [pc,fout,fext]=fileparts(calfile);
               if ~exist([outdirlocal  fout fext],'file')
                   copyfile(calfile,outdirlocal);
               end

               %change yaml parameters
               server=obj.getGlobalSetting('DECODE_server');
               status=webread([server '/status']);
               status.watch_dir=[status.watch_dir filesep];
         
               outpath=['training/' obj.yamlpar.InOut.experiment_out]; 
               yamlp.InOut.calibration_file=strrep(fullfile(outpath,[fout fext]),'\','/');
               %save yaml
               [~,fname]=fileparts(obj.locData.files.file(1).name);
               yamlname=['DECODE_train_' fname '.yaml'];
               yamlpath=[outdirlocal yamlname];
               yamlpath_remote=[outpath yamlname];
               % gpu
               gpustat=webread([server '/status_gpus']);
               [gpus,gpurec]=parsegpustathttp(gpustat);
               yamlp.Hardware.device=gpurec;
               yamlp.Hardware.device_simulation=gpurec;

               yamlp.InOut.experiment_out=strrep(fullfile(outpath),'\','/');
               saveyaml(yamlp,yamlpath);

               url=[server '/submit_training'];
               options = weboptions('RequestMethod', 'post', 'ArrayFormat','json');
               pid = webread(url, 'path_param', yamlpath_remote, options);
               l=length(obj.decodepid);
               obj.decodepid(l+1).pid=pid;
               obj.decodepid(l+1).date=now;

%                http://pc-ries25:8000/docs#/default/status_status_get
%               nvidia-smi
%                 htop
           end

           obj.guihandles.stoplearning.Visible='on';
           obj.guihandles.tensorboard.Visible='on';
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.DialogProcessor(obj);
            obj.createGlobalSetting('DECODE_path','DECODE','The anaconda environmet path of decode (eg. /decode_env/):',struct('Style','dir','String','decode')) 
            obj.createGlobalSetting('DECODE_network_data','DECODE','network directory for DECODE training and fitting',struct('Style','dir','String',' '))
            obj.createGlobalSetting('DECODE_server','DECODE','network directory for DECODE training and fitting',struct('Style','edit','String','http://pc-ries25:8000'))
                 
            yamldefault=[obj.getPar('SettingsDirectory') filesep 'temp' filesep obj.yamldefault];
            if ~exist(yamldefault,'file')
                yamlold=[obj.getPar('SettingsDirectory') filesep 'cameras' filesep 'DECODE_default.yaml'];
                copyfile(yamlold, yamldefault)
            end
            obj.yamlfile=yamldefault;
            tt=uitable(obj.handle);
            tt.Position=obj.guihandles.partablepos.Position;
            tt.Position(4)=tt.Position(4)*7;
            obj.guihandles.parttable=tt;
            obj.yamlpar=ReadYamlSimple(yamldefault);
            obj.yamlfile=yamldefault;
            makejsontable(obj);
        end
    end
end


% function shfile=savesh(js,yamlname,expout)
% shfile=fullfile(js.Connect.local_network_storage,expout,'starttraining.sh');
% remoteyaml=fullfile(js.Connect.remote_network_storage, expout, yamlname);
% command=['ssh ' js.Connect.remote_workstation ' ''nohup ' ...
%     js.Connect.remote_decode_path ' -m decode.neuralfitter.train.live_engine -p ' ...
%     remoteyaml ' -l ' fileparts(remoteyaml) filesep 'runs > foo.out 2> foo.err < /dev/null &'' '];
% fid=fopen(shfile,'w');
% fprintf(fid,command);
% fclose(fid);
% 
% end

% function command=savetbsh(js)
% shfile=fullfile(js.Connect.local_network_storage,js.InOut.experiment_out,'starttensorboard.sh');
% pdecode=fileparts(js.Connect.local_decode_path);
% command=[pdecode '/tensorboard --samples_per_plugin images=100 --port=6008 --logdir=' js.Connect.local_network_storage js.InOut.experiment_out  'runs'];
% end

function makejsontable(obj)
use.SMAP={'set_emitters_per_um2','zrange_nm','tensorboardport'};
use.InOut={'calibration_file','experiment_out'};
use.Simulation={'intensity_mu_sig','lifetime_avg','bg_uniform','PSF_min_mod_max'};
use.Camera={'em_gain','e_per_adu','baseline','read_sigma','spur_noise','px_size'};
use.Hardware={'device','device_simulation'};
% use.Connect={'remote_workstation','local_decode_path','local_network_storage'};
fnu=fieldnames(use);
js=obj.yamlpar;
if ~isfield(js, 'SMAP')
    js.SMAP.set_emitters_per_um2='HD';
    js.SMAP.zrange_nm=[-0 0];
end
if ~isfield(js, 'Connect')
    js.Connect.remote_workstation='';
    js.Connect.remote_decode_path='';
    js.Connect.local_decode_path='';
end

ind=1;
for f=1:length(fnu)
    fnh=use.(fnu{f});
    for k=1:length(fnh)
        tab{ind,1}=fnu{f};
        tab{ind,2}=fnh{k};
        if isfield(js.(fnu{f}),fnh{k})
            vh=js.(fnu{f}).(fnh{k});
        else
            vh='';
        end
        
        if isempty(vh)
            vh='';
        end
        if iscell(vh)
            vh=cell2mat(vh);
        end
        if isnumeric(vh)

            jtype.(fnu{f}).(fnh{k})=size(vh);
            if numel(vh)>1
                vh=reshape(vh,[1, numel(vh)]);
            end
            vh=num2str(vh);
        else
            jtype.(fnu{f}).(fnh{k})='s';
        end
        tab{ind,3}=vh;
        ind=ind+1;
    end
end

% what is a file?
jtype.SMAP.watchFolder='d';
jtype.InOut.calibration_file='3d';
jtype.InOut.experiment_out='d';
obj.jsontypes=jtype;

ta=obj.guihandles.parttable;
ta.Data=tab;
ta.ColumnName={'Class','Parameter','Value'};
ta.ColumnEditable=[false false true];
ta.CellSelectionCallback={@cellselect_callback,obj};
w=ta.Position(3);
ta.ColumnWidth={.1*w, w*0.25,w*0.55};
obj.yamlpar=js;
end

function table2json(obj)
dat=obj.guihandles.parttable.Data;
js=obj.yamlpar;
for k=1:size(dat,1)
    if ~isfield(js.(dat{k,1}),(dat{k,2}))
        js.(dat{k,1}).(dat{k,2})='';
    elseif isnumeric(obj.jsontypes.(dat{k,1}).(dat{k,2})) && ~isempty(dat{k,3})
       vh=str2num(dat{k,3});
       js.(dat{k,1}).(dat{k,2})=reshape(vh,obj.jsontypes.(dat{k,1}).(dat{k,2}));
    else
        js.(dat{k,1}).(dat{k,2})=(dat{k,3});
    end      
end
obj.yamlpar=js;
end

function last3d=get3dcalfile(obj)
 last3d=obj.getPar('cal_3Dfile');
     if isempty(last3d)
         if isfield(obj.locData.files.file.savefit.fitparameters,'MLE_GPU_Yiming')
            calftest=obj.locData.files.file.savefit.fitparameters.MLE_GPU_Yiming.cal_3Dfile;
         else
             calftest=obj.locData.files.file.savefit.fitparameters.MLE_global_spline.cal_3Dfile;
         end
        if exist(calftest,'file')
            last3d=calftest;
        end
     end
end

function cellselect_callback(tabobj,selection,obj)
if isempty(selection.Indices)
    return
end
cl=tabobj.Data{selection.Indices(1),1};
pn=tabobj.Data{selection.Indices(1),2};
type=obj.jsontypes.(cl).(pn);

if strcmp(pn,'set_emitters_per_um2')
    list={'LD','HD','UHD'};
    answ=listdlg('ListString',list,'SelectionMode','single','InitialValue',2);
    if answ>0
    tabobj.Data{selection.Indices(1),3}=list{answ};
    end
end

if ischar(type) && (strcmp(type,'3d')  || strcmp(type,'d')) && selection.Indices(2)<3 %allow for direct edit
   last3d='';
    if strcmp(type,'3d')
        last3d=get3dcalfile(obj);
        if isempty(last3d)
            lf=obj.getPar('lastSMLFile');
            if ~isempty(lf)
                [path, file]=fileparts(lf);
                last3d=[path filesep '*_3dcal.mat'];
            end
        end
    end
    if isempty(last3d)
        last3d=tabobj.Data{selection.Indices(1),3};
    end
    if strcmp(type,'d')
        path=uigetdir(fileparts(last3d));
        file='';
    else
        [file,path]=uigetfile(last3d);
    end
    if path
        tabobj.Data{selection.Indices(1),3}=[path file];
    
        if contains(file,'_3dcal.mat')
            obj.setPar('cal_3Dfile',[path file]);
        end   
    if strcmp(type,'3d')
        obj.yamlpar.InOut.calibration_file=[path file];
        setz(obj);
%         if isempty(obj.yamlpar.InOut.experiment_out)
            setnettrainingpath(0,0,obj);
%         end
        makejsontable(obj)
    end
    end
end
end

function load_callback(a,b,obj,ext)
switch ext
    case '.json'
        [file, path]=uigetfile(obj.jsonfile);
        if file
            obj.jsonfile=[path file];
            jsontxt=fileread(obj.jsondefault);
            obj.jsonstruct=jsondecode(jsontxt);
            makejsontable(obj);
        end          
    case '_3dcal.mat'
        table2json(obj);
        file3dcal=obj.jsonstruct.InOut.calibration_file;
        if ~exist(fileparts(file3dcal),'dir')
            file3dcal=obj.getPar('cal_3Dfile');
        end
        [file, path]=uigetfile(file3dcal);
        if file 
            if   contains(file,'_3dcal.mat')
                obj.jsonstruct.InOut.calibration_file=[path file];
                makejsontable(obj)
                obj.setPar('cal_3Dfile',[path file]);
            end
        end
    case '.yaml'
        [file, path]=uigetfile(obj.yamlfile);
        if file
            obj.yamlpar=ReadYaml([path file]);
            obj.yamlfile=[path file];
            makejsontable(obj);
        end   
end
end
function usecurrent_callback(a,b,obj)
    p=obj.getGuiParameters;
    table2json(obj)
    js=obj.yamlpar;
    
    fields={'filenumber','frame','phot','locprecnm','znm','PSFxnm','locprecznm','numberInGroup','bg','bg2'};
    [locs,~,hroi]=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','grouped');
    if isempty(locs.frame)
        warndlg('please load an SMLM data set and make sure that in the rendered image a sufficient number of localizations are displayed')
        return
    end
    locsu=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','ungrouped');
    f=figure('Visible','off');
    stat=make_statistics2({locs});
    close(f)

    js.Simulation.lifetime_avg=max(1,stat.lifetime.mu-1);      
    js.Simulation.intensity_mu_sig= [1,0.2]*stat.photons.meanphot/js.Simulation.lifetime_avg; %30% variation
    if stat.lifetime.mu<0.8 %if most fluorophores in1 frame only, but we simulate with longer on-time, increase brightness
        js.Simulation.intensity_mu_sig=js.Simulation.intensity_mu_sig*1.3;
    end
    bgminmax=quantile(locsu.bg,[0.01, 0.95]);
    dbg=bgminmax(2)-bgminmax(1);
    bgrange=bgminmax+ [-1, 1]*dbg*0.3;
    bgrange(1)=max(bgrange(1), quantile(locsu.bg,0.0005)*0.9);
    js.Simulation.bg_uniform=bgrange; %set a bit lower to allow for varying background
    
%     PSF for 2D learning
    if ~isempty(locsu.PSFxnm) && isempty(locsu.znm)
        PSFminmax=quantile(locsu.PSFxnm,[0.02, 0.98]);
        PSFmodal=stat.PSFsize.max;
        js.Simulation.PSF_min_mod_max=double([PSFminmax(1) PSFmodal PSFminmax(2)]);
    else
        js.Simulation.PSF_min_mod_max=[];
    end



    fi=obj.locData.files.file(1).info;
    if fi.EMon
        js.Camera.em_gain=fi.emgain;
    else
        js.Camera.em_gain=[];
    end
    js.Camera.e_per_adu=fi.conversion;
    js.Camera.px_size=fi.cam_pixelsize_um([2 1])*1000;
    js.Camera.baseline=fi.offset;
    if js.Camera.em_gain>0 
        js.Camera.read_sigma=74.4;
        disp('read noise set to 74.4 e- for EM gain')
    else
        js.Camera.read_sigma=1.5;
        disp('read noise set to 1.5 e- for sCMOS or non-EM camera')   
    end
    
    if isa(hroi,'imroi') && isvalid(hroi)
        m=hroi.createMask;
        area=sum(m(:))*obj.getPar('sr_pixrec')^2; % in nanometers
        
    else %FoV
        sizeim=obj.getPar('sr_size')*2;
        area=sizeim(1)*sizeim(end);
    end
    cnm=obj.getPar('cam_pixelsize_nm');
    areapix=area/cnm(1)/cnm(end);
    density=length(locs.frame)/areapix/(max(locs.frame)-min(locs.frame)); %in locs per pix^2 per frame
    emitters=density*40*40;
    disp(['emitters in 40x40 per frame: ' num2str(emitters)]);
    js.SMAP.density=density;
    [expdir,js.SMAP.name]=fileparts(obj.locData.files.file(1).name);
    if isempty(js.InOut.calibration_file)
        js.InOut.calibration_file=get3dcalfile(obj);
    end
    obj.yamlpar=js;  
    setz(obj);
    setnettrainingpath(0,0,obj)
end

function setnettrainingpath(a,b,obj)
    p=obj.getGuiParameters;
    js=obj.yamlpar;
    js.Hardware.device='cuda:0';
    js.Hardware.device_simulation='cuda:0';
    if p.trainlocal
        if ismac 
            js.Hardware.device='cpu';
            js.Hardware.device_simulation='cpu';
        end
%         if ~isempty(js.InOut.calibration_file)
            [calfile,calname]=fileparts(js.InOut.calibration_file);
            dates=datestr(now,'yymmdd');
            js.InOut.experiment_out=[calfile '/' dates '_local_' calname '/'];
%         end        
    else
        if ~isempty(js.InOut.calibration_file)
            [~,calname]=fileparts(js.InOut.calibration_file);
            dates=datestr(now,'yymmdd');
            js.InOut.experiment_out=[dates '_' calname '/'];
        end
    end
    obj.yamlpar=js;  
    makejsontable(obj)
end

function savejson_callback(a,b,obj,isdefault)
finalizejson(obj);
if isdefault
    fout=obj.yamldefault;
else
    [~,fname]=fileparts(obj.locData.files.file(1).name);
    fnameo=[obj.yamlpar.InOut.experiment_out filesep 'DECODE_train_' fname '.yaml'];
    
    [file, path]=uiputfile(fnameo);
    if ~file
        disp('no file saved')
        return
    end
    fout=[path file];
end

saveyaml(obj.yamlpar,fout)
end

function saveyaml(yamlpar,fout)
yamlpar=make2Dastig(yamlpar);
yout=rmfield(yamlpar,'SMAP');
WriteYamlSimple(fout, yout);
end

function yout=make2Dastig(yamlpar)
yout=yamlpar;
if ~isempty(yamlpar.Simulation.PSF_min_mod_max)
    % calculate PSF parameters
    PSFmodal=yamlpar.Simulation.PSF_min_mod_max(2);
    PSFmin=yamlpar.Simulation.PSF_min_mod_max(1);
    pixelsize=mean(yamlpar.Camera.px_size);
    w0=PSFmodal(1)*2;
    dz=50; 
    z=-dz+yamlpar.SMAP.zrange_nm(1):dz:yamlpar.SMAP.zrange_nm(2)+dz;
%     zR=pi*4*PSFmodal^2*1.4/600; % nm pi*w0^2+n/lambda, w0=2sigma,
    zR=pi*4*(PSFmin+PSFmodal)^2/4*1.4/600;% nm pi*w0^2+n/lambda, w0=2sigma,

    wz=w0*sqrt(1+(z./zR).^2);
    sigma_pix(1,1,:)=wz/2/pixelsize;

    % calculate volume PSF
    roisize=21;
    n=-(roisize-1)/2:(roisize-1)/2;
    [X,Y]=meshgrid(n);
    PSFvol=exp(-(X.^2+Y.^2)/2./sigma_pix.^2)/pi/2./sigma_pix.^2;

    %calculate cspline coefficients
     coeffr = single(Spline3D_interp(PSFvol));
    % Save PSF
    SXY.cspline.coeff{1}=coeffr;
    SXY.cspline.dz=dz;
    SXY.cspline.z0=ceil(size(PSFvol,3)/2);
    SXY.cspline.x0=round((roisize+1)/2);
    SXY.cspline.mirror=0;

%     parameters.midpoint=ceil(size(PSFvol,3)/2);
%     parameters.mirror=false;
%     parameters.dz=dz;
    parameters.zR=zR;
    parameters.w0=w0;
    fout=[yamlpar.InOut.experiment_out '/Gauss_3dcal.mat'];
    save(fout,'SXY','parameters')
    yout=yamlpar;
    yout.InOut.calibration_file=fout;
    yout.Simulation.emitter_extent{3}{1}=0;
end
end

function finalizejson(obj)
% write some default parameters to json that are not part of the table.
table2json(obj);
js=obj.yamlpar; 

list={'LD','HD','UHD'};
emit=[15 25 50];
density=find(strcmp(js.SMAP.set_emitters_per_um2,list));
js.Simulation.emitter_av=emit(density);
js.Simulation.emitter_extent{3}{1}=round(js.SMAP.zrange_nm(1));
js.Simulation.emitter_extent{3}{2}=round(js.SMAP.zrange_nm(2));
js.InOut.calibration_file=strrep(js.InOut.calibration_file,'\','/');
js.InOut.experiment_out=strrep(js.InOut.experiment_out,'\','/');

obj.yamlpar=js;         
end

function setz(obj)
 [locs,~,hroi]=obj.locData.getloc('znm','layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','grouped');
 if isempty(locs.znm)
     zminmax=[-750, 750];
 else
     z=quantile(locs.znm,[0.01,0.99]); 
     dz=z(2)-z(1);
     zminmax=z+[-1 1]*dz*0.2;
 end
 calf=obj.yamlpar.InOut.calibration_file;
 if ~isempty(calf)
     if ~exist(calf,'file')
         disp('please load PSF calibration file')
         [fn,pn]=uigetfile(calf,'load PSF calibration file');
         calf=[pn fn];
     end
    l=load(calf);
    obj.yamlpar.InOut.calibration_file=calf;
    if isempty(obj.yamlpar.InOut.experiment_out)
        obj.yamlpar.InOut.experiment_out=fileparts(calf);
    end
    if isfield(l,'parameters')
    zr=(l.parameters.fminmax(2)-l.parameters.fminmax(1))*l.parameters.dz/2;
    zminmax(1)=max(zminmax(1),-zr);
    zminmax(2)=min(zminmax(2),zr);
    end
    
 end
 obj.yamlpar.SMAP.zrange_nm=zminmax;
 
end

function stoplearning_callback(a,b,obj)
if obj.fitlocal
    pid = obj.decodeprocess;
    if ispc
        cmd = sprintf('taskkill /F /PID %d', pid);
    else
        cmd = sprintf('kill %d', pid);
    end
    system(cmd)
else
if ~isempty(obj.decodepid) %remote training
    if length(obj.decodepid)>1
        for k=1:length(obj.decodepid)
            ls{k}=[num2str(obj.decodepid(k).pid) ':' datestr(obj.decodepid(k).date)];
        end
        ind=listdlg('PromptString',  "Select process to end",'ListString',ls);
    else
        ind=1;
    end
    url = 'http://pc-ries25:8000/kill';
    options = weboptions('RequestMethod', 'post', 'ArrayFormat','json');
    for k=1:length(ind)    
        try
            webread(url, 'pid', obj.decodepid(ind(k)).pid, options);
             obj.decodepid(ind)=[];
        catch err
            err
        end
    end
   
    if isempty(obj.decodepid)
        obj.guihandles.stoplearning.Visible='off';
    end
else
    obj.tensorboardprocess.stop;
%     [err,gpustat]=system(['ssh ' obj.yamlpar.Connect.remote_workstation ' ''nvidia-smi '' ']);
%     [gpus,gpurec]=parsegpustat(gpustat);
%     for k=1:length(gpus.pid)
%         command=['ssh ' obj.yamlpar.Connect.remote_workstation ' ''kill ' num2str(gpus.pid(k).PID) ''''];
%         system(command);
%     end  
    obj.guihandles.stoplearning.Visible='off';
    obj.guihandles.tensorboard.Visible='off';
    try
        obj.decodeprocess.stop;
    catch err
        err
    end
end
end
end

function tensorboard_callback(a,b,obj)
if ~obj.fitlocal %remote training
   page='http://pc-ries25:8001/#scalars';
else
    page=['http://localhost:' obj.tensorboardprocess.port];
end
web(page,'-browser')
end

function pard=guidef(obj)
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Saves a training file for DECODE';

pard.usecurrent.object=struct('Style','pushbutton','String','Use parameters from current SML file','Callback',{{@usecurrent_callback,obj}});
pard.usecurrent.position=[1,1];
pard.usecurrent.Width=2;

pard.tensorboard.object=struct('Style','pushbutton','String','Tensorboard','Callback',{{@tensorboard_callback,obj}},'Visible','on');
pard.tensorboard.position=[1,3];
pard.tensorboard.Width=1;

pard.stoplearning.object=struct('Style','pushbutton','String','Stop training','Callback',{{@stoplearning_callback,obj}},'Visible','off');
pard.stoplearning.position=[1,4];
pard.stoplearning.Width=1;

pard.loadjson.object=struct('Style','pushbutton','String','Load yaml','Callback',{{@load_callback,obj,'.yaml'}});
pard.loadjson.position=[9,1];
pard.loadjson.Width=1;

pard.savejson_default.object=struct('Style','pushbutton','String','Save as default','Callback',{{@savejson_callback,obj,1}});
pard.savejson_default.position=[9,4];
pard.savejson_default.Width=1;

pard.savejson.object=struct('Style','pushbutton','String','Save yaml','Callback',{{@savejson_callback,obj,0}});
pard.savejson.position=[9,2];
pard.savejson.Width=1;

% pard.load3D_cal.object=struct('Style','pushbutton','String','Load 3d_cal','Callback',{{@load_callback,obj,'_3dcal.mat'}});
% pard.load3D_cal.position=[1,3];
% pard.load3D_cal.Width=1;
% 
% pard.setoutput.object=struct('Style','pushbutton','String','Output directory','Callback',{{@load_callback,obj,'*'}});
% pard.setoutput.position=[1,4];
% pard.setoutput.Width=1;

pard.partablepos.object=struct('Style','text','String',' ');
pard.partablepos.position=[8,1];
pard.partablepos.Width=4;

pard.trainlocal.object=struct('Style','checkbox','String','train locally');
pard.trainlocal.position=[10,1];
pard.trainlocal.Width=1;
end