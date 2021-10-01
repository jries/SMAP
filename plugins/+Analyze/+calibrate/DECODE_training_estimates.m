classdef DECODE_training_estimates<interfaces.DialogProcessor
%     Saves a training file for DECODE
    properties
        yamldefault='DECODE_default.yaml';
        yamlfile
        yamlpar
        jsontypes
        decodeprocess
        tensorboardprocess
    end
    methods
        function obj=DECODE_training_estimates(varargin)  
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=false;
        end
        function out=run(obj,p)
           pdecode=obj.getGlobalSetting('DECODE_path');
           if ~exist(pdecode,'dir')
               warning('Decode not found, please specify in the SMAP/Preferences menu in the Plugin tab.');
               return
           end
           out=[];
           [~,fname]=fileparts(obj.locData.files.file(1).name);
           fnameo=[obj.yamlpar.InOut.experiment_out filesep 'DECODE_train_' fname '.yaml'];
           finalizejson(obj);
           saveyaml(obj,fnameo);
%            https://github.com/brian-lau/MatlabProcessManager
           pcall=[pdecode '/bin/python -m decode.neuralfitter.train.live_engine -p ' fnameo ];
           pm=processManager('command',pcall,'autoStart',false,'workingDir',obj.yamlpar.InOut.experiment_out);
           pm.printStdout=false ;
           pm.printStderr=true;
           pm.wrap=1000;
           pm.pollInterval=10;
           pm.start()
           obj.decodeprocess=pm;

           pcalltb=[pdecode '/bin/tensorboard --samples_per_plugin images=100 --port=6007 --logdir=runs' ];
           ptb=processManager('command',pcalltb,'workingDir',obj.yamlpar.InOut.experiment_out);
           obj.tensorboardprocess=ptb;
           obj.guihandles.stoplearning.Visible='on';
           obj.guihandles.tensorboard.Visible='on';
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)

            initGui@interfaces.DialogProcessor(obj);
                obj.createGlobalSetting('DECODE_path','Plugins','The anaconda environmet path of decode (eg. /decode_env/):',struct('Style','dir','String','decode')) 
        
            yamldefault=[obj.getPar('SettingsDirectory') filesep 'cameras' filesep obj.yamldefault];
            obj.yamlfile=yamldefault;
            tt=uitable(obj.handle);
            tt.Position=obj.guihandles.partablepos.Position;
            tt.Position(4)=tt.Position(4)*7;
            obj.guihandles.parttable=tt;
            obj.yamlpar=ReadYaml(yamldefault);
            obj.yamlfile=yamldefault;
            makejsontable(obj);
        end

    end
end

function makejsontable(obj)
use.SMAP={'set_emitters_per_um2','zrange_nm'};
use.InOut={'calibration_file','experiment_out'};
use.Simulation={'intensity_mu_sig','lifetime_avg','bg_uniform'};
use.Camera={'em_gain','e_per_adu','baseline','read_sigma','spur_noise','px_size'};
use.Hardware={'device','device_simulation'};
fnu=fieldnames(use);

js=obj.yamlpar;
if ~isfield(js, 'SMAP')
    js.SMAP.set_emitters_per_um2='HD';
    js.SMAP.zrange_nm=[-0 0];
end

ind=1;
for f=1:length(fnu)
    fnh=use.(fnu{f});
    for k=1:length(fnh)
        tab{ind,1}=fnu{f};
        tab{ind,2}=fnh{k};
        vh=js.(fnu{f}).(fnh{k});
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
    vjs=js.(dat{k,1}).(dat{k,2});
    if isnumeric(obj.jsontypes.(dat{k,1}).(dat{k,2}))
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
     calftest=obj.locData.files.file.savefit.fitparameters.MLE_GPU_Yiming.cal_3Dfile;
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
        if isempty(obj.yamlpar.InOut.experiment_out)
            obj.yamlpar.InOut.experiment_out=path;
        end
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
    table2json(obj)
    js=obj.yamlpar;
    
    fields={'filenumber','frame','phot','locprecnm','znm','PSFxnm','locprecznm','numberInGroup','bg','bg2'};
    [locs,~,hroi]=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','grouped');
    if isempty(locs.frame)
        warndlg('please load an SMLM data set and make sure that in the rendered image a sufficient number of localizations are displayed')
        return
    end
    locsu=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','ungrouped');
    f=figure;
    stat=make_statistics2({locs});
    close(f)
    
    
    
    js.Simulation.lifetime_avg=stat.lifetime.mu-1;
    js.Simulation.intensity_mu_sig= [1,0.2]*stat.photons.meanphot/js.Simulation.lifetime_avg; %30% variation
    bgminmax=quantile(locsu.bg,[0.01, 0.95]);
    dbg=bgminmax(2)-bgminmax(1);
    bgrange=bgminmax+ [-1, 1]*dbg*0.3;
    bgrange(1)=max(bgrange(1), quantile(locsu.bg,0.0005)*0.9);
    js.Simulation.bg_uniform=bgrange; %set a bit lower to allow for varying background

    
    fi=obj.locData.files.file(1).info;
    js.Camera.em_gain=fi.emgain*fi.EMon;
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
    if ismac
        js.Hardware.device='cpu';
        js.Hardware.device_simulation='cpu';
    else
        js.Hardware.device='cuda:0';
        js.Hardware.device_simulation='cuda:0';
    end
    if isempty(js.InOut.calibration_file)
        js.InOut.calibration_file=get3dcalfile(obj);
    end
%     if ~isempty(js.InOut.calibration_file)
%         js.InOut.experiment_out=fileparts(js.InOut.calibration_file);
%     end
    obj.yamlpar=js;  
    setz(obj);
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

saveyaml(obj,fout)
end

function saveyaml(obj,fout)
yout=rmfield(obj.yamlpar,'SMAP');
WriteYaml(fout, yout)
end

% function savejsonfile(obj,fout)
% jsontxt=jsonencode(obj.jsonstruct);
% jsontxt = strrep(jsontxt, ',', sprintf(',\r'));
% jsontxt = strrep(jsontxt, '[{', sprintf('[\r{\r'));
% jsontxt = strrep(jsontxt, '}]', sprintf('\r}\r]'));
% fid=fopen(fout,'w');
% fprintf(fid,jsontxt);
% fclose(fid);
% end

function finalizejson(obj)
% write some default parameters to json that are not part of the table.
table2json(obj);
js=obj.yamlpar; 

list={'LD','HD','UHD'};
emit=[15 25 50];
density=find(strcmp(js.SMAP.set_emitters_per_um2,list));
js.Simulation.emitter_av=emit(density);
js.Simulation.emitter_extent{3,1}=round(js.SMAP.zrange_nm(1));
js.Simulation.emitter_extent{3,2}=round(js.SMAP.zrange_nm(2));

js.InOut.calibration_file=strrep(js.InOut.calibration_file,'\','/');
js.InOut.experiment_out=strrep(js.InOut.experiment_out,'\','/');

obj.yamlpar=js;
         
end

function setz(obj)
 [locs,~,hroi]=obj.locData.getloc('znm','layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','grouped');
%  obj.yamlpar.SMAP.zrange_nm=[-750, 750];
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
    zr=(l.parameters.fminmax(2)-l.parameters.fminmax(1))*l.parameters.dz/2;
    zminmax(1)=max(zminmax(1),-zr);
    zminmax(2)=min(zminmax(2),zr);
    obj.yamlpar.SMAP.zrange_nm=zminmax;
 end
end

function stoplearning_callback(a,b,obj)
    obj.tensorboardprocess.stop;
    obj.decodeprocess.stop;
    obj.guihandles.stoplearning.Visible='off';
    obj.guihandles.tensorboard.Visible='off';
end

function tensorboard_callback(a,b,obj)
page='http://localhost:6007';
web(page,'-browser')
end

function pard=guidef(obj)
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Saves a training file for DECODE';

pard.usecurrent.object=struct('Style','pushbutton','String','Use parameters from current SML file','Callback',{{@usecurrent_callback,obj}});
pard.usecurrent.position=[1,1];
pard.usecurrent.Width=2;

pard.tensorboard.object=struct('Style','pushbutton','String','Tensorboard','Callback',{{@tensorboard_callback,obj}},'Visible','off');
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

end