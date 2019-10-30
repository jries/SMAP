classdef TrainDeepSMLM<interfaces.DialogProcessor
%     Saves a training file for deepSMLM
    properties
        jsonstruct
        jsondefault='settings/temp/model_deepSMLMdefault.json';
        jsonfile
        jsontypes
    end
    methods
        function obj=TrainDeepSMLM(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
           out=[];
           %set defaults
            outdir=finalizejson(obj);
            js=obj.jsonstruct;
           %make directory
           status=mkdir(outdir);
           %copy 3dcal and set 3dcal path relative
           status2=copyfile(js.InOut.calibration_file,outdir);
           %save jsonfile
           
           fileout=[outdir filesep 'model_' js.SMAP.name '.json'];
           if status && status2
                savejsonfile(obj,fileout)
           end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.DialogProcessor(obj);
            obj.jsonfile=obj.jsondefault;
            tt=uitable(obj.handle);
            tt.Position=obj.guihandles.partablepos.Position;
            tt.Position(4)=tt.Position(4)*7;
            obj.guihandles.parttable=tt;
            jsontxt=fileread(obj.jsondefault);
            obj.jsonstruct=jsondecode(jsontxt);
            makejsontable(obj);
        end

    end
end

function makejsontable(obj)
use.SMAP={'watchFolder','name','set_emitters_per_um2'};
use.InOut={'calibration_file'};
use.Simulation={'density','intensity_mu_sig','lifetime_avg','bg_uniform','bg_perlin_scale','bg_perlin_amplitude'};
use.Camera={'em_gain','read_sigma','spur_noise','e_per_adu','px_size'};
use.Logging={'log_comment'};
use.Scaling={'z_max','phot_max','bg_max'};
fnu=fieldnames(use);

% 
% for f=1:length(fnu)
%     isf.(fnu{f})=false(length(use.(fnu{f})));
% end
% %set if file
% isf.SMAP=[1 1]; %this is a file
% isf.InOut(1)=1; %this is a file

js=obj.jsonstruct;

if ~isfield(js, 'SMAP')
    js.SMAP.watchFolder='';
    js.SMAP.name='';
    js.SMAP.set_emitters_per_um2=0;
end

ind=1;
for f=1:length(fnu)
    fnh=use.(fnu{f});
    for k=1:length(fnh)
        tab{ind,1}=fnu{f};
        tab{ind,2}=fnh{k};
        vh=js.(fnu{f}).(fnh{k});
        if isnumeric(vh)
            jtype.(fnu{f}).(fnh{k})=size(vh);
            if numel(vh)>1
                vh=reshape(vh,[1, numel(vh)]);
            end
        else
            jtype.(fnu{f}).(fnh{k})='s';
        end
        tab{ind,3}=num2str(vh);
        ind=ind+1;
    end
end

% what is a file?
jtype.SMAP.watchFolder='d';
jtype.InOut.calibration_file='rf';
obj.jsontypes=jtype;

ta=obj.guihandles.parttable;
ta.Data=tab;
ta.ColumnName={'Class','Parameter','Value'};
ta.ColumnEditable=[false false true];
ta.CellSelectionCallback={@cellselect_callback,obj};
w=ta.Position(3);
ta.ColumnWidth={.1*w, w*0.25,w*0.55};

obj.jsonstruct=js;

end

function table2json(obj)
dat=obj.guihandles.parttable.Data;
js=obj.jsonstruct;
for k=1:size(dat,1)
    vjs=js.(dat{k,1}).(dat{k,2});
    if isnumeric(vjs)
       vh=str2num(dat{k,3});
       js.(dat{k,1}).(dat{k,2})=reshape(vh,obj.jsontypes.(dat{k,1}).(dat{k,2}));
    else
        js.(dat{k,1}).(dat{k,2})=(dat{k,3});
    end      
end
obj.jsonstruct=js;
end

function cellselect_callback(tabobj,selection,obj)
if isempty(selection.Indices)
    return
end
cl=tabobj.Data{selection.Indices(1),1};
pn=tabobj.Data{selection.Indices(1),2};
type=obj.jsontypes.(cl).(pn);

if ischar(type) && (strcmp(type,'af') || strcmp(type,'rf') || strcmp(type,'d')) && selection.Indices(2)<3 %allow for direct edit
    last3d=obj.getPar('cal_3Dfile');
    if isempty(last3d)
        lf=obj.getPar('lastSMLFile');
        if ~isempty(lf)
            [path, file]=fileparts(lf);
            last3d=[path filesep '*_3dcal.mat'];
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
    end
end
end
function progress_callback(a,b,obj)
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
end

end
function usecurrent_callback(a,b,obj)
    table2json(obj)
    js=obj.jsonstruct;
    
    fields={'filenumber','frame','phot','locprecnm','znm','PSFxnm','locprecznm','numberInGroup','bg','bg2'};
    [locs,~,hroi]=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','grouped');
    if isempty(locs.frame)
        warndlg('please load an SMLM data set and make sure that in the rendered image a sufficient number of localizations are displayed')
        return
    end
    locsu=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','ungrouped');
    stat=make_statistics2({locs});
    
    
    
    js.Simulation.lifetime_avg=stat.lifetime.mu-1;
    js.Simulation.intensity_mu_sig= [1;0.2]*stat.photons.meanphot/js.Simulation.lifetime_avg; %30% variation
    bgminmax=quantile(locsu.bg,[0.05; 0.95]);
    dbg=bgminmax(2)-bgminmax(1);
    js.Simulation.bg_uniform=bgminmax-dbg*0.2; %set a bit lower to allow for varying background
    js.Simulation.bg_perlin_amplitude=dbg*[1; 0.5];
    
    fi=obj.locData.files.file(1).info;
    js.Camera.em_gain=fi.emgain;
    js.Camera.e_per_adu=fi.conversion;
    js.Camera.px_size=fi.cam_pixelsize_um*1000;
    
    z=quantile(locs.znm,[0.1 0.9]);
    js.Scaling.z_max=max(abs(z))*1.3;
    js.Scaling.phot_max=quantile(locsu.phot,.95)*2;
    js.Scaling.bg_max=bgminmax(2)+dbg*2+max(js.Simulation.bg_perlin_amplitude);
    %true simulation lifetime = mean(lifetime_grouped)-1
    %density, based on all
%     [locsall,~,hroi]=obj.locData.getloc({'frame'},'position','roi','grouping','grouped');
%     hroi=obj.getPar('sr_roihandle');
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
    emitters=density*64*64;
    disp(['emitters in 64x64 per frame: ' num2str(emitters)]);
    js.Simulation.density=density;
    [~,js.SMAP.name]=fileparts(obj.locData.files.file(1).name);
    
    obj.jsonstruct=js;  
    makejsontable(obj)
end
function savejson_callback(a,b,obj,isdefault)
finalizejson(obj);
if isdefault
    fout=obj.jsondefault;
else
    [file, path]=uiputfile([obj.jsonfile]);
    if ~file
        disp('no file saved')
        return
    end
    fout=[path file];
end

savejsonfile(obj,fout)
end

function savejsonfile(obj,fout)
jsontxt=jsonencode(obj.jsonstruct);
jsontxt = strrep(jsontxt, ',', sprintf(',\r'));
jsontxt = strrep(jsontxt, '[{', sprintf('[\r{\r'));
jsontxt = strrep(jsontxt, '}]', sprintf('\r}\r]'));
fid=fopen(fout,'w');
fprintf(fid,jsontxt);
fclose(fid);
end

function [outdir]=finalizejson(obj)
% write some default parameters to json that are not part of the table.
table2json(obj);
js=obj.jsonstruct; 
       %overwrite density by set_emitters_per_um2
outdir=[js.SMAP.watchFolder filesep 'train_' datestr(now,'yymmdd') '_' js.SMAP.name];
js.InOut.model_out=['model_' js.SMAP.name '.pt'];   %set obj.jsonstruct.InOut.model_out

if js.SMAP.set_emitters_per_um2>0
    px=js.Camera.px_size;
    density=js.SMAP.set_emitters_per_um2/1000/1000*px(1)*px(2);
    disp(['new density per pixel^2: ' num2str(density)])
    js.Simulation.density=density;
end
obj.jsonstruct=js;
         
end

function pard=guidef(obj)
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Saves a training file for deepSMLM';

pard.progress.object=struct('Style','pushbutton','String','Progress','Callback',{{@progress_callback,obj}});
pard.progress.position=[1,4];
pard.progress.Width=1;

pard.usecurrent.object=struct('Style','pushbutton','String','Use parameters from current SML file','Callback',{{@usecurrent_callback,obj}});
pard.usecurrent.position=[1,1];
pard.usecurrent.Width=2;


pard.loadjson.object=struct('Style','pushbutton','String','Load Json','Callback',{{@load_callback,obj,'.json'}});
pard.loadjson.position=[9,1];
pard.loadjson.Width=1;


pard.savejson_default.object=struct('Style','pushbutton','String','Save as default','Callback',{{@savejson_callback,obj,1}});
pard.savejson_default.position=[9,4];
pard.savejson_default.Width=1;

pard.savejson.object=struct('Style','pushbutton','String','Save Json','Callback',{{@savejson_callback,obj,0}});
pard.savejson.position=[9,2];
pard.savejson.Width=1;

pard.load3D_cal.object=struct('Style','pushbutton','String','Load 3d_cal','Callback',{{@load_callback,obj,'_3dcal.mat'}});
pard.load3D_cal.position=[1,3];
pard.load3D_cal.Width=1;

pard.partablepos.object=struct('Style','text','String',' ');
pard.partablepos.position=[8,1];
pard.partablepos.Width=4;

end