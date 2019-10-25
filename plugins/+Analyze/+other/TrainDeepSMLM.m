classdef TrainDeepSMLM<interfaces.DialogProcessor
%     Saves a training file for deepSMLM
    properties
        jsonstruct
        jsondefault='settings/temp/model_deepSMLMdefault.json';
        isfile
    end
    methods
        function obj=TrainDeepSMLM(varargin)    
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        
        function out=run(obj,p)
           
           
%             makejsontable(obj)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.DialogProcessor(obj);
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
use.InOut={'calibration_file'};
use.Simulation={'density','photon_range','lifetime_avg'};
use.Camera={'em_gain','read_sigma','spur_noise'};

isf.InOut=1;
isf.Simulation=[0 0 0];
isf.Camera=[0 0 0];

fnu=fieldnames(use);
js=obj.jsonstruct;
ind=1;
for f=1:length(fnu)
    fnh=use.(fnu{f});
    for k=1:length(fnh)
        tab{ind,1}=fnu{f};
        tab{ind,2}=fnh{k};
        tab{ind,3}=num2str(js.(fnu{f}).(fnh{k}));
        obj.isfile(ind) = isf.(fnu{f})(k);
        ind=ind+1;
    end
end
ta=obj.guihandles.parttable;
ta.Data=tab;
ta.ColumnName={'Class','Parameter','Value'};
ta.ColumnEditable=[false false true];
ta.CellSelectionCallback={@cellselect_callback,obj};
w=ta.Position(3);
ta.ColumnWidth={.1*w, w*0.15,w*0.65};
end

function table2json(obj)
dat=obj.guihandles.parttable.Data;
js=obj.jsonstruct;
for k=1:size(dat,1)
    vjs=js.(dat{k,1}).(dat{k,2});
    if isnumeric(vjs)
       js.(dat{k,1}).(dat{k,2})=str2num(dat{k,3});
    else
        js.(dat{k,1}).(dat{k,2})=(dat{k,3});
    end      
end
end

function cellselect_callback(tabobj,selection,obj)
if obj.isfile(selection.Indices(1))&&selection.Indices(2)<3 %allow for direct edit
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
    [file,path]=uigetfile(last3d);
    if file
        tabobj.Data{selection.Indices(1),3}=[path file];
        obj.setPar('cal_3Dfile',[path file]);
    end   
end
end
function progress_callback(a,b,obj)
end
function load_callback(a,b,obj,ext)
end
function usecurrent_callback(a,b,obj)
    table2json(obj)
    js=obj.jsonstruct;
    
    fields={'filenumber','frame','phot','locprecnm','znm','PSFxnm','locprecznm','numberInGroup','bg','bg2'};
    locs=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','grouped');
    locsu=obj.locData.getloc(fields,'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','ungrouped');
    stat=make_statistics2({locs});
    
    stat.photons.meanphot
    js.Simulation.lifetime_avg=stat.lifetime.mu-1;
    %true simulation lifetime = mean(lifetime_grouped)-1
    %density, based on all
    [locsall,~,hroi]=obj.locData.getloc({'frame'},'position','roi','grouping','grouped');
    if ~isempty(hroi) && isvalid(hroi)
        m=hroi.createMask;
        area=sum(m(:))*obj.getPar('sr_pixrec')^2; % in nanometers
        cnm=obj.getPar('cam_pixelsize_nm');
    else %FoV
        area=obj.getPar('sr_size')^2*4;
    end
    areapix=area/cnm(1)/cnm(end);
    js.Simulation.density=length(locsall.frame)/areapix/(max(locsall.frame)-min(locsall.frame)); %in locs per pix^2 per frame
    
    makejsontable(obj)
end
function savejson_callback(a,b,obj)
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


pard.savejson.object=struct('Style','pushbutton','String','Save as default','Callback',{{@savejson_callback,obj}});
pard.savejson.position=[9,2];
pard.savejson.Width=1;

pard.load3D_cal.object=struct('Style','pushbutton','String','Load 3d_cal','Callback',{{@load_callback,obj,'_3dcal.mat'}});
pard.load3D_cal.position=[1,3];
pard.load3D_cal.Width=1;

pard.partablepos.object=struct('Style','text','String',' ');
pard.partablepos.position=[8,1];
pard.partablepos.Width=4;

end