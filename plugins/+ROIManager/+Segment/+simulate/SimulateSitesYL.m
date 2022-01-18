classdef SimulateSitesYL<interfaces.DialogProcessor&interfaces.SEProcessor
    %SimulateSites is a localization based simulation engine for SMAP. It
    %uses as an input a list of localizations, a matlab function that
    %returns coordiantes or an image which defines a 2D structure. It
    %returns simulated localizations to SMAP using a realistic model for the
    %photophysics of the dye. Simulated structures are added to the
    %RoiManager
    properties
    end
    methods
        function obj=SimulateSitesYL(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_sitefov','se_cellpixelsize','se_siteroi'};
            obj.history=true;
        end
        function initGui(obj)
            setvisibility(obj);
        end
        function out=run(obj,p)
            ptime = p.time;
            if length(ptime)==1
                runSingleSim(obj,p);
            else
                allTime = ptime(1):ptime(2):ptime(3);
                for i=1:length(allTime)
                    p.time=allTime(i);
                    runSingleSim(obj,p);
                end
            end
            out = [];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function out=runSingleSim(obj,p)  
[locst,possites]=simulatelocs2(p, p.channel(1));
if (p.pro2nd.selection)~= "No pes"
    [locst2,~]=simulatelocs2(p, p.channel(2));
    locst = mergeStruct(locst, locst2);
end

if ~p.savez
   locst=rmfield(locst,{'znm','znm_gt'});
end

    info = [num2str(p.channel(1)) p.selectPro.selection num2str(p.channel(2)) p.pro2nd.selection 'T' num2str(p.time)];
    labels=num2str(p.labeling_efficiency*100,'%2.0f');
    phots=num2str(p.photons,'%3.0f');
    blinks=num2str(p.blinks,'%3.0f');
    filename=['L' labels 'P' phots 'B' blinks 'Info' info];


obj.locData.addfile(['simulated_' num2str(obj.locData.files.filenumberEnd) '_' filename]);
obj.locData.files.file(end).info.simulationParameters=obj.getGuiParameters;
obj.locData.addLocData(locst);
obj.locData.sort('filenumber','frame');
try
initGuiAfterLoad(obj);
catch err
   err
end
se=obj.locData.SE;
cell=interfaces.SEsites;
cell.pos=[mean(locst.xnm) mean(locst.ynm)];
cell.ID=0;
cell.info.filenumber=obj.locData.files.filenumberEnd;
se.addCell(cell);
for k=1:length(possites)
   thissite=interfaces.SEsites;
   thissite.pos=[possites(k).x possites(k).y];
   thissite.info.cell=cell.ID;
   thissite.info.filenumber=cell.info.filenumber;
    % thissite.cellnumber=sitepar.currentcell.number;
%         thissite.number=sitepar.sitelist.cellnumber+1;
    se.addSite(thissite);
end 

%            obj.setPar('SimulateSitesParameters',p);

try
se.currentsite=se.sites(1);
se.currentcell=se.cells(1);
se.currentfile=se.files(1);
se.processors.preview.updateFilelist;
se.processors.preview.updateCelllist;
se.processors.preview.updateSitelist; 
se.processors.preview.nextsite(1)
%            se.processors.plotsite(se.sites(1));
catch err
    err
end

if p.savenow.Value==2
    obj.locData.loc.xnm=obj.locData.loc.xnm_gt;
    obj.locData.loc.ynm=obj.locData.loc.ynm_gt;
    obj.locData.loc.znm=obj.locData.loc.znm_gt;
    obj.locData.loc=rmfield(obj.locData.loc,{'xnm_gt','ynm_gt','znm_gt'});
end
if p.savenow.Value>1
    lastsml=obj.getPar('lastSMLFile');
    if ~isempty(lastsml)
        path=fileparts(lastsml);
    else
        path='';
    end
   [file,path]= uiputfile([path filesep obj.locData.files.file(1).name]);
   if file
        obj.locData.savelocs([path file])
        obj.setPar('lastSMLFile',[path file]);
   end
end
out=[];
end

function load_callback(a,b,obj)
f=obj.getSingleGuiParameter('coordinatefile'); 
[f,p]=uigetfile({'*.*';'*.tif';'*.png';'*.csv';'*.txt';'*.mat'},'Choose coordinate file',f);
if ~f
    return
end
obj.setGuiParameters(struct('coordinatefile',[p f]));
setvisibility(obj)
end
function setvisibility(obj)
f=obj.getSingleGuiParameter('coordinatefile');
[p,fh,ext]=fileparts(f);
switch ext
    case {'.txt','.csv'}
        txt='on';
        tif='off';
    case {'.tif','.png'}
        txt='off';
        tif='on';
    case '.mat'
        l=load([p f]);
        if isfield(l,'image')
            txt='off';
            tif='on';
        else
            txt='on';
            tif='off';            
        end
    case '.m'
        cf=pwd;
        cd(p)
        [~,fh]=fileparts(f);
        l=eval(fh);
        cd(cf);
        if isfield(l,'image')
            txt='off';
            tif='on';
        else
            txt='on';
            tif='off';            
        end
end
obj.guihandles.labeling_efficiency.Visible=txt;
obj.guihandles.t_labelingefficiency.Visible=txt;
obj.guihandles.tif_density.Visible=tif;
obj.guihandles.tif_numbermode.Visible=tif;
obj.guihandles.tif_imagesizet.Visible=tif;
obj.guihandles.tif_imagesize.Visible=tif;
obj.guihandles.load_button.Visible='on';
obj.guihandles.coordinatefile.Visible='on';
obj.guihandles.pro2nd.Visible = 'off';
end

function pesInput_callback(a,b,obj)
if obj.getSingleGuiParameter('usePes')
    txt = 'off';
    tif = 'off';
    pes = 'on';
    obj.guihandles.labeling_efficiency.Visible=txt;
    obj.guihandles.t_labelingefficiency.Visible=txt;
    obj.guihandles.tif_density.Visible=tif;
    obj.guihandles.tif_numbermode.Visible=tif;
    obj.guihandles.tif_imagesizet.Visible=pes;
    obj.guihandles.tif_imagesize.Visible=pes;
    obj.guihandles.load_button.Visible='off';
    obj.guihandles.coordinatefile.Visible='off';
    obj.guihandles.time.Visible=pes;
    obj.guihandles.t_time.Visible=pes;
    obj.guihandles.selectPro.Visible=pes;
    obj.guihandles.pro2nd.Visible='off';
    if isfield(obj.P.par, 'proteinES')
        obj.guihandles.selectPro.String = fieldnames(obj.getPar('proteinES'));
    end
else
    pes = 'off';
    setvisibility(obj)
    obj.guihandles.time.Visible=pes;
    obj.guihandles.t_time.Visible=pes;
    obj.guihandles.selectPro.Visible=pes;
    obj.guihandles.pro2nd.Visible=pes;
    obj.guihandles.setPro2nd.Value = 0;
end
end

function pro2_callback(a,b,obj)
if obj.getSingleGuiParameter('setPro2nd')
    lenOptions=numel(obj.guihandles.selectPro.String(1:end ~= obj.guihandles.selectPro.Value));
    if lenOptions < 1
        obj.guihandles.pro2nd.String = {'No pes'};
    else
        obj.guihandles.pro2nd.String = obj.guihandles.selectPro.String(1:end ~= obj.guihandles.selectPro.Value);
    end
    obj.guihandles.pro2nd.Visible = 'on';
else
    obj.guihandles.pro2nd.Visible = 'off';
end
end

function pard=guidef(obj)

pard.t_selectPro.object=struct('String', 'Protein','Style','text');
pard.t_selectPro.position=[1,1];
pard.t_selectPro.Width=0.5;

pard.setPro2nd.object=struct('String', '2nd', 'Style', 'checkbox', 'Value', 0, 'Callback', {{@pro2_callback,obj}});
pard.setPro2nd.position=[1,2];
pard.setPro2nd.Width=0.5;

pard.selectPro.object=struct('String',{{'No pes'}},'Style','popup');
pard.selectPro.position=[1,1.5];
pard.selectPro.Width=0.5;

pard.pro2nd.object=struct('String',{{'No pes'}},'Style','popup');
pard.pro2nd.position=[1,2.5];
pard.pro2nd.Width=0.5;

pard.coordinatefile.object=struct('String','plugins/+ROIManager/+Segment/hidden/MakeNPCCoordinates.m','Style','edit');
pard.coordinatefile.position=[1,1];
pard.coordinatefile.Width=3;
pard.coordinatefile.TooltipString=sprintf('.txt or .csv file with coordinates, .tif or .png file in which the pixel values are a \n measure for the concentration of the labels or a matlab function that \n returns the position of the lables in form of a structure with the fields .x .y .z');

pard.load_button.object=struct('String','Load','Style','pushbutton','Callback',{{@load_callback,obj}});
pard.load_button.position=[1,4];
pard.load_button.TooltipString=pard.coordinatefile.TooltipString;

pard.usePes.object=struct('String','Use pes','Style','checkbox','Value',0,'Callback',{{@pesInput_callback,obj}});
pard.usePes.position=[2,4];
pard.usePes.TooltipString=pard.coordinatefile.TooltipString;

pard.tif_numbermode.object=struct('String',{{'Density (labels/um^2)','Number of labels'}},'Style','popupmenu');
pard.tif_numbermode.Width=1.5;
pard.tif_numbermode.position=[3,1];

pard.tif_density.object=struct('String','100','Style','edit');
pard.tif_density.position=[3,2.5];
pard.tif_density.Width=0.5;

pard.tif_imagesizet.object=struct('String','Image width (nm)','Style','text');
pard.tif_imagesizet.Width=1.5;
pard.tif_imagesizet.position=[3,3];

pard.tif_imagesize.object=struct('String','300','Style','edit');
pard.tif_imagesize.position=[3,4.5];
pard.tif_imagesize.Width=0.5;

pard.t_labelingefficiency.object=struct('String','Labeling efficiency','Style','text');
pard.t_labelingefficiency.position=[3,1];
pard.t_labelingefficiency.Width=1.5;

pard.labeling_efficiency.object=struct('String','.5','Style','edit');
pard.labeling_efficiency.Width=.5;
pard.labeling_efficiency.position=[3,2.5];

pard.t_time.object=struct('String','Time point','Style','text');
pard.t_time.position=[3,1];
pard.t_time.Width=1.5;

pard.time.object=struct('String','1','Style','edit');
pard.time.Width=.5;
pard.time.position=[3,2.5];

pard.t2.object=struct('String','mean re-activations','Style','text');
pard.t2.position=[4,1];
pard.t2.Width=1.5;

pard.blinks.object=struct('String','1','Style','edit');
pard.blinks.Width=.5;
pard.blinks.position=[4,2.5];


pard.t3.object=struct('String','lifetime (frames)','Style','text');
pard.t3.position=[4,3];
pard.t3.Width=1.5;

pard.lifetime.object=struct('String','1','Style','edit');
pard.lifetime.Width=.5;
pard.lifetime.position=[4,4.5];

pard.t4.object=struct('String','mean number photons','Style','text');
pard.t4.position=[5,1];
pard.t4.Width=1.5;

pard.photons.object=struct('String','2000','Style','edit');
pard.photons.Width=.5;
pard.photons.position=[5,2.5];

pard.t5.object=struct('String','BG per pixel (photons)','Style','text');
pard.t5.position=[5,3];
pard.t5.Width=1.5;

pard.background.object=struct('String','20','Style','edit');
pard.background.Width=.5;
pard.background.position=[5,4.5];

pard.t_channel.object=struct('String', 'Channel','Style','text');
pard.t_channel.position=[6,1];
pard.t_channel.Width=1.5;

pard.channel.object=struct('String', '1','Style','edit');
pard.channel.position=[6,2.5];
pard.channel.Width=0.5;

pard.t_viewType.object=struct('String', 'View','Style','text');
pard.t_viewType.position=[6,3];
pard.t_viewType.Width=1.5;

pard.viewType.object=struct('String', {{'side', 'top'}},'Style','popup','Value', 1);
pard.viewType.position=[6,4.5];
pard.viewType.Width=0.5;

pard.t6.object=struct('String','Number of sites','Style','text');
pard.t6.position=[8,1];
pard.t6.Width=1.5;

pard.randomrot.object=struct('String','Random rotation','Style','checkbox');
pard.randomrot.Width=1;
pard.randomrot.position=[7,1];

pard.randomxy.object=struct('String','Random position (nm):','Style','checkbox');
pard.randomxy.Width=1.5;
pard.randomxy.position=[7,2];

pard.randomxyd.object=struct('String','20','Style','edit');
pard.randomxyd.Width=.5;
pard.randomxyd.position=[7,3.5];

pard.savez.object=struct('String','save z','Style','checkbox','Value',1);
pard.savez.Width=1;
pard.savez.position=[2,3];

pard.numberofsites.object=struct('String','6 5','Style','edit');
pard.numberofsites.Width=.5;
pard.numberofsites.position=[8,2.5];


pard.t7.object=struct('String','number of frames','Style','text');
pard.t7.position=[8,3];
pard.t7.Width=1.5;

pard.maxframes.object=struct('String','3000','Style','edit');
pard.maxframes.Width=.5;
pard.maxframes.position=[8,4.5];

pard.savenow.object=struct('String',{{'No saving','Save ground truth','Save simulated locs'}},'Style','popupmenu');
pard.savenow.Width=2;
pard.savenow.position=[2,1];
pard.plugininfo.type='ROI_Analyze';

pard.plugininfo.description=sprintf(['SimulateSites is a localization based simulation engine for SMAP. It \n'...
    'uses as an input a list of localizations, a matlab function that\n'...
    'returns coordiantes or an image which defines a 2D structure. It \n'...
    'returns simulated localizations to SMAP using a realistic model for the \n'...
    'photophysics of the dye. Simulated structures are added to the RoiManager']);
end

function mergedStruct = mergeStruct(A, B)
mergedStruct=A;
f = fieldnames(A);
    for i = 1:length(f)
        mergedStruct.(f{i}) = [A.(f{i});B.(f{i})];
    end
end