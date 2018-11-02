classdef SEGUISettings< interfaces.SEProcessor
    properties
        SEpreview
    end
    methods
        function obj=SEGUISettings(varargin)
            obj@interfaces.SEProcessor(varargin{:})    
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.SEProcessor(obj);
            set(obj.guihandles.showSE,'Callback',{@make_siteexplorer,obj})
            set(obj.guihandles.redrawall,'Callback',{@redrawall_callback,obj})
            set(obj.guihandles.clearall,'Callback',{@clearall_callback,obj})
        end
    end
end

function make_siteexplorer(data,b,obj)

if obj.getGlobalSetting('SE_autosavecheck')
    if ~obj.getPar('autosavecheck')
        answer=questdlg('Set autosave settings in Parameters... Autosave is disabled. Enable autosave now?');
        if strcmp(answer,'Yes')
            obj.setPar('autosavecheck',1);
        end
    end
end

SEpreview=obj.SE.processors.preview;
    if isempty(SEpreview)||~isvalid(SEpreview.handle)
        obj.SE.processors.SEMainGui.make_siteexplorer;
        SEpreview=obj.SE.processors.preview;
    end
        set(SEpreview.handle,'Visible','on')
        figure(SEpreview.handle)
end

function redrawall_callback(a,b,obj)
% obj.SEpreview.redrawall;
obj.SE.processors.preview.redrawall

end

function clearall_callback(a,b,obj) 
answer=questdlg('Do you really want to clear all sites and cells?','Clear ROIs','No');
if strcmp(answer,'Yes')
obj.SE.processors.preview.clearall;
end
end

function setcurrent(a,b,obj,what)
imax=obj.SE.(['current' what]).image.imax;
obj.setPar(['se_imax_' what],imax(:));
obj.guihandles.(['se_imax_' what]).String=num2str(imax(:)',2);
end

function pard=guidef(obj)

pard.text3.object=struct('String','FoV (nm)','Style','text');
pard.text3.position=[1,2.25];
pard.text3.Width=.75;


pard.text5.object=struct('String','ROI (nm)','Style','text');
pard.text5.position=[1,3.0];
pard.text5.Width=.75;

pard.textdz.object=struct('String','dz (nm)','Style','text');
pard.textdz.position=[1,3.75];
pard.textdz.Width=.75;

pard.text4.object=struct('String','pixelsize (nm)','Style','text');
pard.text4.position=[1,1.5];
pard.text4.Width=.75;

pard.text1.object=struct('String','Site','Style','text');
pard.text1.position=[2,1];
pard.text1.Width=0.5;
pard.text2.object=struct('String','Cell','Style','text');
pard.text2.position=[3,1];
pard.text2.Width=0.5;
% pard.text6.object=struct('String','File','Style','text');
% pard.text6.position=[4,1];

pard.se_sitefov.object=struct('Style','edit','String',500); 
pard.se_sitefov.position=[2,2.25];
pard.se_sitefov.Width=0.75;
pard.se_cellfov.object=struct('Style','edit','String',5000); 
pard.se_cellfov.position=[3,2.25];
pard.se_cellfov.Width=0.75;
pard.se_sitepixelsize.object=struct('Style','edit','String',3); 
pard.se_sitepixelsize.position=[2,1.5];
pard.se_sitepixelsize.Width=0.75;
pard.se_cellpixelsize.object=struct('Style','edit','String',10); 
pard.se_cellpixelsize.position=[3,1.5];
pard.se_cellpixelsize.Width=0.75;
% pard.filepixelsize.object=struct('Style','edit','String',1); 
% pard.filepixelsize.position=[4,3];

pard.se_siteroi.object=struct('Style','edit','String',300); 
pard.se_siteroi.position=[2,3.0];
pard.se_siteroi.Width=0.75;

pard.se_dz.object=struct('Style','edit','String',1000); 
pard.se_dz.position=[2,3.75];
pard.se_dz.Width=0.75;

pard.se_imaxcheck_site.object=struct('Style','checkbox','String','Set Imax for sites to: [ch1 ch2 ...]','Value',0);
pard.se_imaxcheck_site.position=[5,1];
pard.se_imaxcheck_site.Width=2;

pard.se_imax_site.object=struct('Style','edit','String','1','Value',1);
pard.se_imax_site.position=[5,3];

pard.se_imax_site_current.object=struct('Style','pushbutton','String','<- Current','Callback',{{@setcurrent,obj,'site'}});
pard.se_imax_site_current.position=[5,4];

pard.se_imaxcheck_cell.object=struct('Style','checkbox','String','Set Imax for cells to: [ch1 ch2 ...]','Value',0);
pard.se_imaxcheck_cell.position=[6,1];
pard.se_imaxcheck_cell.Width=2;

pard.se_imax_cell.object=struct('Style','edit','String','1','Value',1);
pard.se_imax_cell.position=[6,3];

pard.se_imax_cell_current.object=struct('Style','pushbutton','String','<- Current','Callback',{{@setcurrent,obj,'cell'}});
pard.se_imax_cell_current.position=[6,4];


pard.se_rotate.object=struct('Style','checkbox','String','rotate','Value',0);
pard.se_rotate.position=[8,1];

pard.se_drawboxes.object=struct('Style','checkbox','String','draw boxes','Value',1);
pard.se_drawboxes.position=[8,2];

pard.se_drawsideview.object=struct('Style','checkbox','String','draw sideview','Value',1);
pard.se_drawsideview.position=[8,3];



pard.redrawall.object=struct('Style','pushbutton','String','redraw all','Value',0);
pard.redrawall.position=[10.5,4];
pard.redrawall.TooltipString='redraw all sites (inclusive files and cells) or selected sites if more than one site seltected';

pard.clearall.object=struct('Style','pushbutton','String','clear all','Value',0);
pard.clearall.position=[4,4];

pard.showSE.object=struct('Style','pushbutton','String','show ROI manager','Value',0);
pard.showSE.position=[11,1];
pard.showSE.Height=2;
pard.showSE.Width=1.5;

pard.outputParameters={'se_sitefov','se_cellfov','se_sitepixelsize','se_cellpixelsize','se_siteroi','se_drawboxes','se_rotate','se_imax_site','se_imaxcheck_cell','se_imax_cell','se_imaxcheck_site','se_dz','se_drawsideview'};
end