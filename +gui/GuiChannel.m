classdef GuiChannel< interfaces.LayerInterface
    properties
        shift
        colorrange
        quantilestore=-4;
        imax
        synchronizedFields
        rec_addpar = {'mingaussnm','mingausspix','gaussfac','par.gamma'};
        rec_addparval;
        name
    end
    methods
        function obj=GuiChannel(varargin)
            obj@interfaces.LayerInterface(varargin{:});
            obj.processorgui=false;
            obj.guiPar.par=[];
            obj.guiPar.mincall=[];
            obj.guiPar.maxcall=[];
            obj.guiPar.srmodes={'normal','z','field'};
            obj.outputParameters={'ch_filelist','filelistfilter','channels','layercheck','rendermode','render_colormode','renderfield',...
                'groupcheck','imaxtoggle','imax','lut','remout','shift','shiftxy_min','shiftxy_max','layer','colorrange',...
                'znm_min','znm_max','frame_min','frame_max','scalex','scaley'};
            obj.guiselector.show=true;
            obj.propertiesToSave={'rec_addparval'};
            
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        
        function set.rec_addparval(obj,par)
            
            obj.rec_addparval=par;
            if ~isempty(par)
            layerp=obj.getPar(obj.layerprefix);
            layerp=copyfields(layerp,par);
            obj.setPar(obj.layerprefix,layerp) 
            end
        end

        function makeGui(obj)
            makeGui@interfaces.DialogProcessor(obj);
            hmain=obj.getPar('filterpanel');
            units=hmain.Units;
            hmain.Units='pixels';
            pp=hmain.Position;
            hmain.Units=units;
            
            %filtertable gui
            hfigf=uipanel('Parent',hmain,'Units','pixel','Position',[0 200 pp(3) 130]);
            filtergui=gui.GuiFilterTable(hfigf,obj.P);
            filtergui.layer=obj.layer;
            filtergui.attachLocData(obj.locData);
            filtergui.makeGui;
            obj.children.filterTableGui=filtergui;

            %histogram gui
            hfigh=uipanel('Parent',hmain,'Units','pixel','Position',[0 0 pp(3) 200]);

            histgui=gui.GuiHistSlider(hfigh,obj.P);
            histgui.layer=obj.layer;
            histgui.attachLocData(obj.locData);
            histgui.makeGui;
            
            obj.children.histgui=histgui;

               %menu to detach Format Gui
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f); 
            hfigf.UIContextMenu=c;
            hfigh.UIContextMenu=c;
            
            m4 = uimenu(c,'Label','detach','Callback',{@detach_callback,obj,hmain});           
        end
        
        function initGui(obj)
            h=obj.guihandles;
            fn=fieldnames(h);
            for k=1:length(fn)
                h.(fn{k}).Callback={@parchanged_callback,obj,fn{k}};
            end
            
            callobj=obj;
            h.ch_filelist.Callback={@callobj.filelist_callback,1};

            h.rendermode.Callback={@rendermode_callback,obj};
            h.imaxtoggle.Callback={@imaxtoggle_callback,obj};
            h.render_colormode.Callback={@render_colormode_callback,obj}; 
            h.renderfield.Callback={@renderfield_callback,obj};  
            h.default_button.Callback=@obj.default_callback;
            h.defaultsave_button.Callback=@obj.default_callback;
            h.externalrender.Callback={@renderpar_callback,obj};
            h.copyalllayers_button.Callback={@copyalllayer_callback,obj};
            
            guimodules=obj.getPar('menu_plugins');
            fn=fieldnames(guimodules.Analyze.render);
            h.externalrender.String=[{'Select'},fn];
            
            cfields={ 'colorfield','String', 'locprecnm','znm','PSFxnm','locprecznm','frame','shiftxy'};
            obj.synchronizedFields={ 'colorfield', 'locprecnm','znm','PSFxnm','locprecznm','frame','shiftxy'};
            for k=1:length(cfields)
                 h.([cfields{k} 'b']).Callback={@obj.updatefields_callback,cfields{k}};
                 h.([cfields{k} '_min']).Callback={@obj.updatefields_callback,cfields{k}};
                 h.([cfields{k} '_max']).Callback={@obj.updatefields_callback,cfields{k}};
            end
             h.imaxb.Callback={@obj.updatefields_callback,'imax'};
             h.imax_min.Callback={@obj.updatefields_callback,'imax'};
%                  h.([cfields{k} '_max']).Callback={@obj.updatefields_callback,cfields{k}};
            
            obj.addSynchronization([obj.layerprefix 'layercheck'],h.layercheck,'Value',{@callobj.updatelayer_callback,'layercheck'})   
            h.layercheck.Callback={@updatelayercheck,obj};
            h.parbutton.Callback={@renderpar_callback,obj};
            
%             cross-layer synch
             callobj=obj;
            obj.addSynchronization('filelist_short',obj.guihandles.ch_filelist,'String',{@callobj.filelist_callback,1})
            obj.addSynchronization('filenumber',[],[],{@callobj.filenumber_callback})
            
            h.filelistfilter.Callback=@obj.filelistfilter_callback;
            obj.addSynchronization([obj.layerprefix 'filelistfilter'],h.filelistfilter,'Value',{@callobj.filelistfilter_callback})
           
           
            obj.addSynchronization([obj.layerprefix 'selectedField'],[],[],{@callobj.selectedField_callback})
            obj.addSynchronization('locFields',[],[],{@callobj.setcolorfield_callback})
            obj.addSynchronization([obj.layerprefix 'groupcheck'],h.groupcheck,'String')
            obj.addSynchronization([obj.layerprefix 'channels'],h.channels,'String')
            obj.addSynchronization([obj.layerprefix 'shiftxy_min'],h.shiftxy_min,'String',{@callobj.updatefields_callback,'shiftxy'})
            obj.addSynchronization([obj.layerprefix 'shiftxy_max'],h.shiftxy_max,'String',{@callobj.updatefields_callback,'shiftxy'})
            obj.addSynchronization(['frame_min'],[],[],{@callobj.updateframes_callback})
            obj.addSynchronization(['frame_max'],[],[],{@callobj.updateframes_callback})
%             obj.addSynchronization([obj.layerprefix 'scalex'],h.scalex,'String')
%             obj.addSynchronization([obj.layerprefix 'scaley'],h.scaley,'String')
            obj.guihandles=h;

            recpar=renderpardialog(obj.rec_addparval,1);
            p=obj.getAllParameters;
            layerp=copyfields(p,recpar);
            obj.setPar(obj.layerprefix,layerp);
            
            setvisibility(0,0,obj)
            obj.updateLayerField;

            obj.updatelayer_callback;
            updatelayercheck(h.layercheck,0,obj);
            obj.setfiltergray;
        end
        function filelistfilter_callback(obj,a,b)
            sf={'filenumber',[],[],obj.getSingleGuiParameter('filelistfilter')};
            obj.setPar('selectedField',sf,'layer',obj.layer);
        end
        function updateframes_callback(obj,a,b)
            p.frame_min=obj.getPar('frame_min');
            p.frame_max=obj.getPar('frame_max');
            obj.setGuiParameters(p);
        end
        function filenumber_callback(obj,a,b,c)            
            fn=obj.getPar('filenumber');
            obj.guihandles.ch_filelist.Value=fn;
            obj.filelist_callback;
            notify(obj.P,'sr_render');
        end
           
        function filelist_callback(obj,x,y,z)
                znm=obj.locData.getloc('znm');
                if isfield(znm,'znm')&&any(znm.znm~=0) % loical: dont update hist slider
                     obj.updatefields_callback(0,0,'znm',[],false)
                     obj.updatefields_callback(0,0,'locprecznm',[],false)
                     obj.updatefields_callback(0,0,'PSFxnm',false,false)
                else
                    obj.updatefields_callback(0,0,'PSFxnm',[],false)
                    obj.updatefields_callback(0,0,'znm',false,false)
                    obj.updatefields_callback(0,0,'locprecznm',[],false)
                end                   
                obj.updatefields_callback(0,0,'locprecnm',true,true)
                
                resetframe=obj.getGlobalSetting('resetframefilter');
                if ~isempty(resetframe) && resetframe
                    obj.updatefields_callback(0,0,'frame',false,false)
                end
%              sf={'filenumber',[],[],true};
%              obj.setPar('selectedField',sf,'layer',obj.layer);
             obj.setGuiParameters(struct('filelistfilter',true));
             obj.updateLayerField;
             setvisibility(0,0,obj)
             obj.setfiltergray;
        end
        
        function setcolorfield_callback(obj)
            if isempty(obj.locData.loc), return; end
            if isfield(obj.locData.loc,'colorfield') && length(obj.locData.loc.colorfield)==length(obj.locData.loc.frame)
                return
            end
            numlocs=length(obj.locData.loc.frame);
            
            obj.locData.loc.colorfield=single((1:numlocs)'/numlocs*2-0.5);
            if ~isempty(obj.locData.grouploc)
                numlocs=length(obj.locData.grouploc.frame);
                obj.locData.grouploc.colorfield=single((1:numlocs)'/numlocs);
            end
        end
         
        function updateLayerField(obj,field,value)
            if nargin<2 %update all
                layerp=obj.getPar(obj.layerprefix);
            	p=obj.getAllParameters;
                layerp=copyfields(layerp,p,fieldnames(layerp));
                obj.setPar(obj.layerprefix,layerp);
                fn=obj.guihandles.ch_filelist.Value;
                sf={'filenumber',fn,fn,p.filelistfilter};
                obj.setPar('selectedField',sf,'layer',obj.layer);
            else             
                if nargin<3
                    value=obj.getSingleGuiParameter(field);
                end
                %reads field from gui, writes it in P.par.layerX_.(field)
                try
                p=obj.getPar(obj.layerprefix);
                p.(field)=value;
                obj.setPar(obj.layerprefix,p);
                catch err
                end
            end
        end
        
        function updatelayer_callback(obj,a,b,field)
            layeron=obj.getPar([obj.layerprefix 'layercheck']);
            layerp=obj.getPar(obj.layerprefix);
            layerp.layercheck=layeron;
            recpar=renderpardialog(obj.rec_addparval,1);
%             p=obj.getAllParameters;
            layerp=copyfields(layerp,recpar);
            obj.setPar(obj.layerprefix,layerp);
        end
        
        function setlayer(obj,layer)
            if obj.layer==layer
                state='on';            
            else
                state='off';
            end
            obj.children.filterTableGui.handle.Visible=state;
            obj.children.histgui.handle.Visible=state;
        end
        
        function selectedField_callback(obj)
            sfield=obj.getPar('selectedField','layer',obj.layer);
            field=sfield{1};
            
            minf=num2str(sfield{2});
            maxf=num2str(sfield{3});
            p.([field '_min'])=minf;
            p.([field '_max'])=maxf;
            obj.setGuiParameters(p);
            
            obj.updateLayerField([field '_min'],sfield{2})
            obj.updateLayerField([field '_max'],sfield{3})
            
            allfields={'PSFxnm','znm','locprecznm','locprecnm','frame'};
            if any(strcmp(allfields,field))
                if ~isempty(sfield{4})
                    if sfield{4}(1)
                        color=[1 1 1]*.94;
                    else
                        color=[1 1 1]*.7;
                    end
                    obj.guihandles.([field '_min']).BackgroundColor=color;
                    obj.guihandles.([field '_max']).BackgroundColor=color;
                end
            end
            if strcmp(field,'filenumber')
                if ~isempty(sfield{4})
                    obj.setGuiParameters(struct('filelistfilter',sfield{4}(1)))
%                     if sfield{4} %do filtering
%                         filenumber=min(max(1,round(sfield{2})),length(obj.guihandles.ch_filelist.String));
%                         obj.guihandles.ch_filelist.Value=filenumber;
%                     end
                end
            end
        end
        function setfiltergray(obj)
            cf={'locprecnm','znm','PSFxnm','locprecznm','frame'};
            s=obj.getPar('filtertable','layer',obj.layer);
            if isempty(s)
                return
            end
            
            fields=s(:,1);
            for k=1:length(cf)
                fieldh=cf{k};
                ind=find(strcmp(fields,fieldh));
                if ~isempty(ind)
                    filterstate=s{ind(1),7};
                    if filterstate
                        color=[1 1 1];
                    else
                        color=[1 1 1]*.7;
                    end
                    obj.guihandles.([fieldh '_min']).BackgroundColor=color;
                    obj.guihandles.([fieldh '_max']).BackgroundColor=color;
                end
            end
        end
        function updatefields_callback(obj,f1, data, field,filteron,updatehist)
            if nargin <6
                updatehist=true;
            end
              
            if nargin<4||isempty(field)
                field=f1;
                filteron=[];
            
            elseif nargin<5
                filteron=true;
            end
            obj.updateLayerField([field '_min'])
            obj.updateLayerField([field '_max'])
            if strcmp(field,'colorfield');
                filteron=false;
            end
            sf={field,obj.getSingleGuiParameter([field '_min']),obj.getSingleGuiParameter([field '_max']),filteron,updatehist};
            obj.setPar('selectedField',sf,'layer',obj.layer)
            obj.selectedField_callback;
        end
        
        function default_callback(obj,callobj,b)
            deffile=[obj.getPar('SettingsDirectory') filesep 'temp' filesep 'Channel_default.mat'];
            fh=getParentFigure(obj.handle);
            modifiers = get(fh,'currentModifier');
            if ismember('shift',modifiers)||strcmpi(callobj.String,'save');
                %save
                defaultChannelParameters=getChannelParameters(obj);
                
%                 sdfdy
%                 defaultChannelParameters=obj.getGuiParameters;
%                 globalParameters=obj.getLayerParameters(obj.layer);
%                 globalParameters=myrmfield(globalParameters,{'guimodules','mainGuihandle','mainGui','filterpanel','ov_axes','guiFormat',...
%                     'sr_image','sr_imagehandle','sr_figurehandle','sr_axes','ROI_lineannotation_handle_1','ROI_lineannotation_handle_2'});
% %                 save(deffile,'defaultChannelParameters','globalParameters');
                save(deffile,'defaultChannelParameters');
                disp('default parameters saved.');
                obj.status('default parameters saved.');
            else
                if exist(deffile,'file')
                    l=load(deffile);
                    p=l.defaultChannelParameters;
                    
                    p=rmfield(p,{'ch_filelist','render_colormode','renderfield'});
                    setChannelParameters(obj,p);
                    
%                     obj.setGuiParameters(p);
                    obj.status('default parameters restored.');
                    %update
%                     obj.updateLayerField; 
                    
                    %update filter table
                    %filter again
                    
                else
                    disp('no default configuration found at: /settings/temp/Channel_default.mat')
                end
            end
        end
        
        function setGuiParameters(obj,varargin)          
            setGuiParameters@interfaces.GuiModuleInterface(obj,varargin{:});  
            if length(varargin)>1&&varargin{2}==true
%             %update fields to histogram
                obj.updateLayerField;
                
            end
            if isfield(varargin{1},'rec_addparval')
                    obj.rec_addparval=varargin{1}.rec_addparval;
            end
%             obj.updatelayer_callback;
        end
        
        function setvisibility(obj,field)
            if nargin==1
                setvisibility(0,0,obj);
            else
                setvisibility(0,0,obj,field);
            end
        end

    end
end


function p=getChannelParameters(obj)
p=obj.getGuiParameters;
p2=obj.getLayerParameters(obj.layer);
p=copyfields(p,p2,obj.rec_addpar);
end

function setChannelParameters(obj,p)
obj.setGuiParameters(p);
obj.updateLayerField;
%copy additional paramters
layerp=obj.getPar(obj.layerprefix);
layerp=copyfields(layerp,p);
obj.setPar(obj.layerprefix,layerp,obj.rec_addpar) 
guirender=obj.getPar('mainGui').children.guiRender.children;
thisgui=guirender.(['Layer' num2str(obj.layer)]);

thisfield=thisgui.children.histgui.getGuiParameters.field;
if isempty(thisfield)
    syncf=obj.synchronizedFields;
else
syncf=setdiff(obj.synchronizedFields,thisfield);
updatefields_callback(obj,0, 0, thisfield,[])
end
for k=1:length(syncf)
    updatefields_callback(obj,0, 0, syncf{k},[])
%     obj.setPar(obj.synchronizedFields{k})
end


%refilter?
end
function updatelayercheck(object,b,obj)

    lo=obj.getPar('sr_layerson');
    lo(obj.layer)=object.Value;
    obj.setPar('sr_layerson',lo);
    obj.setPar([obj.layerprefix 'layercheck'],object.Value);
    if obj.layer>6
        notify(obj.P,'sr_render');
    end
% end
end
function imaxtoggle_callback(object,data,obj)
if(object.Value)
    object.String='quantile';
    obj.guihandles.imax_min.String=num2str(obj.quantilestore);
else
    object.String='Imax';
    obj.quantilestore=str2num(obj.guihandles.imax_min.String);
    try
    imax=obj.locData.layer(obj.layer).images.finalImages.imax;
    catch
        imax=1;
    end
    obj.guihandles.imax_min.String=num2str(imax);
end
obj.updateLayerField('imaxtoggle');
obj.updateLayerField('imax_min');
end

function copyalllayer_callback(object,data,obj)

excludefields={'ch_filelist','channels','layercheck'};
phere=getChannelParameters(obj);
phere=rmfield(phere,excludefields);
thislayer=obj.layer;
numlayers=obj.getPar('numberOfLayers');
mainGui=obj.getPar('mainGui');
guirender=mainGui.children.guiRender.children;
% thisgui=guirender.(['Layer' num2str(thislayer)]);
% 
% thisfield=thisgui.children.histgui.getGuiParameters.field;
% filtertable=thisgui.children.filterTableGui.getGuiParameters;
% 
% fn=fieldnames(obj.P.par);
% indlh=(strfind(fn,['layer' num2str(thislayer) '_'],'ForceCellOutput',true));
% indfield=~cellfun(@isempty,indlh);
% layerhfields=fn(indfield);
% 
% filterh=obj.locData.layer(thislayer).filter;
% gfilterh=obj.locData.layer(thislayer).groupfilter;

for k=1:numlayers
    if k==thislayer
        continue
    end
%     %copy GUI
    guilayer=guirender.(['Layer' num2str(k)]);
    setChannelParameters(guilayer,phere);
%     guilayer.setGuiParameters(phere);
%     
%     layergui=guirender.(['Layer' num2str(k)]);
%     %copy gui parameters
%     for l=1:length(layerhfields)
%         newfield=strrep(layerhfields{l},['layer' num2str(thislayer) '_'],['layer' num2str(k) '_']);
%         pp=obj.getPar(layerhfields{l});
%         obj.setPar(newfield,pp);
% %         obj.P.par.(newfield).content=obj.P.par.(layerhfields{l}).content; %keep callbacks etc 
%     end
%     
%     %replace by guilayer.updateLayerField??
%     
%     %update filtertable
%     
%     layergui.children.filterTableGui.setGuiParameters(filtertable);
%     %update histslider
% %     layergui.children.histgui.selectedField_callback(thisfield);
%     
%     %copy filter
%     obj.locData.layer(k).filter=filterh;
%     obj.locData.layer(k).groupfilter=gfilterh;
    

end

end

function render_colormode_callback(object,data,obj)
p=obj.getAllParameters;
switch p.render_colormode.selection
    case {'normal'}
        cmin=0;
        cmax=1;
        obj.setcolorfield_callback
    case 'z'
        if isfield(obj.locData.loc,'znm')
        obj.locData.loc.colorfield=single(obj.locData.loc.znm);
        obj.locData.grouploc.colorfield=single(obj.locData.grouploc.znm);
        cmin=-300;
        cmax=300;
        else 
            disp('no z information')
            return
        end
    case 'field'
        renderfield_callback(object, 0,obj)
        return
%         if ~isempty(obj.colorrange)&&length(obj.colorrange.mincall)>=p.renderfield.Value
%             cmin=obj.colorrange.mincall(p.renderfield.Value);
%             cmax=obj.colorrange.maxcall(p.renderfield.Value);
%         else
%             cmin=0;
%             cmax=1;
%         end
    otherwise %tiff file
            cmin=0;
            cmax=1;
end
% obj.setPar([obj.layerprefix 'colorfield_min'],num2str(cmin),'String')
% obj.setPar([obj.layerprefix 'colorfield_max'],num2str(cmax),'String')
if obj.getSingleGuiParameter('colorauto')
obj.guihandles.colorfield_min.String=num2str(cmin);
obj.guihandles.colorfield_max.String=num2str(cmax);
end

obj.updateLayerField('render_colormode');
obj.updateLayerField('colorfield_min');
obj.updateLayerField('colorfield_max');
 setvisibility(0,0,obj)
end

function renderfield_callback(object, handle,obj)
if isempty(obj.locData.loc)
     setvisibility(0,0,obj)
    return
end
p=obj.getSingleGuiParameter('renderfield');
field=p.selection;
v=obj.locData.getloc(field,'layer',obj.layer).(field);
if ~isempty(v) && obj.getSingleGuiParameter('colorauto')
obj.locData.loc.colorfield=single(obj.locData.loc.(field));
obj.locData.grouploc.colorfield=single(obj.locData.grouploc.(field));
q=myquantilefast(v,[0.02,0.98]);
dx=10^floor(log10(abs(q(2))/100));
minv=round(q(1)/dx)*dx;
maxv=round(q(2)/dx)*dx;            
obj.colorrange.mincall(p.Value)=minv-dx;
obj.colorrange.maxcall(p.Value)=maxv+dx;
obj.guihandles.colorfield_min.String=num2str(minv);
obj.guihandles.colorfield_max.String=num2str(maxv);
% render_colormode_callback(object,0,obj)      
end
obj.updateLayerField('render_colormode');
obj.updateLayerField('renderfield');
obj.updateLayerField('colorfield_min');
obj.updateLayerField('colorfield_max');
 setvisibility(0,0,obj)
end

function renderpar_callback(object,data,obj)
p=obj.getGuiParameters;

switch p.rendermode.selection
case 'Other'
    if strcmp(p.externalrender.selection,'Select')
        return
    end
    rendermodules=obj.getPar('rendermodules');
    
    if length(rendermodules)>=obj.layer && ~isempty(rendermodules{obj.layer})&&isvalid(rendermodules{obj.layer})
        delete(rendermodules{obj.layer}.handle)
        delete(rendermodules{obj.layer})
    end
    renderer=p.externalrender.selection;
    guiplugins=obj.getPar('menu_plugins');
    pluginpath=guiplugins.Analyze.render.(renderer).module;
    module=plugin(pluginpath{1:3});
    if strcmp(object.Style,'popupmenu')
        vis='off';
    else
        vis='on';
    end
    module.handle=figure('MenuBar','none','Toolbar','none','Name',pluginpath{end},'Visible',vis);
    module.attachPar(obj.P);
    module.attachLocData(obj.locData);
    ph.Vrim=100;
    ph.Xrim=10;
    module.setGuiAppearence(ph)
    module.makeGui;
    module.guihandles.showresults.Value=0;
    module.layer=obj.layer;
    rendermodules{obj.layer}=module;
    obj.setPar('rendermodules',rendermodules);
    
otherwise
     layerp=obj.getPar(obj.layerprefix);
     recp=obj.rec_addparval;
    [par,settings]=renderpardialog(recp);
    if ~isempty(par)
        obj.rec_addparval=par;
       layerp=copyfields(layerp,par);
       obj.setPar(obj.layerprefix,layerp)       
%        if par.copyall
%            nl=obj.getPar('numberOfLayers');
%            for k=1:nl
%                layerp=obj.getPar('','layer',k);
%                layerp=copyfields(layerp,settings);
%                obj.setPar('',layerp,'layer',k)
%            end
%        end
    end
end
end

function rendermode_callback(a,b,obj)
setvisibility(a,b,obj,'rendermode');
render_colormode_callback(a,b,obj)
end

function setvisibility(a,b,obj,field)
hgui=obj.guihandles;
if hgui.render_colormode.Value>length(hgui.render_colormode.String)
    set(hgui.render_colormode,'Value',1);
%     disp('gui channel line 370')
end
p=obj.getAllParameters;


% z data?
fh={'PSFxnmb','PSFxnm_min','PSFxnm_max','znmb','znm_min','znm_max','locprecznmb','locprecznm_min','locprecznm_max',...
    'channels','text1','renderfield','groupcheck','frameb','frame_min','frame_max'...
    'locprecnmb','locprecnm_min','locprecnm_max','intensitycoding','intensitytxt','colortxt','renderparameter'};
%tiff
if strcmp(p.rendermode.selection,'tiff')||strcmp(p.rendermode.selection,'raw') %obj.fileinfo(fileselect).istiff
    obj.fieldvisibility('off',fh);
%     controlVisibility(hgui,fh,'off')
else
    obj.fieldvisibility('on',fh);
%     controlVisibility(hgui,fh,'on')
    zh={'znmb','znm_min','znm_max','locprecznmb','locprecznm_min','locprecznm_max','shiftxy_z'};
    znoth={'PSFxnmb','PSFxnm_min','PSFxnm_max'};
    znm=obj.locData.getloc({'znm'},'layer',obj.layer,'position','all');
    
    if ~isempty(znm)&&~isempty(znm.znm)&&any(znm.znm~=0)
        obj.fieldvisibility('on',zh);
%         controlVisibility(hgui,zh,'on')
%         controlVisibility(hgui,znoth,'off')
        obj.fieldvisibility('off',znoth);
    else
        obj.fieldvisibility('off',zh);
        obj.fieldvisibility('on',znoth);
%         controlVisibility(hgui,zh,'off')
%         controlVisibility(hgui,znoth,'on')
    end
end

if strcmp(p.rendermode.selection,'Other')
    obj.fieldvisibility('on','externalrender');
%     controlVisibility(hgui,'externalrender','on')
else
    obj.fieldvisibility('off','externalrender');
%     controlVisibility(hgui,'externalrender','off')
end

%render_colormode 
hgui.render_colormode.Value=min(length(hgui.render_colormode.String),hgui.render_colormode.Value);
switch lower(p.render_colormode.selection)
case 'normal'
    obj.fieldvisibility('off','renderfield');
%     controlVisibility(hgui,{'renderfield'},'off')
case 'z'
    obj.fieldvisibility('on',{'remout','cb','colorfield_min','colorfield_max'});
    obj.fieldvisibility('off','renderfield');
%     controlVisibility(hgui,{'remout','cb','colorfield_min','colorfield_max'},'on')
%     controlVisibility(hgui,{'renderfield'},'off')
case 'field'
    obj.fieldvisibility('on',{'renderfield','remout','cb','colorfield_min','colorfield_max'});
%      controlVisibility(hgui,{'renderfield','remout','cb','colorfield_min','colorfield_max'},'on')
end

% render or image
switch lower(p.rendermode.selection)
case {'hist','gauss','dl','other','constgauss'}
    if isstruct(obj.locData.loc)
        set(hgui.renderfield,'String',fieldnames(obj.locData.loc));
        hgui.renderfield.Value=min(hgui.renderfield.Value,length(hgui.renderfield.String));
    end
    set(hgui.render_colormode,'String',{'normal','z','field'});
    if hgui.render_colormode.Value>length(hgui.render_colormode.String)
        set(hgui.render_colormode,'Value',1);
        disp('gui channel line 370')
    end
    hgui.render_colormode.TooltipString='normal: intenisty coded. z: z-coded. Field: select field for color-coding';
    
case {'tiff','raw'}
    if strcmp(p.rendermode.selection,'tiff');
        form='tif';
        hgui.render_colormode.TooltipString='Name of tiff file (loaded with File -> add)';
    else
        form='raw';
        hgui.render_colormode.TooltipString='raw camera frames: frame number';
    end
    s={};
    sind=1;
    file=hgui.ch_filelist.Value;

    s{1}='empty';
    if ~isempty(obj.locData.files.file)&&isfield(obj.locData.files.file(file),form)
        for k=1:length(obj.locData.files.file(file).(form))
            if strcmp(form,'tif')
               if ~isempty(obj.locData.files.file(file).tif(k).info)
                    [path,f,ext]=fileparts(obj.locData.files.file(file).tif(k).info.name);
                    s{sind}=f;
                    sind=sind+1;
               end
            else
                s{sind}=num2str(obj.locData.files.file(file).raw(k).frame);
                sind=sind+1;
            end
        end

        set(hgui.render_colormode,'String',s);
        set(hgui.render_colormode,'Value',min(hgui.render_colormode.Value,length(s)));
        set(hgui.colorfield_min,'String','0')
        set(hgui.colorfield_max,'String','1')
        updateLayerField(obj,'colorfield_min')
        updateLayerField(obj,'colorfield_max')
        obj.updateLayerField('render_colormode');
    end
end

if contains(p.rendermode.selection,'constGauss')
    hgui.renderparameter.Visible='on';
else
    hgui.renderparameter.Visible='off';
end


if nargin>3
    updateLayerField(obj,field)
end
end




function [paro,settings]=renderpardialog(par,default)
if nargin==0 || isempty(par) || ~isfield(par,'mingaussnm')
%     if isempty(obj.rec_addparval)
    par.mingaussnm=3;
    par.mingausspix=.7;
    par.gaussfac=0.4;
    par.gamma=1;
    par.copyall=false;
%     else
%         par=obj.rec_addparval;
%     end
end
if nargin<2||~default
[settings, button] = settingsdlg(...
    'Description', 'Render Parameters',... 
    'title'      , 'Par',...  
    {'min sigma Gauss (nm)';'mingaussnm'}, par.mingaussnm,...
    {'min sigma Gauss (pix)';'mingausspix'}, par.mingausspix,...
    {'factor Gauss';'gaussfac'}, par.gaussfac,...
    {'gamma image';'gamma'}, par.gamma); %,...
%  {'copy to all cahnnels';'copyall'},false);

if strcmpi(button,'ok')
    par=addFields(par,settings);
    paro=par;
else
    paro=[];
end
else
    paro=par;
end
paro.copyall=false;
end


function parchanged_callback(object,event,obj,field)
obj.updateLayerField(field)
end


function pard=guidef(obj)
pard.layercheck.object=struct('Style','checkbox','String','','Value',1);
pard.layercheck.position=[1,1];
pard.layercheck.Width=0.2;
pard.layercheck.TooltipString='switch layer on and off';
            
pard.ch_filelist.object=struct('Style','popupmenu','String',{'File'});
pard.ch_filelist.position=[1,1.2];
pard.ch_filelist.Width=1.9;
pard.ch_filelist.TooltipString='which file (loc or image) to display';

pard.filelistfilter.object=struct('Style','checkbox','String','','Value',1);
pard.filelistfilter.position=[1,3.1];
pard.filelistfilter.Width=0.2;
pard.filelistfilter.TooltipString='Filter files';
pard.filelistfilter.Optional=true;

pard.text1.object=struct('Style','text','String','Ch');
pard.text1.position=[1,3.4];
pard.text1.Width=0.2;

pard.channels.object=struct('Style','edit','String','0 1');
pard.channels.position=[1,3.6];
pard.channels.Width=0.6;
pard.channels.TooltipString='channels to display. Use a b c and a:c notation';

pard.groupcheck.object=struct('Style','checkbox','String','group','Value',1);
pard.groupcheck.position=[1,4.2];
pard.groupcheck.Width=.6;
pard.groupcheck.TooltipString='use grouped or ungrouped locs';

w1=.6;
w2=0.9;
w3=0.9;
p1=1;
p2=1.5;
p3=2.4;
pard.rendertxt.object=struct('Style','text','String','Renderer:');
pard.rendertxt.position=[2,p1];
pard.rendertxt.Width=w1;  
pard.rendertxt.TooltipString='how to render image. DL is diffraction limited';


pard.rendermode.object=struct('Style','popupmenu','String',{{'hist','Gauss','constGauss','DL','tiff','raw','Other'}},'Value',2);
pard.rendermode.position=[2,p2];
pard.rendermode.Width=w2;  
pard.rendermode.TooltipString='how to render image. DL is diffraction limited';

pard.renderparameter.object=struct('Style','edit','String','','Visible','off');
pard.renderparameter.position=[2,p3];
pard.renderparameter.Width=w3;
pard.renderparameter.TooltipString='numeric parameter for renderer. constGauss: size of Gauss (nm). Leave empty for automatic size determination';
% pard.renderparameter.Optional=true;

pard.externalrender.object=struct('Style','popupmenu','String','empty');
pard.externalrender.position=[2,p3];
pard.externalrender.Width=w3;
pard.externalrender.TooltipString='External renderer';

pard.intensitytxt.object=struct('Style','text','String','Intensity:');
pard.intensitytxt.position=[3,p1];
pard.intensitytxt.Width=w1;  
pard.intensitytxt.TooltipString='How to normalize every localization: normal: Integral=1, photons: total photons, blinks: number of connected localizations';
pard.intensitytxt.Optional=true;

pard.intensitycoding.object=struct('Style','popupmenu','String',{{'normal','photons','blinks','√photons','√blinks'}});
pard.intensitycoding.position=[3,p2];
pard.intensitycoding.Width=w2;
pard.intensitycoding.TooltipString=pard.intensitytxt.TooltipString;
pard.intensitycoding.Optional=true;

pard.parbutton.object=struct('Style','pushbutton','String','render par');
pard.parbutton.position=[3,p3];
pard.parbutton.Width=w3;
pard.parbutton.TooltipString='Additional render paramters';
pard.parbutton.Optional=true;


pard.tiftxt.object=struct('Style','text','String','Image:');
pard.tiftxt.position=[4,p1];
pard.tiftxt.Width=w1;  

pard.colortxt.object=struct('Style','text','String','Colormode:');
pard.colortxt.position=[4,p1];
pard.colortxt.Width=w1*1.1;  


pard.render_colormode.object=struct('Style','popupmenu','String',{obj.guiPar.srmodes}); 
pard.render_colormode.position=[4,p2];
pard.render_colormode.Width=w2;
pard.render_colormode.TooltipString=sprintf('normal: intensity coded, z: z color coded, \n param: select field which to color code');

pard.colortxt.TooltipString=pard.render_colormode.TooltipString;

pard.renderfield.object=struct('Style','popupmenu','String','field');            
pard.renderfield.position=[4,p3];
pard.renderfield.Width=w3;
pard.renderfield.TooltipString='field to color code';



pard.luttxt.object=struct('Style','text','String','LUT:');
pard.luttxt.position=[5,p1];
pard.luttxt.Width=w1;  
pard.luttxt.TooltipString='how to render image. DL is diffraction limited';

pard.lut.object=struct('Style','popupmenu','String',{mymakelut});
pard.lut.position=[5,p2];
pard.lut.Width=w2;
pard.lut.TooltipString='select the lookup table';

pard.remout.object=struct('Style','checkbox','String','restrict');
pard.remout.position=[5,p3];
pard.remout.Width=w3*.7;
pard.remout.TooltipString=sprintf('if checked: remove loclizations outside lut. \n If unchecked: set them to maximum color');
pard.remout.Optional=true;

pard.lutinv.object=struct('Style','checkbox','String','Inv','Value',0);
pard.lutinv.position=[5,p3+w3*.6];
pard.lutinv.Width=w3*.5;
pard.lutinv.TooltipString=sprintf('if checked: invert the colormap');
pard.lutinv.Optional=true;


pard.colorfieldb.object=struct('Style','pushbutton','String','c range');
pard.colorfieldb.position=[6,p1];
pard.colorfieldb.Width=.6;
pard.colorfieldb.TooltipString=sprintf('Range of values to fill the lookup table (LUT). Usually 0 and 1. \n If fields (or z) are used for color coding: range of these values mapped to LUT.');

pard.colorfieldb.Optional=true;

pard.colorfield_min.object=struct('Style','edit','String','0');
pard.colorfield_min.position=[6,1.6];
pard.colorfield_min.Width=.5;
pard.colorfield_min.TooltipString=pard.colorfieldb.TooltipString;
pard.colorfield_min.Optional=true;

pard.colorfield_max.object=struct('Style','edit','String','1');
pard.colorfield_max.position=[6,2.1];
pard.colorfield_max.Width=.5;
pard.colorfield_max.TooltipString=pard.colorfieldb.TooltipString;
pard.colorfield_max.Optional=true;

pard.colorauto.object=struct('Style','checkbox','String','auto','Value',1);            
pard.colorauto.position=[6,2.6];
pard.colorauto.Width=0.6;
pard.colorauto.TooltipString='automatically set color range';
pard.colorauto.Optional=true;

pard.imaxb.object=struct('Style','pushbutton','String','contrast');
pard.imaxb.position=[7,1];
pard.imaxb.Width=.6;
pard.imaxb.TooltipString=sprintf('Contrast');


pard.imaxtoggle.object=struct('Style','togglebutton','String','quantile','Value',1);
pard.imaxtoggle.position=[7,1.6];
pard.imaxtoggle.Width=.5;
pard.imaxtoggle.TooltipString='toggle absolute intensity maximum (Imax) or quantile';
pard.imaxtoggle.Optional=true;

pard.imax_min.object=struct('Style','edit','String','-3.5');
pard.imax_min.position=[7,2.1];
pard.imax_min.Width=.5;
pard.imax_min.TooltipString=sprintf('absolut intensity or \n quantile of the pixels that are not saturated (0<q<1, typically 0.999) or  \n v for q=1-10^(v), v<0, typically -3.5');

pard.imax_max.object=struct('Style','edit','String','-3.5','Visible','off');
pard.imax_max.position=[7,2.2];

p4=3.4;
pard.filtertxt.object=struct('Style','text','String','Filter fields:');
pard.filtertxt.position=[2,p4];
pard.filtertxt.Width=2;  
pard.filtertxt.TooltipString='';

    
pard.locprecnmb.object=struct('Style','pushbutton','String','locprec');
pard.locprecnmb.position=[3,3.4];
pard.locprecnmb.Width=.6;


pard.locprecnm_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.locprecnm_min.position=[3,4];
pard.locprecnm_min.Width=.5;

pard.locprecnm_max.object=struct('Style','edit','String','30');  
pard.locprecnm_max.position=[3,4.5];
pard.locprecnm_max.Width=.5;

pard.znmb.object=struct('Style','pushbutton','String','z');
pard.znmb.position=[4,3.4];
pard.znmb.Width=.6;

pard.znm_min.object=struct('Style','edit','String','-500','BackgroundColor',[1 1 1]*.7);
pard.znm_min.position=[4,4];
pard.znm_min.Width=.5;

pard.znm_max.object=struct('Style','edit','String','500');  
pard.znm_max.position=[4,4.5];
pard.znm_max.Width=.5;
      
pard.PSFxnmb.object=struct('Style','pushbutton','String','PSF xy');
pard.PSFxnmb.position=[4,3.4];
pard.PSFxnmb.Width=.6;
% pard.PSFxnmb.Optional=true;            
pard.PSFxnm_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.PSFxnm_min.position=[4,4];
pard.PSFxnm_min.Width=.5;
% pard.PSFxnm_min.Optional=true;
pard.PSFxnm_max.object=struct('Style','edit','String','175');  
pard.PSFxnm_max.position=[4,4.5];
pard.PSFxnm_max.Width=.5;
% pard.PSFxnm_max.Optional=true;

pard.locprecznmb.object=struct('Style','pushbutton','String','locprec z');
pard.locprecznmb.position=[5,3.4];
pard.locprecznmb.Width=.6;
pard.locprecznmb.Optional=true;

pard.locprecznm_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.locprecznm_min.position=[5,4];
pard.locprecznm_min.Width=.5;
pard.locprecznm_min.Optional=true;

pard.locprecznm_max.object=struct('Style','edit','String','45');  
pard.locprecznm_max.position=[5,4.5];
pard.locprecznm_max.Width=.5; 
pard.locprecznm_max.Optional=true;

pard.frameb.object=struct('Style','pushbutton','String','frame');
pard.frameb.position=[6,3.4];
pard.frameb.Width=.6;
pard.frameb.Optional=true;

pard.frame_min.object=struct('Style','edit','String','0','BackgroundColor',[1 1 1]*.7);
pard.frame_min.position=[6,4];
pard.frame_min.Width=.5;
pard.frame_min.Optional=true;

pard.frame_max.object=struct('Style','edit','String','inf');  
pard.frame_max.position=[6,4.5];
pard.frame_max.Width=.5;
pard.frame_max.Optional=true;

pard.shiftxyb.object=struct('Style','pushbutton','String','shift xyz');
pard.shiftxyb.position=[7,3.4];
pard.shiftxyb.Width=.6;
pard.shiftxyb.TooltipString='Shift reconstructed image by this value (nm). Useful to correct for chromatic aberrations.';
pard.shiftxyb.Optional=true;

pard.shiftxy_min.object=struct('Style','edit','String','0');
pard.shiftxy_min.position=[7,4];
pard.shiftxy_min.Width=1/3;
pard.shiftxy_min.TooltipString=pard.shiftxyb.TooltipString;
pard.shiftxy_min.Optional=true;
pard.shiftxy_max.object=struct('Style','edit','String','0');
pard.shiftxy_max.position=[7,4+1/3];
pard.shiftxy_max.Width=1/3;
pard.shiftxy_max.TooltipString=pard.shiftxyb.TooltipString;
pard.shiftxy_max.Optional=true;

pard.shiftxy_z.object=struct('Style','edit','String','0');
pard.shiftxy_z.position=[7,5-1/3];
pard.shiftxy_z.Width=1/3;
pard.shiftxy_z.TooltipString=pard.shiftxyb.TooltipString;
pard.shiftxy_z.Optional=true;



% pard.scalext.object=struct('Style','text','String','Scale X:');
% pard.scalext.position=[8,3.4];
% pard.scalext.Width=.6;
% pard.scalext.TooltipString='Scale image in X/Y by this factor. Use it for image deformations e.g. due to cylindrical lens.';
% pard.scalext.Optional=true;
% 
% pard.scalex.object=struct('Style','edit','String','1');
% pard.scalex.position=[8,4];
% pard.scalex.Width=1/3;
% pard.scalex.TooltipString=pard.scalext.TooltipString;
% pard.scalex.Optional=true;
% 
% pard.scaleyt.object=struct('Style','text','String','Y:');
% pard.scaleyt.position=[8,4+1/3];
% pard.scaleyt.Width=1/3;
% pard.scaleyt.TooltipString=pard.scalext.TooltipString;
% pard.scaleyt.Optional=true;
% 
% pard.scaley.object=struct('Style','edit','String','1');
% pard.scaley.position=[8,4+2/3];
% pard.scaley.Width=1/3;
% pard.scaley.TooltipString=pard.scalext.TooltipString;
% pard.scaley.Optional=true;

pard.copyalllayers_button.object=struct('Style','pushbutton','String','-> all L');
pard.copyalllayers_button.position=[9,3.4];
pard.copyalllayers_button.Width=.6;
pard.copyalllayers_button.TooltipString='Copy these parameters to all layers';
pard.copyalllayers_button.Optional=true;
pard.default_button.object=struct('Style','pushbutton','String','Default');
pard.default_button.position=[9,4];
pard.default_button.Width=.6;
pard.default_button.TooltipString='Reset to default. ';
pard.default_button.Optional=true;
pard.defaultsave_button.object=struct('Style','pushbutton','String','Save');
pard.defaultsave_button.position=[9,4.6];
pard.defaultsave_button.Width=.4;
pard.defaultsave_button.TooltipString='Save default. ';
pard.defaultsave_button.Optional=true;
%%%put in again
% pard.layercolorz.object=struct('Style','checkbox','String','layers same c/z');
% pard.layercolorz.position=[7,3.8];
% pard.layercolorz.Width=1.2; 

end

function detach_callback(a,b,obj,handle)
f=figure('MenuBar','none','Toolbar','none');
f.Units='pixel';
handle.Units='pixel';

handle.Position(1)=0;
handle.Position(2)=0;
handle.Parent=f;
f.Position(3:4)=handle.Position(3:4);
handle.Tag='detached';
% if strcmp(handle.Tag,'OV')
%     obj.ovdetached=true;
% end
end
