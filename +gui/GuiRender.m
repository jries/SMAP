classdef GuiRender< interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties
        numberOfLayers=1;
        multilayerfig;
        temproi
        layernames
    end

    methods
        function obj=GuiRender(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:}) 
            obj.propertiesToSave={'layernames'};
        end
 
        function makeGui(obj)            
            %make channel tabs
            h.layertab=uitabgroup(obj.handle,'SelectionChangedFcn',{@selectLayer_callback,obj});
            obj.adjusttabgroup(h.layertab);
            h.tab_layer1=uitab(h.layertab,'Title',['Layer' num2str(1)],'Tag','Layer1');
            
            h.tab_addlayer=uitab(h.layertab,'Title','+');   
            h.reconstruct=uicontrol(obj.handle,'Units','pixels','Position',[7, 6,150,35],'String','Render','Tag','reconstructbutton','FontSize',obj.guiPar.fontsize*1.5,...
                'Callback',@obj.render_callback);
            
            %formatGui
            hmain=obj.getPar('mainGuihandle');
            pp=hmain.Position;
            h.formatgui=uipanel('Parent',hmain,'Units','pixel','Position',[pp(3)*0.8 400 pp(3)*0.2 350]);
            h.filterpanel=uipanel('Parent',hmain,'Units','pixel','Position',[0 400 430 350],'Visible','off');
            obj.setPar('filterpanel',h.filterpanel);
            formatgui=gui.GuiFormat(h.formatgui,obj.P);
            formatgui.attachLocData(obj.locData);
            formatgui.makeGui;
            obj.setPar('guiFormat',formatgui);
            
            obj.children.guiFormat=formatgui;
            obj.setPar('numberOfLayers',1);
            h.layer_1=obj.addlayer(h.tab_layer1,1);
            obj.setPar('layernames',{'Layer1'});
            obj.guihandles=h;
            
            addlistener(obj.P,'sr_render',@obj.render_notify);
            addlistener(obj.P,'sr_display',@obj.display_notify);  
            
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
            h.layertab.UIContextMenu=c;
            
            m3 = uimenu(c,'Label','add layer','Callback',{@menu_callback,obj});
            m1 = uimenu(c,'Label','remove layer','Callback',{@menu_callback,obj});
            m4 = uimenu(c,'Label','rename layer','Callback',{@menu_callback,obj});
            if ispc
                posmen='tlo';
                shiftmen=[10 0];
            else
                posmen='tli';
                shiftmen=[10 -10];
            end
            makemenuindicator(h.layertab,posmen,shiftmen);
        end
        
        function setGuiParameters(obj,p,setchildren,setmenulist)
%             add layers if needed
            %remove all layers
            fn=fieldnames(p.children);
            currentlayers=fieldnames(obj.children);
            for k=2:length(currentlayers)-1
                obj.removelayer(k)
            end
            if isfield(p,'layernames')
                obj.guihandles.layertab.SelectedTab.Title=p.layernames{1};
            end
            indl=2;
%             layertabnames={obj.guihandles.layertab.Children(:).Title};
            for k=2:length(fn)
                if length(fn{k})>5&&strcmp(fn{k}(1:5),'Layer') && ~strcmp(fn{k}(6:end),'1')
%                     if ~any(strcmp(layertabnames,fn{k}))
                        eventdata.NewValue.Title='+';
                        selectLayer_callback(obj.guihandles.layertab,eventdata,obj)
                        if isfield(p,'layernames')
                            obj.guihandles.layertab.SelectedTab.Title=p.layernames{indl};
                            indl=indl+1;
                        end
%                     end
                end
            end
            setGuiParameters@interfaces.GuiModuleInterface(obj,p,setchildren,setmenulist);            
        end
        
        function hpanel=addlayer(obj,handle,k)
             hpanel=uipanel(handle,'Unit','pixels','Position',obj.guiPar.tabsize2);        
             layer=gui.GuiChannel(hpanel,obj.P);
             tag=handle.Tag;
             layer.layer=k;
             layer.attachLocData(obj.locData);
             p.Vsep=3;p.FieldHeight=30;
             layer.setGuiAppearence(p);
             layer.makeGui;
             if isfield(obj.children,'Layer1')
                 pold=getGuiParameters(obj.children.Layer1);
                 s=obj.getPar('filtertable','layer',1);
                 obj.setPar('filtertable',s,'layer',layer.layer);
                 layer.setGuiParameters(pold);
                 layer.setvisibility;
                 layer.setfiltergray;
             end
             layer.updateLayerField;
             obj.children.(tag)=layer;
             layer.name=tag;
%              obj.locData.layer(k).filter=[];
%              obj.locData.layer(k).groupfilter=[];
             obj.locData.filter([],k);           
        end
        
        function display_notify(obj,lp,eventdata)
            
            rf=obj.getPar('fastrender');
            if isempty(rf)
                rf=false;
            end
            
%             if ~rf
%             obj.status('display')
%             end
            
            guiformat=obj.getPar('guiFormat');
            proi=guiformat.roiset;
            if (proi.isvalid)
                obj.temproi=proi;
            end
            
            hfig=obj.getPar('sr_figurehandle');
            if ~isvalid(hfig)
                hfig=figure('Name','Reconstructed superresolution image');
                obj.setPar('sr_figurenumber',hfig.Number);
            end
%             displayer=Displayer(obj.locData);
%             pk=obj.getLayerParameters(1,displayer.inputParameters);
            pk=obj.getLayerParameters(1,displayerSMAP);
            cs=obj.getPar('sr_colorbarthickness');
            if ~isempty(cs)
                pk.sr_colorbarthickness=cs;
            end
%             displayer.setParameters(pk);
            [finalImage,sr_imagehandle]=displayerSMAP(obj.locData.layer,pk);
            

%             [finalImage,sr_imagehandle]=displayer.displayImage(obj.locData.layer);
            if ~rf
            obj.setPar('sr_imagehandle',sr_imagehandle);
            obj.setPar('sr_image',finalImage);
            end
            
            if ~rf && ~isempty(obj.temproi)&&obj.temproi.isvalid
             guiformat.roiset(obj.temproi);
            end
             
        end
        
        function draw(obj)
%             obj.status('draw')
            lp=obj.locData;
%             drawer=Drawer(obj.locData);
            layerson=obj.getPar('sr_layerson');
            for k=1:obj.numberOfLayers
                pk=obj.getLayerParameters(k,drawerSMAP);
                if layerson(k)
                    lp.layer(k).images.finalImages=drawerSMAP(lp.layer(k).images.srimage,pk); 
%                     drawer.setParameters(pk);
%                     lp.layer(k).images.finalImages=drawer.drawImage(lp.layer(k).images.srimage);            
                end
            end
            notify(obj.P,'sr_display')
%             obj.status('draw done')
        end
        function render_callback(obj,object,eventdata)
            hfig=obj.getPar('sr_figurehandle');
            if ~isvalid(hfig)
                f=obj.getPar('sr_figurenumber');
                fh=figure(f);
                fh.Name='Reconstructed superresolution image';
                obj.setPar('sr_figurenumber',f);
            end
            hfig=obj.getPar('sr_figurehandle');
             figure(hfig);
            obj.render_notify;
        end
        function render_notify(obj,object,eventdata)
            
%             drawnow
            lp=obj.locData;
            extraspace=150;
            pos=obj.getPar('sr_pos');
            sizesr=obj.getPar('sr_size');
%             xmin=pos(1)-sizesr(1)-extraspace;xmax=pos(1)+sizesr(1)+extraspace;
%             ymin=pos(2)-sizesr(2)-extraspace;ymax=pos(2)+sizesr(2)+extraspace;
%             renderer=Renderer(obj.locData);
            isfast=obj.getPar('fastrender');
            if isempty(isfast)
                isfast=false;
            end
            if ~isfast
            obj.status('render') 
            drawnow
            end
            layeron=obj.getPar('sr_layerson');
            for k=1:obj.numberOfLayers
                pk=obj.getLayerParameters(k,renderSMAP);
%                 pk=obj.getLayerParameters(k,renderer.inputParameters);
                if layeron(k)
%                     indin=[];
                    if isfast
                        pk.groupcheck=false;
%                         pk.sr_pixrec=pk.sr_pixrec*2;

                            posh=[pk.sr_pos(1) pk.sr_pos(2) pk.sr_size(1)*2 pk.sr_size(2)*2];
                            fields={'xnm','ynm'};
                            if strcmp(pk.rendermode.selection,'Gauss')
                                pk.rendermode.selection='constGauss';
                                pk.rendermode.Value=3;
                                    calculatesigma=true;
                                    fields{end+1}='locprecnm';
                            else 
                                calculatesigma=false;

                            end
%                             if strcmp(pk.rendermode.selection,'DL')
%                                 fields{end+1}='P';
%                             end
                            switch (pk.render_colormode.selection)
                                case 'z'
                                    fields{end+1}='znm';
                                case 'field'
                                    fields{end+1}=pk.renderfield.selection;
                            end
                            switch pk.intensitycoding.selection
                                case {'blinks','√blinks'}
                                    fields{end+1}='numberInGroup';
                                case {'photons','√photons'}
                                    fields{end+1}='phot';
                            end
%                             {'xnm','ynm','znm','locprecnm','PSFxnm','phot',pk.renderfield.selection}
                            locD=obj.locData.getloc(fields,'layer',k,'position',posh);
                            if calculatesigma
                                
                                pk.renderparameter=max(pk.sr_pixrec*pk.mingausspix,myquantilefast(locD.locprecnm,0.4,10000));
                            end

%                             maxlfast=2e5;
%                             if length(locD.xnm)>maxlfast
%                                 indin=false(size(locD.xnm));
%                                 gap=ceil(length(locD.xnm)/maxlfast);
%                                 indin(1:gap:end)=true;
%                                 locD=copystructReduce(locD,indin);
%                             end
                            ld=interfaces.LocalizationData;
                            ld.files=obj.locData.files;
                            
                            ld.attachPar(obj.P);
                            ld.loc=locD;
                            ld.layer(k)=ld.layer(1);
%                             ld.grouploc=locD;
                            locDat=ld;
                    else
                        obj.locData.filter('channel',k,'inlist',pk.channels) 
                        locDat=obj.locData;
                    end

                    if strcmp(pk.rendermode.selection,'Other')
                        modules=obj.getPar('rendermodules');
                        if length(modules)<k || isempty(modules{k})
                            warndlg('please select external renderer first')
                            return
                        end
                        lp.layer(k).images.srimage=modules{k}.render(obj.locData,pk);
                    else
                        lp.layer(k).images.srimage=renderSMAP(locDat,pk,k);        
%                     lp.layer(k).images.srimage=renderer.render(k);  
                    end
                end
            end
            obj.draw;
            if ~isfast
            obj.status('render done')  
            end
        end
    
        function removelayer(obj,number)
            layername=['Layer' num2str(number)];
            tabname=obj.children.(layername).name;
            lon=obj.getPar('sr_layerson');
            lon(number)=[];
            lon(end+1)=0;
            obj.setPar('sr_layerson',lon);
            
            %         numold=obj.numberOfLayers;
%         
         obj.setPar(['layer' num2str(number) '_layercheck'],0);
         if length(obj.locData.layer)>=number
            obj.locData.layer(number)=[];
         end
         deletechildren(obj.children.(layername));
         chn=fieldnames(obj.children);
         for k=number+1:length(chn)-1
             obj.children.(['Layer' num2str(k-1)])=obj.children.(['Layer' num2str(k)]);
         end
         obj.children=rmfield(obj.children,chn{end});
         
         
         tab=findobj(obj.guihandles.layertab,'Title',tabname);
        delete(tab);
         obj.numberOfLayers=obj.numberOfLayers-1;
         obj.setPar('numberOfLayers',obj.numberOfLayers);
         tab=obj.guihandles.layertab.Children(1);
         obj.guihandles.layertab.SelectedTab=tab;
        updatelayernames(obj)
        end
    end
end

function selectLayer_callback(tabgroup,eventdata,obj)
if isprop(eventdata.NewValue,'Tag')
    layer=(eventdata.NewValue.Tag);
else
    layer=[];
end
layertitle=(eventdata.NewValue.Title);
if strcmp(layertitle,'+')
    newlayernumber=obj.numberOfLayers+1;
    tag=['Layer' num2str(newlayernumber)];
    title=tag;
    while any(strcmp({obj.guihandles.layertab.Children.Title},title))
        title=[title 'n'];
    end
    
      h.(['tab_layer' num2str(newlayernumber)])=uitab(obj.guihandles.layertab,'Title',title,'Tag',tag);
       h.(['layer_' num2str(newlayernumber)])=obj.addlayer(  h.(['tab_layer' num2str(newlayernumber)]),newlayernumber);
    obj.numberOfLayers=newlayernumber;
    obj.setPar('numberOfLayers',newlayernumber);
    s=1:length(tabgroup.Children);
    s(end-1)=s(end);
    s(end)=s(end)-1;
    tabgroup.Children=tabgroup.Children(s);
    tabgroup.SelectedTab=tabgroup.Children(end-1); 
    number=newlayernumber;
    selectedlayer=['Layer' num2str(number)];
    obj.setPar('layercheck',true,'Value','layer',selectedlayer)
else
    number=find(strcmp({obj.guihandles.layertab.Children.Title},eventdata.NewValue.Title));
    selectedlayer=['Layer' num2str(number)];
%     selectedlayer=sscanf(layer,'Layer%d');
end
for k=1:length(tabgroup.Children)-1
%     layernumber=tabgroup.Children(k).Title;
%     layertag=tabgroup.Children(k).Tag;
    layertag=['Layer' num2str(k)];
    obj.children.(layertag).setlayer(number);
end
updatelayernames(obj)
end

function menu_callback(callobj,b,obj)
switch callobj.Label
    case 'add layer'
        tab=findobj(obj.guihandles.layertab,'Title','+');
        obj.guihandles.layertab.SelectedTab=tab;
        eventdata.NewValue=tab;
        selectLayer_callback(obj.guihandles.layertab,eventdata,obj)
    case 'remove layer'
        %delete
        if obj.numberOfLayers<2
            disp('one layer required, cannot delete layer 1');
             return
        end
        number=find(strcmp({obj.guihandles.layertab.Children.Title},obj.guihandles.layertab.SelectedTab.Title));
%         tag=obj.guihandles.layertab.SelectedTab.Tag;
        obj.removelayer(number);
        return
%         numold=obj.numberOfLayers;
%         
%         name=['Layer',num2str(numold)];
%         p=obj.getPar('sr_layerson');
%         p(numold)=0;
%         obj.setPar('sr_layerson',p);
%         obj.setPar(['layer' num2str(numold) '_layercheck'],0);
%         obj.locData.layer(numold)=[];
%         deletechildren(obj.children.(name));
%         tab=findobj(obj.guihandles.layertab,'Title',['Layer',num2str(numold)]);
%         delete(tab);
%         obj.numberOfLayers=obj.numberOfLayers-1;
%         obj.setPar('numberOfLayers',obj.numberOfLayers);
% %         tab=findobj(obj.guihandles.layertab,'Title','Layer1');
%         tab=obj.guihandles.layertab.Children(1);
%         obj.guihandles.layertab.SelectedTab=tab;
    case 'rename layer'
        currentname=obj.guihandles.layertab.SelectedTab.Title;
        tag=obj.guihandles.layertab.SelectedTab.Tag;
        newname=inputdlg('New layer name: ','rename layer',1,{currentname});
        if isempty(newname)
            return
        end
        newname=newname{1};
        obj.children.(tag).name=newname;
        obj.guihandles.layertab.SelectedTab.Title=newname;
        updatelayernames(obj)
        
        
end
end


function updatelayernames(obj)
k=1;
names={};
while isfield(obj.children,['Layer' num2str(k)])
    if isvalid(obj.children.(['Layer' num2str(k)]))
        names{k}=obj.children.(['Layer' num2str(k)]).name;
    end
            k=k+1;
end
obj.setPar('layernames',names);
obj.layernames=names;
end        