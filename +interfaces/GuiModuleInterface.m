classdef GuiModuleInterface<interfaces.GuiParameterInterface
    %interface to provide GUI functionality to modules. 
    %The GUI is described with a simple syntax, its appearence is freely
    %configurable
    % extends GuiParameterInterface to deal with GUI specific parameters
    % which are connected to uicontrols
    %distinction between GuiParameterInterface and GuiModuleInterface not
    %very clear, could have been one class.
    properties
        handle %figure or panel where the GUI is created
        guiPar %structure which defines global parameters for gui appearence (e.g. font size)
        
        %modules created inside a module should be put here. 
        %Many functions recursively search children to e.g. retrieve and
        %restore all gui settings or change the appearence
        children
        plugininfo %space to put info about plugin. 
        excludeFromSave={};% uicontrols with this name are not saved or read out by getGuiParameters
        propertiesToSave={}; %{fileds} obj.(field) is saved/restored with getGuiParameters.
        pluginpath %{path1,path2,filename}: where plugin was stored before creation
        guihandles %handles to gui objects are stored here. ONLY flat structure. All handles should be added, otherwise they will not be resized etc.
        simplegui=false;
        guiselector=struct('position',[],'show',false);
       
    end
    methods
         function obj=GuiModuleInterface(varargin)
             %varargin= obj.handle, obj.P (global parameter object)
            if nargin>0
                obj.handle=varargin{1};
            end
            if nargin>1
                obj.attachPar(varargin{2});
            else
                obj.attachPar(interfaces.ParameterData);
            end
            
            %PC-Mac differences
            if ispc
                 guiPar.fontsize=10;
                 guiPar.FieldHeight=26;
                 guiPar.tabsize1=[0    -1  546 342];
                 guiPar.tabsize2=[0    -1  540 311];    
                 guiPar.Vsep=3;
                 guiPar.Xrim=3;
                 guiPar.Vrim=2;  
            else
                guiPar.fontsize=15;
                guiPar.FieldHeight=25;
                guiPar.tabsize1=[1    0  527 324];
                guiPar.tabsize2=[0    -1  521 292];
                guiPar.Vsep=1;
                guiPar.Xrim=2;
                guiPar.Vrim=0; 
            end
            
            %initialze guiPar
              
            guiPar.Xsep=1;
            guiPar.Vsep=1;
            
            guiPar.Xpos=1;
            guiPar.Vpos=1;

            if ishandle(obj.handle)
                hpos=obj.handle.Position;
                guiPar.FieldWidth=(hpos(3)-3*guiPar.Xsep-2*guiPar.Xrim)/4;
            end
            obj.guiPar=guiPar;
         end
        
        function pout=getGuiParameters(obj,getchildren,onlyedit)
            %gets all GUI parameters (editable parts of uicontrols etc) in
            %a parsed format. 
            %p=getGuiParameters(getchildren)
            %if getchildren=true (optional): get paraemters also from children
            %if onlyedit: only editable fields are returned.
            if nargin<2
                getchildren=false;
            end
            if nargin<3
                onlyedit=false;
            end
            h=obj.guihandles;
            p=[];
            if ~isempty(h)
                fn=fieldnames(h);
                for k=1:length(fn)
                    vh=obj.getSingleGuiParameter(fn{k},onlyedit);  
%                     if ~isempty(vh)
                            p.(fn{k})=vh;             
%                     end
                end
            end
            psave=obj.propertiesToSave;
            for k=1:length(psave)
                p.(psave{k})=obj.(psave{k});
            end
            p.classname=class(obj);
            if isprop(obj,'pluginpath')
            p.pluginpath=obj.pluginpath;
            end
            pout=p;

            %children
            if getchildren
                if isstruct(obj.children)
                    guichildren=fieldnames(obj.children);
                    for k=1:length(guichildren)
                        if isvalid(obj.children.(guichildren{k}))
                            ph=obj.children.(guichildren{k}).getGuiParameters(true,onlyedit);
                        if ~isempty(ph)
                            pout.children.(guichildren{k})=ph;
                        end
                        end
                    end
                end
            end        
        end
        
        function par=getSingleGuiParameter(obj,field,onlyedit)
            %get parsed value of uicontrol
            % p=getSingleGuiParameter(guifield). 
            % uicontrol needs to be stored as obj.guihandles.(field)
            if nargin <3
                onlyedit=false;
            end
            try
            hfn=obj.guihandles.(field);
%             if onlyedit
%                 switch hfn.Style
%                     case {'edit','popupmenu','listbox','checkbox','togglebutton'}
%                         par=obj.handle2value(hfn);
%                     otherwise
%                         par=[];
%                 end
%             else
            par=obj.handle2value(hfn,onlyedit);
            catch
                par=[];
            end
%             end
        end
        function switchvisibleall(obj)
            h=obj.guihandles;
            fn=fieldnames(h);
            
            for k=1:length(fn)
                
                if isfield(h,fn{k})&&isprop(h.(fn{k}),'Callback') && length(h.(fn{k}).Callback)>1 &&contains(func2str(h.(fn{k}).Callback{1}),'switchvisible')
                    fnc=h.(fn{k}).Callback;
                    fnc{1}(h.(fn{k}),0,fnc{2});
                end
            end
        end
        function switchvisible(obj,control,data,p,callbackfnc)
            val=control.Value;
            for k=1:length(p)
                if any(p(k).value==val)
                    off=p(k).off;
                    for l=1:length(off)
                         if isfield(obj.guihandles,off{l})
                        hh=obj.guihandles.(off{l});
                        hh.Visible='off';
                        if isprop(hh,'Callback') && ~isempty(hh.Callback) && iscell(hh.Callback) && contains(func2str(hh.Callback{1}),'switchvisible')
                            fnc=hh.Callback;
                            fnc{1}(hh,0,fnc{2});
                        end
                         end
%                         obj.guihandles.(off{l}).Visible='off';
                    end
                    on=p(k).on;
                    for l=1:length(on)
                        if isfield(obj.guihandles,on{l})
                            hh=obj.guihandles.(on{l});
                            if ~(myisfield(obj.guidef.(on{l}),'Optional') && obj.guidef.(on{l}).Optional &&obj.simplegui)
                            hh.Visible='on';
                            end
                            if isprop(hh,'Callback') && ~isempty(hh.Callback) &&  iscell(hh.Callback) && contains(func2str(hh.Callback{1}),'switchvisible')

                            %call switchvisible for all children to nest
                            %functions
                                fnc=hh.Callback;
                                fnc{1}(hh,0,fnc{2});
                            end
                        end
                        
                    end                    
                end
            end
            if nargin>4 && ~isempty(callbackfnc)
                callbackfnc{1}(callbackfnc{2:end});
            end
            
        end
        
        function fieldvisibility(obj,varargin)
            %sets the gui state (simple/advanced) and sets visibility of
            %certain fields
            %Arguments:
            %'guistate' 'simple'/'advanced' 1/0
            %'on', {fields}
            %'off', {fields}
            if nargin<2
                p.on={};p.off={};p.guistate=obj.simplegui;
            else
            p=fieldvisibiltyparser(varargin);
            end
            pard=obj.guidef;
            if ~isstruct(pard)
                return
            end
            fn=fieldnames(pard);
            oldstate=obj.simplegui;
            switch p.guistate
                case {'S','s','simple',1,true}
                    obj.simplegui=true;
                    obj.guihandles.simplegui.String='v';
                    obj.guihandles.simplegui.Value=1;
                    
                    for k=1:length(fn)
                        if isfield(pard.(fn{k}),'Optional')&&pard.(fn{k}).Optional==true
                            if oldstate==false %if was advanced
                                obj.guihandles.(fn{k}).UserData.visibleSMAP=obj.guihandles.(fn{k}).Visible;
                            end
                            obj.guihandles.(fn{k}).Visible='off';
                        end
                    end
%                     if previeously advanced: write visible status into
%                     guihandles.field.smapVisible for optional parameters,
%                     switch off those parameters
                case {'A','a','advanced',0,false}
%                     switch on all optional fields with
%                     smapVisible=visible
                    obj.simplegui=false;
                    obj.guihandles.simplegui.String='-';
                    obj.guihandles.simplegui.Value=0;
                    for k=1:length(fn)
                        if isfield(pard.(fn{k}),'Optional')&&pard.(fn{k}).Optional==true
                            if isfield(obj.guihandles.(fn{k}).UserData,'visibleSMAP')       
                                 obj.guihandles.(fn{k}).Visible=obj.guihandles.(fn{k}).UserData.visibleSMAP;
                            elseif oldstate==true %only change if from simple to advanced
                                obj.guihandles.(fn{k}).Visible='on';
                            end   
                        end
                    end
            end

            if ~iscell(p.on)
                p.on={p.on};
            end
            for k=1:length(p.on)
                if isfield(obj.guihandles,p.on{k})
                    if ~obj.simplegui||~(isfield(pard.(p.on{k}),'Optional')&&pard.(p.on{k}).Optional==true) %not an optional parameter
                        obj.guihandles.(p.on{k}).Visible='on';
                    end
                    obj.guihandles.(p.on{k}).UserData.visibleSMAP='on';
                end
            end
            if ~iscell(p.off)
                p.off={p.off};
            end            
            for k=1:length(p.off)
                if isfield(obj.guihandles,p.off{k})      
                    obj.guihandles.(p.off{k}).Visible='off';
                    obj.guihandles.(p.off{k}).UserData.visibleSMAP='off';
                end
            end            
%             for all p.on: smapvisible='on'. if optional: switch on if state=advanced. 
%             for all p.off: smapvisible, visible='off';
        end
        function setGuiParameters(obj,p,setchildren,setmenulist)
            %sets parameters in GUI uicontrols
            %setGuiParameters(p,setchildren)
            % if setchildren=true (optional): set also gui parameters in
            % children
            if isempty(p)
                return
            end
            if nargin<3
                setchildren=false;
            end
            if nargin<4
                setmenulist=true;
            end
            if isstruct(p)
                fn=fieldnames(p);
                phere=p;
                h=obj.guihandles;
                for k=1:length(fn)
                    if isfield(h,fn{k})&&isprop(h.(fn{k}),'Style')&&~strcmp(h.(fn{k}).Style,'text')&&~any(ismember(obj.excludeFromSave,fn))                        
                        
                        hs=obj.value2handle(phere.(fn{k}),h.(fn{k}));                      
%                         if (strcmp(h.(fn{k}).Style,'popupmenu'))
%                             htmp.Value=hs.Value;
%                             if (iscell(hs.String)&&hs.Value>length(hs.String)||(~iscell(hs.String)&&hs.Value>size(hs.String,1)))
%                                 htmp.Value=1;
%                             end
%                             hs=htmp;                   
%                         end
                        if ~setmenulist && strcmp(h.(fn{k}).Style,'popupmenu')&&isprop(h.(fn{k}),'String')
                            hs=myrmfield(hs,'String');
                            hs.Value=min(hs.Value, length(h.(fn{k}).String));
                        end
                        h.(fn{k})=copyfields(h.(fn{k}),hs);
%                     elseif strcmp(fn{k},'globaltable')
                    elseif isfield(h,fn{k}) && isa(h.(fn{k}),'matlab.ui.control.Table') %Table    
                        hs=obj.value2handle(phere.(fn{k}),h.(fn{k})); 
                        h.(fn{k})=copyfields(h.(fn{k}),hs);
                    end
%                     if strcmp(fn{k},'isscmos')
                        if isfield(h,fn{k})&&isprop(h.(fn{k}),'Callback') && length(h.(fn{k}).Callback)>1 &&contains(func2str(h.(fn{k}).Callback{1}),'switchvisible')
                            fnc=h.(fn{k}).Callback;
                            fnc{1}(h.(fn{k}),0,fnc{2});
                        end
%                     end
                end
                
                psave=obj.propertiesToSave;
                for k=1:length(psave)
                    if isfield(p,psave{k})
                        obj.(psave{k})=p.(psave{k});
                    end
                end
                %put output parameters to P
                if isstruct(obj.guihandles)
                fo=intersect(intersect(fn,obj.outputParameters),fieldnames(obj.guihandles));
                for k=1:length(fo)
                    obj.updateGuiParameter(0,0,fo{k})
                end
                end
                
            elseif iscell(p) %handle, value
                hs=obj.value2handle(p{2},p{1});
                
                copyfields(p{1},hs);
            end
            
            if setchildren&&isfield(p,'children')
                fn=fieldnames(p.children);
                for k=1:length(fn)
                    if isfield(obj.children,fn{k})
                        child=obj.children.(fn{k});
                        pchild=p.children.(fn{k});
                        try
                            child.setGuiParameters(pchild,true,setmenulist);
                        catch err
                            child
                            err
                        end
                    
                    end
                end
            end
%             obj.executecallback('switchvisible');
%             obj.initializeGuiParameters;
        end

        function p=getAllParameters(obj,inputParameters,editonly)
            % gets Gui Parameters (without children) and inputParameters
            % p=getAllParameters(inputParameters)
            % inputParameters can be omitted: then obj.inputParameters are
            % used
            if nargin<2||isempty(inputParameters)
                inputParameters=obj.inputParameters;
            end
            if nargin<3
                editonly=false;
            end
            if any(strcmp(inputParameters,'layers'))
                for k=1:obj.getPar('numberOfLayers')
                    inputParameters{end+1}=['layer' num2str(k) '_'];
                end
            end
            p=getAllParameters@interfaces.ParameterInterface(obj,inputParameters);
            p=copyfields(p,obj.getGuiParameters(false,editonly));
            
        end
        
        function p=getLayerParameters(obj,layeri,inputParameters)
            %get parameters of a specific layer as specified in input
            %Parameters
            %p=getLayerParameters(layer,inputParameters)
            % if inputParameters empty: use obj.inputParameters
            %layer is layer number or vector of layers. If layer is empty:
            %use all layers.
            %p is cell array of structures
            if nargin<3
                inputParameters=fieldnames(obj.P.par);
            end
            if nargin<2
                layeri=[];
            end
            if isempty(layeri)
                layer=1:obj.getPar('numberOfLayers');
            else
                layer=layeri;
            end
            pall=obj.getAllParameters(inputParameters,false);
            for k=1:length(layer)
                p{k}=pall;
                lp=obj.getPar('','layer',layer(k));
                p{k}=copyfields(p{k}, lp);
%                 p{k}=copyfields(p{k}, lp.rec_addpar);
            end
            if length(layeri)==1
                p=p{1};
            end
        end

        function resize(obj,factor)
            % resizes all GUI and children GUIs
            % usually called from figure.SizeChangeCallback (or similar)
            
            if myisfield(obj,'guihandles') && ~isempty(obj.guihandles)
            fn=fieldnames(obj.guihandles);
            for k=1:length(fn)
                try
                obj.guihandles.(fn{k}).FontSize=obj.guihandles.(fn{k}).FontSize*factor;
                catch err
                end
                try
                    if strcmpi(obj.guihandles.(fn{k}).Units,'pixels')
                        obj.guihandles.(fn{k}).Position=obj.guihandles.(fn{k}).Position*factor;
                    end
                catch err
                end
                try
                    if isa(obj.guihandles.(fn{k}),'matlab.ui.control.Table')
                    obj.guihandles.(fn{k}).ColumnWidth=num2cell([obj.guihandles.(fn{k}).ColumnWidth{:}]*factor);
                    end
                catch err
                end
            end
            end
            if ~isempty(obj.children)
                ch=fieldnames(obj.children);
                for k=1:length(ch)
                    obj.children.(ch{k}).resize(factor);
                end
            end    
        end
        
        function  adjusttabgroup(obj,htg)
            %adjusts width of second tabgroup on mac
            if ispc
            else
                htg.Units='pixel';
                htg.Position(1)=htg.Position(1)-8;
                htg.Position(3)=htg.Position(3)+16;
                htg.Position(2)=htg.Position(2)-12;
                htg.Position(4)=htg.Position(4)+16;
                htg.Units='normalized';
            end
        end
        
        function setGuiAppearence(obj,p)
            % sets guiPar: (e.g. font size, field height etc)
            %setGuiAppearence(p): p structure with any of
            % Vrim=0; vertical rim (space above and below)
            % Xrim=0; horizontal rim (space right and left)
            % Xsep=1; horizontal space between controls
            % Vsep=1; vertical space between controls
            
            % Xpos=1; horizontal position of GUI (in GUI-units)
            % Vpos=1; vertical GUI position
            % fontsize 
            % FieldWidth: x-extension of 1 GUI unit (minus Xrim) in pixels
            %usually not set but calculated from size of
            % figure/panel in obj.handle
            % FieldHeight: y-extenson of 1 GUI unit: height of uicontrol, 
            
            fn=fieldnames(p);
            for k=1:length(fn)
                obj.guiPar.(fn{k})=p.(fn{k});
            end        
        end 

       function anyoptional=makeGui(obj,guidef)
           % renders the GUI according to guidef, then calls obj.initGui.
           % if guidef not passed on: calls obj.guidef (that is the usual
           % way of defining a GUI)
            if nargin==1
                guidef=obj.guidef;
            end
            
            if ~isempty(obj.handle)
                obj.handle.Units='pixels';
                hpos=obj.handle.Position;
                obj.guiPar.FieldWidth=(hpos(3)-3*obj.guiPar.Xsep-2*obj.guiPar.Xrim)/4;
            end
            if isstruct(guidef)
                
                guiPar=obj.guiPar;
                if isfield(guidef,'locselector') && any(guidef.locselector)
                    guidef=copyfields(guidef,locselectordefault);
                    if isfield(guidef,'syncParameters')
                        guidef.syncParameters{end+1}={'filelist_short_ext','selector_filelist','String'};
                    else
                        guidef.syncParameters={{'filelist_short_ext','selector_filelist',{'String'}}};
                    end
                    fl=obj.getPar('filelist_short_ext');
                    if ~isempty(fl)
                        if isfield(fl,'String')
                            guidef.selector_filelist.object.String=fl.String;
                        else
                            guidef.selector_filelist.object.String=fl;
                        end
                    end
                    offstr={'off','on'};
                    if length(guidef.locselector)==3
                        guidef.selector_filelist.object.Visible=offstr{guidef.locselector(1)+1};
                        guidef.selector_filter.object.Visible=offstr{guidef.locselector(2)+1};
                        guidef.selector_pos.object.Visible=offstr{guidef.locselector(3)+1};
                    end
%                     guidef.syncParameters %extend 
                end
                
                %tabs handling
                if isfield(guidef,'tab')
                    tabgroup=uitabgroup(obj.handle);
                    tabgroup.Position(4)=tabgroup.Position(4)*(1-guiPar.Vrim/obj.handle.Position(4));
                    guiPar.Vrim=guiPar.Vrim+20;
                    fn=fieldnames(guidef.tab);
                    for k=1:length(fn)
                        obj.guihandles.(fn{k})=uitab(tabgroup,'Title',guidef.tab.(fn{k}));
                    end
                    istab=true;
                    maintab=obj.guihandles.(fn{1});
                    guidef=rmfield(guidef,'tab');
                else
                    istab=false;
                end
                
                allFields=fieldnames(guidef);
                anyoptional=false;
                synchronizeguistate=obj.getPar('synchronizeguistate');
                
                %help file
                pp=obj.pluginpath;
                if iscell(pp) && length(pp)==3
                    helpfile=[pp{1} '.' pp{2} '.' pp{3} '.txt'];
                else
                    helpfile='';
                end
                
                for k=1:length(allFields) 
                    thisField=guidef.(allFields{k});
                    if strcmp(allFields{k},'syncParameters')
                        obj.syncParameters=thisField;
                    elseif strcmp(allFields{k},'inputParameters')
                        obj.inputParameters=thisField;
                    elseif strcmp(allFields{k},'outputParameters')
                        obj.outputParameters=thisField;
                    elseif strcmp(allFields{k},'plugininfo')
                        obj.plugininfo=thisField;
                    elseif strcmp(allFields{k},'locselector')
                        %do nothing here
                    elseif strcmp(allFields{k},'helpfile')
                        helpfile=thisField;
                    elseif isstruct(thisField) && ~isempty(obj.handle) %results name
                        if ~isfield(thisField,'object') || ~isfield(thisField.object,'Style')
                            allFields{k}
                            thisField
                            str=['guidef definition is incomplete in classe:' class(obj)];
                            warning(str)
                            continue
                            
                        end
                        if isfield(thisField,'Optional')&&thisField.Optional
                            anyoptional=true;
                        end
                        h=thisField.object;
%                         h=uicontrol(obj.handle,thisField.object);
                        h.FontSize=guiPar.fontsize;
                        h.Units='pixels';
%                         set(h,'Units','pixels')

                        if isfield(thisField,'Width')
                            widthf=thisField.Width;
                        else
                            widthf=1;
                        end

                        if isfield(thisField,'Height')
                            heightf=thisField.Height;
                        else
                            heightf=1;
                        end

                        switch h.Style
                            case {'pushbutton','togglebutton'}
                                hadjust=4;
                            case 'popupmenu'
                                hadjust=0;
                            otherwise
                                hadjust=0;
                        end

                        h.Position=[(guiPar.FieldWidth)*(thisField.position(2)-1+guiPar.Xpos-1)+guiPar.Xrim, ...
                            hpos(4)-(thisField.position(1)+guiPar.Vpos-1)*guiPar.FieldHeight-guiPar.Vrim-hadjust/2-guiPar.Vsep-3, ...
                           guiPar.FieldWidth*widthf-2*guiPar.Xsep,...
                            guiPar.FieldHeight*heightf-guiPar.Vsep+hadjust];
                        
                        
                        if strcmpi(h.Style,'text')
                            h.HorizontalAlignment='left';
                        end
                        
                        if istab
                            if isfield(thisField,'tab')
                                parenth=obj.guihandles.(thisField.tab);
                            else
                                parenth=maintab;
                            end
                        else
                            parenth=obj.handle;
                        end
                        hg=uicontrol(parenth,h);
                        %bug: sometimes string does not get passed on
                        if isfield(h,'String')
                            hg.String=h.String;
                        end
                        
                        obj.guihandles.(allFields{k})=hg;
                        thisField=myrmfield(thisField,{'Width','Height','position','object','load'});
                        remaining=fieldnames(thisField);
%                         remaining=setdiff(fieldnames(thisField),{'Width','Height','position','object','load'});
                        for kr=1:length(remaining)
                            if isprop(hg,remaining{kr})
                                hg.(remaining{kr})=thisField.(remaining{kr});
                            end
                        end
                        
                        if isfield(thisField,'uimenu')
                            if iscell(thisField.uimenu)
                            makemenuindicator(hg,thisField.uimenu{1},thisField.uimenu{2})
                            else
                                makemenuindicator(hg,thisField.uimenu);
                            end
                        end
%                         if isfield(thisField,'TooltipString')
%                             h.TooltipString=thisField.TooltipString;
%                         end
                    end
                    if anyoptional && synchronizeguistate
                        obj.addSynchronization('globalGuiState',[],[],@obj.setglobalguistate);
                    end
                       
                end 
                % read help file and write tool tips 
                settingsdir=obj.getPar('SettingsDirectory');
                helpfilep=[fileparts(settingsdir) filesep 'Documentation' filesep 'help' filesep helpfile];
                if strcmp(helpfilep(1),filesep)
                    helpfilep(1)=[];
                end
                maxwidth=60;             
                if ~isempty(helpfile) && exist(helpfilep,'file')
                    [description,tooltips,interpreter]=parsehelpfile(helpfilep);
                    obj.plugininfo.description=(description);
                    obj.plugininfo.descriptioninterpreter=interpreter;
                    if ~isempty(tooltips)
                        fnt=fieldnames(tooltips);
                        for tt=1:length(fnt)
                            if isfield(obj.guihandles,fnt{tt})
                                strtt=tooltips.(fnt{tt});
                                strWrapped=mytextwrap(strtt,maxwidth,10);
                                obj.guihandles.(fnt{tt}).Tooltip=sprintf(strWrapped);
                            end
                        end
                    end
                else
%                     writehelpfile(helpfile,guidef)
                end          
            end
              
            obj.initGui; %exchanged 
            obj.setSyncParameters;
            obj.initializeGuiParameters;
       end
        
       function initGui(obj) 
           % implement if needed in subclass. Called after makeGui
           %dummy, in case not implemented
       end
       
       function status(obj,txt)
           %set the status in the status line
           numchar=65;
           if length(txt)>numchar
               txt2{1}=txt(1:numchar);
               txt2{2}=txt(numchar+1:end);
           else
               txt2=txt;
           end
           obj.setPar('status',txt2);
       end
       function setglobalguistate(obj,a,b)
           simplestate=obj.getPar('globalGuiState');
               obj.fieldvisibility('guistate',simplestate)
       end
       function p=guidef(obj)
           %overwrite in module when defining a GUI.
           %p.fieldname.object=struct('Style','edit','String','x',...)
           %    defines uicontrol with parameters from struct
           %p.fieldname.position=[row,column] defines position of uicontrol
           %    in GUI coordinates (usually four units horizontally, and one
           %    unit vertically given by obj.guiPar.FieldHeight. These need not
           %    be integer
           %p.fieldname.Width, p.fieldname.Height: Widht and Height in GUI
           %    coordinates
           %p.fieldname.TooltipString : add tooltip
           %p.syncParameters={'parameterName',field,syncmode}, field is
           %    position of uicontrol as in obj.guihandles.(field). Adds
           %    synchronization of uicontrol with global parameters via
           %    interfaces.GuiParameterInterface.addSynchronization
           %p.inputParamters: can be defined here as cell array of char
           %p.outputParameters: same
           %p.plugininfo.name, .description: longer text which describes the module
           
           p=[];
       end
       function info=info(obj)
           %returns an info structure. With defaults if empty.
            info=obj.plugininfo;
            if isempty(info)||~isfield(info,'name')
               try
                name=obj.pluginpath{end};
               catch
                   name=obj.subpluginpath{end};
               end
                [~,file]=fileparts(name);
                if ~isempty(file)
                    name=file;
                end
                info.name=name;
            end
            if ~isfield(info,'description')
                try
                [infopath,infofile]=fileparts(plugincell2path(obj.pluginpath));
                descfile=[infopath filesep infofile '.info'];
%                 if exist(descfile,'file')
                    fid=fopen(descfile);
                    idx=1;
                    line=fgetl(fid);
                    txt{idx}=line;
                    while ischar(line)
                        line=fgetl(fid);
                        idx=idx+1;
                        txt{idx}=line;
                    end
                    fclose(fid);
                    info.description=txt;
%                 else
%                     info.description=info.name;
%                 end
                catch
                    info.description=info.name;
                end
            end
            obj.plugininfo=info;
       end
       
       function setnormalizedpositionunits(obj)
          fn=fieldnames(obj.guihandles);
          for k=1:length(fn)
              if isprop(obj.guihandles.(fn{k}),'Units')
                  obj.guihandles.(fn{k}).Units='normalized';
              end
          end
          if isempty(obj.children)
              return
          end
          fn=fieldnames(obj.children);
          for k=1:length(fn)
              obj.children.(fn{k}).setnormalizedpositionunits;
          end
       end
       
       function makeinfobutton(obj,position)
           if nargin<2
               position='sw';
           end
           h=obj.handle;
           units=h.Units;
           h.Units='pixels';
           if isnumeric(position)
               pos(1:2)=position(1:2);
               pos(3:4)=[20 20];
           else
               switch position
                   case 'sw' %right lower
                       pos=[h.Position(3)-17,1,15,20];
                   case 'nw' %right upper
                       pos=[h.Position(3)-17,h.Position(4)-21,15,20];
                   case 'guiselector'
                       pos=[h.Position(3)-32,h.Position(4)-21,15,20];
               end
           end
           obj.guihandles.infobutton=uicontrol(h,'Style','pushbutton','String','i',...
               'Position',pos,'Callback',@obj.showinfo_callback);
           h.Units=units;
           obj.guihandles.infobutton.Tooltip='Show information for this plugin';
       end
       function showinfo_callback(obj,a,b)
           obj.showinfo;
       end
       function showinfo(obj, hp)
%         warnid='MATLAB:strrep:InvalidInputType';
%         warnstruct=warning('off',warnid);
%         obj.guihandles.showresults.Value=1;
%         showresults_callback(obj.guihandles.showresults,0,obj)
%         ax=obj.initaxis('Info');
%         hp=uifigure;
          if nargin<2
            hp=figure('MenuBar','none','Toolbar','figure');
            hp.Color='w';
            pos=[0.03,0.03,.9,.95];
            fs=obj.guiPar.fontsize;
          else
              pos=[0 0 .7 1];
              fs=12;

          end
%         hp=ax.Parent;
        
%         delete(ax);

             %  htxt=uicontrol(hp,'Style','text','Units','normalized','Position',[0,0,.9,1],...
        %      'FontSize',obj.guiPar.fontsize,'HorizontalAlignment','left','Max',100); 
        td=obj.info.description;
         if ~iscell(td)
          txt=strrep(td,char(9),' ');
          txt=strrep(txt,'\n',newline);
         else
             txt=td;
         end
         
         if isfield(obj.plugininfo,'descriptioninterpreter')
             interpreter=obj.plugininfo.descriptioninterpreter; 
         else
             interpreter='none';
         end
         
          htxt=annotation(hp,'textbox',pos,...
             'FontSize',fs,'HorizontalAlignment','left',...
             'BackgroundColor','w','FitBoxToText','off','EdgeColor','w',...
             'String',txt,'Interpreter',interpreter);
          htxt.Position=pos;
%          htxt.String=txt;
        %   htxt.Position=[0 0 1 1];
%           warning(warnstruct);
          h=uicontrol(hp,'Style','pushbutton','Units','normalized','Position',[0.9,0.0,.1,.05],'String','Edit','Callback',{@edit_callback,obj});
%           h=uibutton( hp,'push','Text','Edit','Position',[hp.Position(3)-40,hp.Position(4)-30,30,20],'ButtonPushedFcn',{@edit_callback,obj});
        end

    end
    
    methods (Access=private)
        function setSyncParameters(obj)
               for k=1:length(obj.syncParameters)
                   sp=obj.syncParameters{k};
                   if ~isempty(sp{2})&&ischar(sp{2})
                   h=obj.guihandles.(sp{2});
                   else 
                       h=[];
                   end
                   obj.addSynchronization(sp{1},h,sp{3:end});
                   
               end

               for k=1:length(obj.outputParameters)
                   po=obj.outputParameters{k};
%                    if ~isfield(obj.P.par,po) %not already attached
                       if isfield(obj.guihandles,po)&& isempty(obj.guihandles.(po).Callback) %is part of gui and has no callback
                           obj.guihandles.(po).Callback={@obj.updateGuiParameter,po};
                           obj.updateGuiParameter(0,0,po);
                       end     
%                    end
               end

        end
    
       function initializeGuiParameters(obj)
           po=obj.outputParameters;
           for k=1:length(po)
               if isfield(obj.guihandles,po{k})
               obj.updateGuiParameter(0,0,po{k})
               end
           end 
           pi=obj.inputParameters;
           for k=1:length(pi)
               if ~isempty(obj.guihandles)&&isfield(obj.guihandles,pi{k})
                   v=obj.getPar(pi{k});
                   obj.setGuiParameters(struct(pi{k},v))
               end
           end 
           
           pg=obj.getGuiParameters;
           ps=obj.syncParameters;
           for k=1:length(ps)
               hn=ps{k}{2};
               
               v=obj.getPar(ps{k}{1});
               if isfield(pg,ps{k}{2})
                   ph=pg.(ps{k}{2}); 
                   
                   if isstruct(ph)
                       syncv=ps{k}{3};
                       for s=1:length(syncv)
                            ph=copyfields(ph,v,syncv{s});
                       end

                   end
               else
                    ph=v;
               end
                 
                 
               if ishandle(hn)
                   
                   obj.setGuiParameters({hn,ph})                   
               elseif isfield(obj.guihandles,hn)
%                    v=obj.getPar(ps{k}{1});
                   obj.setGuiParameters(struct(ps{k}{2},ph))
               end
           end 
       end
       
      function updateGuiParameter(obj,a,b,field)
           %writes a single parameter from the GUI as specified by field
           %updateGuiParameter(a,b,field).
           %into the global parameter structure. Usually in a callback of
           %obj.guihandles.(field). Only needed when there is an explicit
           %callback. If no callback is defined but field is among
           %obj.outputParameters, this callback is created in
           %setsyncParameters
           %a,b are space holders, because updateGuiParameters is used as
           %callback.
           p=obj.getSingleGuiParameter(field);
           if isa(obj,'interfaces.LayerInterface')
              obj.setPar(field,p,'layer',obj.layer)
           else
           obj.setPar(field,p)
           end
      end
      

%        function executecallback(obj,callback)
%            if ~isstruct(obj.guihandles)
%                return
%            end
%         fn=fieldnames(obj.guihandles);
%         for k=1:length(fn)
%             h=obj.guihandles.(fn{k});
%             if isfield(h,'Callback') && ~isempty(h.Callback) && iscell(h.Callback)
%                 cfnc=h.Callback{1};
%                 if contains(func2str(cfnc),callback)
%                     cfnc(h,0,h.Callback{2});
%                 end
%             end
%         end

%        end
   end
end

function edit_callback(a,b,obj)
basedir=fileparts(obj.getPar('SettingsDirectory'));
outdir=[basedir filesep 'Documentation' filesep 'help' filesep];
if isempty(basedir)
    outdir(1)=[];
end
g=obj.guidef;
if isfield(g,'helpfile')
    helpfile=g.helpfile;
else
    pp=obj.pluginpath;
    if iscell(pp) && length(pp)==3
        helpfile=[pp{1} '.' pp{2} '.' pp{3} '.txt'];
    else
        disp('this is not a plugin')
        return
    end
end

if ~exist([outdir helpfile],'file')
    writehelpfile([outdir helpfile],obj.guidef);
end
open([outdir helpfile])

end

function pres=fieldvisibiltyparser(args)
% fields{end+1}='all';
p = inputParser;   
p.KeepUnmatched=true;
addParameter(p,'guistate','nd',@(x) any(myvalidatestring(x,{'simple','advanced','s','a','S','A','nd'})));
addParameter(p,'off',{});
addParameter(p,'on',{});

parse(p,args{:});
pres=p.Results;

end

function pard=locselectordefault
pard.selector_filelist.object=struct('String','all','Style','popupmenu');
pard.selector_filelist.position=[1,1];
pard.selector_filelist.Width=.7;

pard.selector_filter.object=struct('String',{{'layers','all u','all g','layer1','layer2','layer3'}},'Style','popupmenu'); 
pard.selector_filter.position=[1,1.7];
pard.selector_filter.Width=.7;

pard.selector_pos.object=struct('String',{{'Roi','FoV','all'}},'Style','popupmenu');
pard.selector_pos.position=[1,2.4];
pard.selector_pos.Width=.6;
% pard.filelist_ext.TooltipString=sprintf('you can define a tooltip string. \n This tip is displayed when you hover the mouse on the control');
           
end