classdef GuiParameterInterface<interfaces.ParameterInterface
    %provides functionality for global parameter handling related to GUI
    %parameters such as synchronization
    properties
       syncParameters %cell array of {field, handle of uicontrol, syncmode ='String|Value|otherproperty'}. Set in guidef.
      
    end
    methods
        function addSynchronization(obj,field,handle,syncmode,changecallback)
             %add synchronization to global parameters
            %field: globalParameterName
            %handle: handle to guiobject to be synchronized
            %addSynchronization(parametername,handle to synchronize to,syncmode,changecallback)
            %syncmode='String|Value|otherproperty': what to synchronize
            %changecallback: function handle to function which is called
            %when parameter is changed. This is similar to events and
            %listeners.
            %Instead of calling addSynchronization you can define in
            %obj.guidef: pard.syncParameters={{'globalParameterName','guiobject2',{'String','Value'},{@changecallback,obj}},...};
            if nargin<5
                changecallback={};
            end
            if ~iscell(changecallback)
                changecallback={changecallback};
            end
            hstruc.changecallback=changecallback;
            if ~isempty(handle) %|| ~isempty(changecallback)
                hstruc.isGuiPar=true;
            else
                hstruc.isGuiPar=false;
            end
            hstruc.synchronize=true;
            hstruc.handle=handle;
            hstruc.syncmode=syncmode;
            hstruc.obj=obj;
            if ~isempty(handle) %hstruc.isGuiPar
                handle.Callback={@obj.field_callback,field,handle.Callback};
            end
            
            hstruc.content=obj.handle2value(handle);
            if ~isfield(obj.P.par,field)
                obj.P.par.(field)=hstruc;
            else
                lf=length(obj.P.par.(field));
                posnew=lf+1;
                outind=false(lf+1,1);
                for k=1:lf
                    phx=obj.P.par.(field)(k);
                    if ~isvalid(phx.obj)
                        outind(k)=true;
                    end
                    if phx.obj==hstruc.obj &&((isempty(handle) && isempty(phx.handle) )|| (~isempty(handle) && phx.handle==handle))
                        posnew=k;
                        break
                    end
                end
                obj.P.par.(field)(posnew)=copyfields(obj.P.par.(field)(1),hstruc);
%                 obj.P.par.(field)(posnew).content=obj.P.par.(field)(1).content; %xxx 
%                 obj.P.par.(field)(posnew)=hstruc;
                obj.P.par.(field)(outind)=[];
            end
            isgp=[obj.P.par.(field).isGuiPar];
            guip=any(isgp);
            indchange=find(isgp~=guip);
%             for k=1:length(obj.P.par.(field))
            for k=1:length(indchange)
                obj.P.par.(field)(indchange(k)).isGuiPar=guip;
            end  
            found=false;
            for k=1:length(obj.syncParameters) %dont attach many times
                if (strcmp(obj.syncParameters{k}{1},field))
                    found=true;
                    break;
                end
            end
            if ~found
                obj.syncParameters{end+1}={field,handle,syncmode,changecallback};
            end
            
            %initialize: not neeede, done before?
%             v=obj.getPar(field);
%             if ~isempty(v)
%                 hh=obj.value2handle(v,handle);
%                 obj.setfields(field,hh);
%             end
        
        end
        function setPar(obj,field,varargin)
            %writes global parameter. If synchronized: update uicontrols
            %and execute callbacks
            %setPar(parameter,newvalue,'Property','Value')
            %parameter: name of global parameter
            %newvalue: string or numeric or handle structure or
            %'Property','Value' pairs for handle.
            %Properties:
            %'layer',layer: access layer parameters. layerN_ used as prefix to
            %parametername
            %if  not GUI parameter: calls interfaces.ParameterInterface

            ind=find(strcmp(varargin,'layer'));
            if ~isempty(ind)
                prefix=['layer' num2str(varargin{ind+1}) '_'];
                varargin(ind:ind+1)=[];
                field=[prefix field];
            end
            par=obj.P.par;
            if isfield(par,field)
                if par.(field)(1).isGuiPar
                    if length(varargin)==1
                        for k=1:length(par.(field))
                            if ishandle(par.(field)(k).handle)
                                break
                            end
                        end
                        
                        handle=obj.value2handle(varargin{1},par.(field)(k).handle);
                    else
                        for k=1:2:length(varargin)
                            handle.(varargin{k+1})=varargin{k};
                        end
                    end
                    obj.setfields(field,handle);
                else
                    obj.setfields(field,varargin{1});
                end
                if any([par.(field)(:).synchronize])
                    obj.changecallback(field)
                end
            else
                setPar@interfaces.ParameterInterface(obj,field,varargin{1})
            end
        end
        function value=getPar(obj,field,varargin)
%             global SMAPparameters
            %returns global parameter. 
            %value=getPar(parameter,'Property','Value')
            %parameter: name of global parameter
            %'Property','Value' pairs for handle.
            %Properties:
            %'layer',layer: access layer parameters. layerN_ used as prefix to
            %parametername
            %if  not GUI parameter: calls interfaces.ParameterInterface
            ind=find(strcmp(varargin,'layer'));
            if ~isempty(ind)
                prefix=['layer' int2str(varargin{ind+1}) '_'];
                varargin(ind:ind+1)=[];
                field=[prefix field];
            end
%             po=obj.P;
%             par=po.par;
%             par=SMAPparameters;
            par=obj.P.par;
            if isfield(par,field)
%                 tested=true;
                if par.(field)(1).isGuiPar
                    if ~isempty(varargin)&&~isempty(varargin{1})
                        for k=1:length(par.(field))
                            value{k}=par.(field)(k).handle.(varargin{1});
                        end
                    else
                        %find valid handle
                        %notempty
                        hhere=[];
                        content=[];
                        for k=1:length(par.(field))
                            if ~isempty(par.(field)(k).handle)&&isvalid(par.(field)(k).handle)
                                hhere=par.(field)(k);
                                break
                            end
                            content=par.(field)(k).content;
                        end
                        if ~isempty(hhere)
                            value=obj.handle2value(hhere.handle);
                        else
                            value=content;
                        end
                    end
                else
                    value=par.(field).content; %from parameterInterface
%                     value=getPar@interfaces.ParameterInterface(obj,field,tested);
                end
            else
                value=[];
               
            end      
        end
        
        function updateObjectParameters(obj)
            %writes obj.outputParameters which are uicontrol objects in gui to global parameters
            updateObjectParameters@interfaces.ParameterInterface(obj);
            pgui=obj.getGuiParameters;
            for k=1:length(obj.outputParameters)
                if isfield(pgui,obj.outputParameters{k})
                    obj.setPar(obj.outputParameters{k},pgui.(obj.outputParameters{k}));
                end
            end
        end   
        
        function par=handle2value(obj,hfn,onlyedit)
            %parses uicontrol handles. 
            %TODO: remove from here, use as function
            %par=handle2value(handle)
            %'edit','text': String or str2num(String) if numeric
            % binary controls: value
            %popupmenu, listbox: structure with h.String,h.Value,
            %.selection=h.String{h.Value}
            %Table: Data    
            if nargin <3
                onlyedit=false;
            end
            par=[];
            if (isa(hfn,'matlab.ui.control.UIControl')&&isvalid(hfn)) || (isstruct(hfn)&&isfield(hfn,'Style')&&(isfield(hfn,'String')||isfield(hfn,'Value')))
                style=hfn.Style;
                switch style
                    case {'edit','text','file','dir'}
                        if onlyedit && ~strcmp(style,'edit')
                            par=[];
                            return
                        end
                        
                        st=hfn.String;
                        if isempty(st)
                            par='';
                            return
                        end
%                         v=str2num(st);
%                         testv=str2double(st(1));
                        if iscell(st)|| (~((st(1)>='0'&&st(1)<='9') || st(1)=='-' || st(1)=='I' || st(1)=='i' || st(1)=='.'))
                            par=hfn.String;
                        else
                            if contains(st,',') || contains(st,' ') || contains(st,':') || contains(st,';')
                                v=str2num(st);
                            else
                                v=str2double(st);
%                                 if isnan(v)
%                                     st
%                                     v=str2num(st); 
%                                 end
                            end
  
                            if isempty(v) || any(isnan(v))
                                 par=hfn.String;
                            else     
                                par=1*v;
                            end
                        end
                    case {'togglebutton','checkbox','slider', 'pushbutton','radiobutton'}
                        if onlyedit && strcmp(style,'pushbutton')
                            par=[];
                            return
                        end
                        v=double(hfn.Value);
                        par=1*v;
                    case {'popupmenu','listbox'}    
                            par.String=hfn.String; 
                            par.Value=hfn.Value;
                            sstring=hfn.String;
                        if iscell(sstring)
                            if hfn.Value>0
                                par.selection=sstring{hfn.Value};
                            else
                                par.selection='';
                            end
                        else
                            s=size(sstring);
                            if s(1)>1
                                ss=sstring(min(hfn.Value,s(1)),:);
                                par.selection=strtrim(ss);
                            else
                                par.selection=hfn.String;
                            end
                        end
                        if onlyedit
                            par=myrmfield(par,'String');
                        end
                end
            elseif isa(hfn,'matlab.ui.control.Table')
                par.Data=hfn.Data;
            else
% %                 hfn
            end
        end
        function handle=value2handle(obj,v,hin)
             %convert value to handle structure. 
            %TODO: remove from here, use as function
            %handle=value2handle(value,hin)
            %hin: sample handle structure used as template
            %'edit','text': String or str2num(String) if numeric
            % binary controls: value
            %popupmenu, listbox: structure with h.String,h.Value,
            %.selection=h.String{h.Value}
            %Table: Data
            handle=[];
            if (isa(hin,'matlab.ui.control.UIControl')&&isvalid(hin))|| (isstruct(hin)&&isfield(hin,'Style')&&(isfield(hin,'String')||isfield(hin,'Value')))
            switch hin.Style
                case {'edit','file'} %file:own creation for global settings
                    if isnumeric(v)
                        if length(v)>1
                            fs='%g,';
                        elseif abs(v)>100
                            fs='%6.0f';
                        else
                            fs=5;
                        end
                        strg=num2str(v,fs);
                        if length(v)>1
                            strg(end)=[];
                        end
                        handle.String=strg;
                    else
                        handle.String=v;
                    end
                case 'text'
                    if ~iscell(v)
                    handle.String=num2str(v);
                    else
                        handle.String=(v);
                    end
                case {'pushbutton','checkbox','togglebutton'}
                    if isnumeric(v)||islogical(v)
                        handle.Value=v;
                    elseif isstruct(v)
                        handle=v;
                    elseif ischar(v)
                        handle.String='v';
                    end

                case {'popupmenu','listbox'}
                    if isfield(v,'String')
                    handle.String=v.String;
                    handle.Value=v.Value;
                    else
                        handle.String=v;
                        handle.Value=1;
                    end
                case 'slider'
                    handle.Value=max(min(v,hin.Max),hin.Min);
            end
            elseif isa(hin,'matlab.ui.control.Table')
                handle.Data=v.Data;
            end

        end
        function createGlobalSetting(obj,field,category,description,structure)
            global SMAP_globalsettings
            if ~isfield(obj.P.globalSettings,field) %don't overwrite current settings
                obj.P.globalSettings.(field).object=structure;
    %             obj.P.globalSettings.(field).name=name;
                obj.P.globalSettings.(field).category=category;
                obj.P.globalSettings.(field).description=description;
                SMAP_globalsettings=obj.P.globalSettings;
                obj.saveGlobalSettings;
                
            end
            %save to global parameters
        end
        function setGlobalSetting(obj,field,value)
            
            h=obj.value2handle(value,obj.P.globalSettings.(field).object);
            obj.P.globalSettings.(field).object=copyfields(obj.P.globalSettings.(field).object,h);
            obj.saveGlobalSettings;
            %save to global parameters
        end
        function p=getGlobalSetting(obj,field)
            p=obj.handle2value(obj.P.globalSettings.(field).object);
        end
        function saveGlobalSettings(obj)
            global SMAP_globalsettings
            file=[obj.getPar('maindirectory') filesep obj.P.globalSettingsFile];
            writestruct(file,obj.P.globalSettings);
            SMAP_globalsettings=obj.P.globalSettings;
        end
        function loadGlobalSettings(obj)
            global SMAP_globalsettings
            file=[obj.getPar('maindirectory') filesep obj.P.globalSettingsFile];
            obj.P.loadGlobalSettings(file);
            SMAP_globalsettings=obj.P.globalSettings;
        end

    end
    
    methods (Access=private)      
        function setfields(obj,field,handle)
%             global SMAPparameters 
            hstruc=obj.P.par.(field);
            if hstruc(1).isGuiPar
                for k=1:length(hstruc)
                    if ~isempty(hstruc(k).handle)&&isvalid(hstruc(k).handle)
                        syncmode=hstruc(k).syncmode;
                        if ~iscell(syncmode)
                            syncmode={syncmode};
                        end
                        for l=1:length(syncmode)
                            
                            if ~isempty(handle) && (isprop(handle,syncmode{l})||isfield(handle,syncmode{l}))
                                hstruc(k).handle.(syncmode{l})=handle.(syncmode{l});
                            end
                            switch hstruc(k).handle.Style
                                case {'popupmenu','listbox'}
                                    hstruc(k).handle.Value=min(hstruc(k).handle.Value,length(hstruc(k).handle.String));
                            end
                        end
                    end
                end
            else
                for k=1:length(hstruc)
                    hstruc(k).content=handle;
                end

%             SMAPparameters.(field)=hstruc; %Hack to improve performance. its the same as:
                obj.P.par.(field)=hstruc;
            end
        end
        
        function changecallback(obj,field)
            hstruc=obj.P.par.(field);
            badstruc=false(length(hstruc));
            for k=1:length(hstruc)
                if ~isempty(hstruc(k).changecallback)
                    cb=hstruc(k).changecallback;
                    callobj=hstruc(k).obj;
                    if ~isvalid(callobj) %object has been deleted
                        badstruc(k)=true;
                    else
                        if ~isa(obj,class(hstruc(k).obj))&& isvalid(callobj.handle) %or ==? of identical object, in case of copy??
                            feval(cb{:})
                        end
                    end
                end
            end
            if any(badstruc)
                hstruc(badstruc)=[];
                obj.P.par.(field)=hstruc;
            end
        end

        function field_callback(obj,handle,event,field,oldcallback)
            obj.setfields(field,handle);
            %determine callbacks
            docallbacks=true;
            sm=obj.P.par.(field)(1).syncmode;
            if ~isempty(sm)   
                if strcmp(handle.Style,'listbox')&&any(strcmp(sm,'String'))
                    docallbacks=false;
                end
            end
            if obj.P.par.(field)(1).synchronize&&docallbacks
                obj.changecallback(field)
            end
            if ~isempty(oldcallback)
                if ~iscell(oldcallback)
                    oldcallback={oldcallback};
                end
                feval(oldcallback{1},handle,event,oldcallback{2:end})
            end
        end
    end
end