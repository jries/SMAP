classdef Workflow<interfaces.DialogProcessor
    properties %(Access=private)
        modules={};
        graphfigure=[];
        infofigure=[];
        description='';
        makesimplegui=false;
        nirvana;
        %.module
        %.tag
        %.path:{'directory','.'
        %.inputmodules: 
        %.inputpanel
    end
    properties
        numberOfModules=0;
        startmodule
    end
    methods
        function obj=Workflow(varargin)
            obj@interfaces.DialogProcessor(varargin{:})
            if isempty(obj.handle)||~isvalid(obj.handle)
                obj.handle=figure;
                obj.handle.Visible='off';
                delete(obj.handle.Children);
            end
            if isempty(obj.nirvana)
                obj.nirvana=figure;
                obj.nirvana.Visible='off';
            end
            obj.outputParameters={'description'};
            
        end
        function initGui(obj)
            f=getParentFigure(obj.handle);
            c=uicontextmenu(f);
            obj.guihandles.modulelist.UIContextMenu=c;
            m0 = uimenu(c,'Label','info','Callback',{@listmenu_callback,obj});
            m01 = uimenu(c,'Label','add','Callback',{@listmenu_callback,obj});
            m4 = uimenu(c,'Label','add starter','Callback',{@listmenu_callback,obj});
            m02 = uimenu(c,'Label','move up','Callback',{@listmenu_callback,obj});
            m03 = uimenu(c,'Label','move down','Callback',{@listmenu_callback,obj});
            m1 = uimenu(c,'Label','remove','Callback',{@listmenu_callback,obj});
            m2 = uimenu(c,'Label','rename','Callback',{@listmenu_callback,obj});
            m3 = uimenu(c,'Label','replace','Callback',{@listmenu_callback,obj});
            m4 = uimenu(c,'Label','clear workflow','Callback',{@listmenu_callback,obj});
            
        end
        function tag=addModule(obj,varargin)
            %module, (tag), 'handle',h,'input',i,'parameters',guiPar
            p=parseinputs(varargin);
            %input parser
            %handle, tag (optional),
            %inputmodules={'module,output,moudle,output,...}, if no output,
            %output=1;
            if isa(p.module,'interfaces.WorkflowModule')|| isa(p.module,'interfaces.DialogProcessor')
                module=p.module;
            elseif iscell(p.module)
                modulepath=p.module;
                module=plugin(p.module{:});
                
                if ~(isa(module,'interfaces.WorkflowModule')|| isa(module,'interfaces.DialogProcessor'))
                    tag=[];
                    return
                end
                if ~isa(module,'interfaces.WorkflowModule') %regular plugin
                    module.pluginpath=p.module;
                    module=interfaces.plugin4workflow;
                    module.subpluginpath=p.module;
                    ph.Vpos=0;
                    p.parameters=copyfields(ph,p.parameters);
                end
                module.attachPar(obj.P);
                module.attachLocData(obj.locData);
            end
            if isempty(p.tag)
                tag=modulepath{end};
                k=1;
                newtag=tag;
                while ~isempty(obj.tag2index(newtag))
                    newtag=[tag int2str(k)];
                    k=k+1;
                end
                tag=newtag;
            else
                tag=p.tag;
            end
            if ~isempty(p.handle)
                module.handle=p.handle;
            end
            input={};
            if ~isempty(p.input)
                input=input2inputstruct(p.input);
            else
                for k=1:module.inputChannels
                    input{k}.tag=tag;
                    input{k}.outputchannel=1;
                end
            end
            modulestruc=struct('module',module,'tag',tag,'path',{modulepath},'inputmodules',{input});
            obj.numberOfModules=obj.numberOfModules+1;
            modules=obj.modules;
            modules{obj.numberOfModules}=modulestruc;
            obj.modules=modules;
            module.parent=obj;
            if ~isempty(p.parameters)
                module.setGuiAppearence(p.parameters);
            end
        end
        function module=module(obj,tag)
            if ~isnumeric(tag)
                tag=obj.tag2index(tag);
            end
            if ~isempty(tag)
                module=obj.modules{tag}.module;
            else
                module=[];
            end
        end
        function clear(obj)
            for k=1:length(obj.modules)
                delete(obj.modules{k}.module);
            end
            obj.modules={};
            obj.numberOfModules=0;
            obj.guiPar.Vpos=1;
            delete(obj.handle.Children);
            
        end
        function load(obj,fn,pGui,writeParameters,overwrite)
           if nargin<5
                overwrite=false;
            end
            if nargin<4 || isempty(writeParameters)
                writeParameters=true;
            end
            if nargin<3|| isempty(pGui)
                pGui=[];
            end
            obj.pluginpath=fn;
            obj.clear;
            obj.makeGui;
            loaded=load(fn);
            for k=1:length(loaded.modules)
                m=loaded.modules(k);
                obj.addModule(m.path,m.tag,'input',m.inputmodules);
                if overwrite
                    obj.module(m.tag).excludeFromSave={};
                end
            end
            if ~isempty(pGui) %set position parameters
                if isfield(pGui,'all')
                    pall=pGui.all;
                    pGui=myrmfield(pGui,'all');
                    if isfield(pall,'guistate')
                        switch pall.guistate
                            case {'s','simple',true,1}
%                                 obj.simplegui=true;
                                obj.makesimplegui=true;
                            otherwise
                                obj.makesimplegui=false;
                        end
                    end
                else
                    pall=[];
                end
                fn=fieldnames(pGui); 
                for k=1:length(fn)
                   
                    if isfield(pGui.(fn{k}),'handle')
%                         fn{k}
% %                         idx=obj.tag2index(fn{k})
% 
%                         obj.modules{idx}
                        modh=obj.module(fn{k});
                        if isempty(modh)
                            continue
                        end
%                         modh.handle=pGui.(fn{k}).handle;
                        modh.sethandle(pGui.(fn{k}).handle);
                    end 
                    p=pGui.(fn{k});
                    p=copyfields(p,pall);
%                     obj.module(fn{k}).setGuiAppearence(pGui.(fn{k}))
                    if ~strcmp( fn{k},'all')
                        obj.module(fn{k}).setGuiAppearence(p)
                    end
                end
            end
            obj.addAllModulesToGui;
            if writeParameters
                obj.setGuiParameters(loaded.parameters,true);
            end
            if isfield(loaded,'startmodule')
                obj.startmodule=loaded.startmodule;
            end
            if isfield(loaded,'description')
                obj.description=loaded.description;
            end
            obj.fieldvisibility;
            obj.connectModules;
             for k=1:length(obj.modules)
                 if isa(obj.modules{k}.module,'interfaces.WorkflowModule')
                 obj.modules{k}.module.updateGui; %call after loading
                 end
             end
        end
        function save(obj,fn,descriptiononly)
            if nargin<2||isempty(fn)
                fn=obj.pluginpath;
            end
            if nargin<3
                descriptiononly=false;
            end
            if descriptiononly
                load(fn);
                description=obj.description;
                save(fn,'modules','parameters','startmodule','fileformat','description');
                return
            end
            for k=1:length(obj.modules)
                mk=obj.modules{k};
                modules(k).tag=mk.tag;
                modules(k).path=mk.module.pluginpath;
                if isempty(modules(k).path)
                     modules(k).path=mk.module.subpluginpath;
                end
                modules(k).inputmodules=mk.inputmodules;
            end
            parameters=obj.getGuiParameters(true);
            startmodule=obj.startmodule;
            fileformat.name='workflow';
            description=obj.description;
            save(fn,'modules','parameters','startmodule','fileformat','description');
            
        end
        function run(obj,tag)
            if nargin<2||~ischar(tag)||isempty(tag)
                tag=obj.startmodule;
            end
            obj.connectModules;
            if isempty(obj.getPar('loc_preview'))
                obj.setPar('loc_preview',false)
            end
            obj.module(tag).clearinitialize;
            obj.module(tag).initialize;
            p=obj.module(tag).getAllParameters;
            data=interfaces.WorkflowData;
            data.eof=true;
            obj.module(tag).run(data,p);
        end
        function info=info(obj)
            
            [~,info.name]=fileparts(obj.pluginpath);
            info.description=obj.description;
        end
        
        function initialize(obj) %dummyfunction, initialize part of run
        end
        function connectModules(obj)
%             remove all connections
            for k=1:length(obj.modules)
                obj.modules{k}.module.outputModules=[];
            end
            for k=1:length(obj.modules)
                mod=obj.modules{k}.module;
                input=obj.modules{k}.inputmodules;
                for ii=1:length(input)
                    ih=input{ii};
                    if obj.tag2index(ih.tag)~=k
                        mod.setInputModule(ii,obj.module(ih.tag),ih.outputchannel);
                    elseif ~mod.isstartmodule
                        disp(['set input of ' obj.modules{k}.tag ', channel ' num2str(ii)]);
                    end
                end
            end
            if isempty(obj.startmodule)
                obj.startmodule=obj.modules{1}.tag;
            end
        end
        function setGuiPosition(obj,tag,varargin)
            set(obj.module(tag).handle,'Visible','on',varargin{:});
            %'name' 'parameter' pairs, directly put to tag.handle.
            %also for main gui
        end
        function addAllModulesToGui(obj,moduletag)
            for k=1:length(obj.modules)
                try
                    obj.addModuleToGui(obj.modules{k}.tag);
                catch err
                    disp('error loading workflow. Check workflow definition file');
                    obj.modules{k}.module
                    err.rethrow;
                end
            end
            obj.setinputlist;
        end
        function addModuleToGui(obj,moduletag)       
            idx=obj.tag2index(moduletag);
            fh=obj.guiPar.FieldHeight;
            %update list
            obj.updateModuleList;
            obj.guihandles.modulelist.Value=idx;
 
            hpos=obj.handle.Position;
            %position of panel
            poslist=obj.guihandles.modulelist.Position;
            pospanel(1)=poslist(1)+poslist(3)+4;
            pospanel(3)=hpos(3)-pospanel(1);
            pospanel(2)=poslist(2)-fh;
            pospanel(4)=poslist(4)+fh;
            
            posinput=pospanel;
            posinput(2)=pospanel(2)+pospanel(4);
            posinput(4)=2*fh+4;
            pGUI=obj.guiPar;
            pGUI.Vpos=1;
            pGUI.fontsize=obj.guiPar.fontsize-2;
            pGUI.FieldHeight=fh;
            
            thistag=obj.modules{idx}.tag;
            module=obj.modules{idx}.module;
%              module
            if isempty(module.handle)||~isvalid(module.handle)
                module.handle=uipanel(obj.handle,'Units','pixels','Position',pospanel); %later: real gui
                module.setGuiAppearence(pGUI);
              
            end
            obj.guihandles.([thistag '_modulepanel'])=module.handle;
            obj.modules{idx}.inputpanel=uipanel(obj.handle,'Units','pixels','Position',posinput); 
            obj.guihandles.([thistag '_inputpanel'])=obj.modules{idx}.inputpanel;
%             try
        
                module.simplegui=obj.makesimplegui;
              module.makeGui;  
                %shift controls if optional and makesimplegui
                if obj.makesimplegui
                    pd=module.guidef;
                    fn=fieldnames(pd);
                    for k=1:length(fn)
                        if isfield(pd.(fn{k}),'Optional')&&pd.(fn{k}).Optional
                            module.guihandles.(fn{k}).Parent=obj.nirvana;
                        end
                    end
                end
                
                tooltip='';
                
                obj.makeinputlist(thistag,obj.modules{idx}.inputpanel,module.inputChannels,module.isstartmodule,obj.modules{idx}.module.inputchanneldescription);
                obj.children.(thistag)=module;
%             catch err
%                 obj.modules(idx)=[];
%                 obj.numberOfModules=length(obj.modules);
%                 err
%             end

        end
        function pard=guidef(obj)
            pard.visualizeconnectionsbutton.object=struct('Style','pushbutton','String','Graph','Callback',@obj.graph);
            pard.visualizeconnectionsbutton.position=[2,1.4];
            pard.visualizeconnectionsbutton.Width=0.4;
            
            pard.wf_info.object=struct('Style','pushbutton','String','Info','Callback',@obj.info_callback);
            pard.wf_info.position=[2,1];
            pard.wf_info.Width=0.4;
            
%             pard.addmodulebutton.object=struct('Style','pushbutton','String','Add','Callback',@obj.add_callback);
%             pard.addmodulebutton.position=[2,1];
%             pard.addmodulebutton.Width=0.4;

%             pard.upbutton.object=struct('Style','pushbutton','String','^','Callback',{{@obj.move_callback,-1}});
%             pard.upbutton.position=[2,1.4];
%             pard.upbutton.Width=0.2;
%             pard.downbutton.object=struct('Style','pushbutton','String','v','Callback',{{@obj.move_callback,+1}});
%             pard.downbutton.position=[2,1.6];
%             pard.downbutton.Width=0.2;
            
            pard.loadbutton.object=struct('Style','pushbutton','String','Load','Callback',@obj.load_callback);
            pard.loadbutton.position=[1,1.];
            pard.loadbutton.Width=0.4;
            pard.loadbutton.object.TooltipString='press shift before clicking to load workflow without setting the saved parameters. Make sure SMAP is in focus before.';
            
            pard.savebutton.object=struct('Style','pushbutton','String','Save','Callback',@obj.save_callback);
            pard.savebutton.position=[1,1.4];
            pard.savebutton.Width=0.4;
        
            pard.modulelist.object=struct('Style','listbox','Callback',@obj.moduleselect_callback);
            pard.modulelist.position=[10.5,1];
            pard.modulelist.Height=8.5;
            pard.modulelist.Width=0.8;
            
%             pard.clearbutton.object=struct('Style','togglebutton','String','Clear','Callback',@obj.clear_callback);
%             pard.clearbutton.position=[1,1];
%             pard.clearbutton.Width=0.4;
        end
        function graph(obj,object,b)
            for k=1:length(obj.modules)
                modulenames{k}=obj.modules{k}.tag;
            end
            for k=1:length(modulenames)
                modulenames{k}=[num2str(k) '. ' modulenames{k}];
            end
            nodebox=char(ones(1,length(modulenames))*'s');
            nodecolor=char(ones(1,length(modulenames))*'b');
            
            output=[];input=[];
            for k=1:obj.numberOfModules
                mod=obj.modules{k};
                for ii=1:length(mod.inputmodules)
                    oind=obj.tag2index(mod.inputmodules{ii}.tag);
                    if isempty(oind)||oind==k
                        l=length(modulenames);
                        modulenames{l+1}='X';
                        nodebox(l+1)='r';
                        nodecolor(l+1)='r';
                        oind=l+1;
                    end
                    output(end+1)=oind;
                    input(end+1)=k;
                end
               
            end
            edgecolor=repmat('k',1,length(input));
            ntest=1:obj.numberOfModules;
            notconnected=setdiff(setdiff(ntest,output),input);
            output(end+1:end+length(notconnected))=notconnected;
            input(end+1:end+length(notconnected))=notconnected;
            edgecolor(end+1:end+length(notconnected))='w';
            if isempty(obj.graphfigure)||~isvalid(obj.graphfigure)
                obj.graphfigure=figure;
            end

            figure(obj.graphfigure);
            delete(obj.graphfigure.Children);
            plot_graph(modulenames,output,input,'-fontsize',obj.guiPar.fontsize,'-shape',nodebox,'-color',nodecolor,...
                '-edgeColor',edgecolor);
        end
        
        function showinfo(obj,edit)           
            txt=obj.description;
            if isempty(obj.infofigure)||~isvalid(obj.infofigure)
                obj.infofigure=figure;
            end
            delete(obj.infofigure.Children);
            f=obj.infofigure;
            he=uicontrol('Style','edit','units','normalized','Position',[0 0.65 1 .35],'String',txt,'Max',100,'Parent',f,'HorizontalAlignment','left');
            hplugin=uicontrol('Style','edit','units','normalized','Position',[0 0.1 1 .55],'String',txt,'Max',100,'Parent',f,'HorizontalAlignment','left');
            b2=uicontrol('Style','pushbutton','units','normalized','Position',[0.75 0 .25 .1],'String','Cancel','Callback',{@buttoncallback},'Parent',f);
            if edit
                b1=uicontrol('Style','pushbutton','units','normalized','Position',[0 0 .25 .1],'String','Save changes','Parent',f,'Callback',{@buttoncallback});
            else
%                 he.Enable='inactive';
            end    
%             hplugin.Enable='inactive';
            function buttoncallback(object,b)
                if ~strcmp(object.String,'Cancel')
                    obj.description=he.String;
                    obj.save;
                end
                close(f)
            end  
            for k=1:length(obj.modules)
                info=obj.modules{k}.module.info;
                txtp{k}=[num2str(k) '. ' info.name ': ' info.description];
            end
            hplugin.String=txtp;
        end
        
        function fieldvisibility(obj,varargin)
            for k=1:obj.numberOfModules
                mh=obj.modules{k}.module;
%                 mh.simplegui=obj.simplegui;
                mh.fieldvisibility(varargin{:});
            end
        end
    end
    methods (Access=private)      
        function load_callback(obj,a,b)
            fh=getParentFigure(obj.handle);
            modifiers = get(fh,'currentModifier');
            writeParameters=~ismember('shift',modifiers);
            if ~writeParameters
                disp('shift pressed, gui parameters not loaded')
            end
            
            [f,p]=uigetfile(['settings' filesep 'workflows' filesep '*.mat']);
            if f
                obj.load([p f],[],writeParameters)
            end
        end
        function save_callback(obj,a,b)
            [f,p]=uiputfile(['settings' filesep 'workflows' filesep '*.mat']);
            if f
                obj.save([p f])
            end
        end

        function idx=tag2index(obj,tag)
            idx=[];
            for k=1:length(obj.modules)
                if strcmp(tag,obj.modules{k}.tag)
                    idx=k;
                    break
                end
            end
        end
        function moduleselect_callback(obj,object,data)
            %add: only switch on and off if at default handle.
            for k=1:(obj.numberOfModules)
                mh=obj.module(k).handle;
                obj.modules{k}.inputpanel.Visible='off';
                mh.Visible='off';
            end
            mh=obj.module(object.Value).handle;
            obj.modules{object.Value}.inputpanel.Visible='on';
%             if mh.Parent==obj.handle
                mh.Visible='on';
%             end
        end
        function makeinputlist(obj,tag,inputhandle,inputChannels,isstartmodule,tooltips)       
            pos=inputhandle.Position;
            fs=obj.guiPar.fontsize;
            fh=obj.guiPar.FieldHeight;
            fw=(pos(3)-5)/max(1,inputChannels+double(isstartmodule));
            vpos=pos(4)-fh-obj.guiPar.Vsep-3;
            for k=1:inputChannels
                mpos(1)=(k-1)*fw+obj.guiPar.Xrim;
                mpos(2)=vpos;
                mpos(3)=fw;
                mpos(4)=fh;
                mpos2=mpos;
                mpos2(2)=mpos(2)-fh+2;
                mpos2(4)=fh-2;
                mpos2(3)=fw/2;
                mpos2(1)=mpos(1)+fw/4;
                h=uicontrol('Parent',inputhandle,'Style','popupmenu','Position',mpos,...
                    'String','1','FontSize',fs,'Callback',{@obj.input_callback,tag,k});
                
                h2=uicontrol('Parent',inputhandle,'Style','edit','Position',mpos2,...
                    'String','1','FontSize',fs,'Callback',{@obj.output_callback,tag,k});
                if length(tooltips)>=k
                    h.TooltipString=tooltips{k};
                    h2.TooltipString=tooltips{k};
                end
                obj.guihandles.([tag '_outputchannel' int2str(k)])=h2;
                obj.guihandles.([tag '_input' int2str(k)])=h;
                
            end
            if isstartmodule
                
                mpos(1)=(inputChannels)*fw+obj.guiPar.Xrim;
                mpos(2)=vpos;
                mpos(3)=fw;
                mpos(4)=fh;
                
                mpos2=mpos;
                mpos2(2)=mpos(2)-fh+2;
                mpos2(4)=fh-2;
                mpos2(3)=fw/2;
                mpos2(1)=mpos(1)+fw/4;
                
                h=uicontrol('Parent',inputhandle,'Style','pushbutton','Position',mpos,...
                    'String','Set this module as start module','FontSize',fs,'Callback',{@obj.setstart_callback,tag});
                obj.guihandles.([tag '_setstart' ])=h;
                h=uicontrol('Parent',inputhandle,'Style','pushbutton','Position',mpos2,...
                    'String','Run workflow','FontSize',fs,'Callback',{@run_callback,obj,tag});
                obj.guihandles.([tag '_run' ])=h;
            end
            
            
        end
        function setstart_callback(obj,a,b,tag)
            obj.startmodule=tag;
        end
        function input_callback(obj,object,b,tag,k)
            idx=obj.tag2index(tag);
           
                obj.modules{idx}.inputmodules{k}.tag=object.String{object.Value};
      
        end
        function output_callback(obj,object,b,tag,k)
            idx=obj.tag2index(tag);
            obj.modules{idx}.inputmodules{k}.outputchannel=str2double(object.String);
        end
        function updateModuleList(obj)       
            mn={};
            for k=1:(obj.numberOfModules)
                mn{k}=obj.modules{k}.tag;
            end
            obj.guihandles.modulelist.String=mn;
            obj.guihandles.modulelist.Value=min(obj.guihandles.modulelist.Value,length(mn));
        end
        function setinputlist(obj)
            %test if value too large, then set to default
            ml=obj.guihandles.modulelist.String;
            for k=1:obj.numberOfModules
                mh=obj.modules{k}.module;
                tag=obj.modules{k}.tag;
                mstr=ml;mstr{k}='select Input';
                 inputms=obj.modules{k}.inputmodules;
                for in=1:mh.inputChannels
                    
                    hh=obj.guihandles.([tag '_input' int2str(in)]);
                    hh.String=mstr;
%                     if length(mstr)>hh.Value
%                         hh.Value=k;
%                     end
                    
                    idx=obj.tag2index(inputms{in}.tag);
                    if isempty(idx)
                        idx=k;
                    end
                    hh.Value=idx;
                end
            end
        end
        
        function info_callback(obj,a,b)
            obj.showinfo(true)
        end
    end
end

function inputstr=input2inputstruct(input)
% 'tag',(channel)
ind=0;
for k=1:length(input)
    sth=input{k};
    if ischar(sth)
        ind=ind+1;
        inputstr{ind}.tag=sth;
    elseif isnumeric(sth)
        inputstr{ind}.outputchannel=sth;
    elseif isstruct(sth)
        ind=ind+1;
        inputstr{ind}=sth;
    end      
end
end

function pres=parseinputs(args)
%tag,mhandle,inputmoduletags
p = inputParser;   
p.KeepUnmatched=true;
addRequired(p,'module');
addOptional(p,'tag',[],@ischar);
addParameter(p,'handle',[]);
addParameter(p,'input',[]);
addParameter(p,'parameters',[]);
parse(p,args{:});
pres=p.Results;
end

function run_callback(object,b,obj,tag)
obj.run(tag);
end

function listmenu_callback(object,b,obj)
switch object.Label
    case 'remove'
        idx=obj.guihandles.modulelist.Value;
        delete(obj.modules{idx}.inputpanel);
        delete(obj.modules{idx}.module.handle);
        obj.modules(idx)=[];
        obj.numberOfModules=length(obj.modules);
        obj.updateModuleList;
        obj.setinputlist;
    case 'rename'
        idx=obj.guihandles.modulelist.Value;
        oldtag=obj.modules{idx}.tag;
        newtag=inputdlg('new tag for module','rename',1,{oldtag});
        
        if ~isempty(newtag)
            newtag=newtag{1};
            obj.modules{idx}.tag=newtag;
            obj.addModuleToGui(newtag);
            obj.updateModuleList;
            obj.setinputlist;
        end  
    case 'add starter'
         tag=obj.addModule({'WorkflowModules','Loaders','workflowstarter'});
         obj.addModuleToGui(tag);
         obj.modules=obj.modules([end 1:end-1]);
         obj.updateModuleList;
         obj.setinputlist;
    case 'replace'
         idx=obj.guihandles.modulelist.Value;
         oldmodule=obj.modules{idx};
%          obj.add_callback(0,0);
         listmenu_callback(struct('Label','add'),b,obj)
         obj.modules{idx}=obj.modules{end};
         obj.modules(end)=[];
         obj.numberOfModules=length(obj.modules);
         delete(oldmodule.inputpanel);
        delete(oldmodule.module.handle);
         obj.updateModuleList;
         obj.setinputlist;
         obj.guihandles.modulelist.Value=idx;
    case 'info'
         idx=obj.guihandles.modulelist.Value;
         module=obj.modules{idx};
         info=module.module.info;
         desc=info.description;
         if iscell(desc)
             desc=desc{1};
         end
         mymsgbox(desc,info.name);
%         txth= findobj(hm,'Type','text','-depth',3);
%         txth.FontSize=14;
    case 'add'
        plugins=plugin;
%         field=browsefields(plugins,'WorkflowModules');
        field=browsefields(plugins,'WorkflowModules',[],[],true);
        if ~isempty(field)
          tag=obj.addModule(field2cell(field));
          if ~isempty(tag)
          obj.addModuleToGui(tag);
          obj.setinputlist;
          else
              disp('could not add module. It is not a WorkflowModule or DialogProcessor.')
          end

        end
    case 'move up'
        dir=-1;
            select=obj.guihandles.modulelist.Value;
            lenlist=length(obj.guihandles.modulelist.String);
            newpos=select+dir;
            if newpos<=lenlist&&newpos>0
                obj.modules([select,newpos])=obj.modules([newpos,select]);
            end
            obj.guihandles.modulelist.Value=newpos;
            obj.updateModuleList;
            
            obj.setinputlist;
    case 'move down'
        dir =1;
            select=obj.guihandles.modulelist.Value;
            lenlist=length(obj.guihandles.modulelist.String);
            newpos=select+dir;
            if newpos<=lenlist&&newpos>0
                obj.modules([select,newpos])=obj.modules([newpos,select]);
            end
            obj.guihandles.modulelist.Value=newpos;
            obj.updateModuleList;
            
            obj.setinputlist;
    case 'clear workflow'
            obj.clear;
            obj.makeGui;
        
end
end


