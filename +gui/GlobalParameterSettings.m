classdef GlobalParameterSettings < interfaces.GuiModuleInterface
    properties

    end
    methods
        function obj=GlobalParameterSettings(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:});
            obj.makeGui;
        end
        function makeGui(obj)
            if isempty(obj.handle)||~isvalid(obj.handle)
                obj.handle=figure('MenuBar','none','Toolbar','none');
            else
                delete(obj.handle.Children);
            end
            numlines=8;
            ap=obj.P.globalSettings;
            fn=fieldnames(ap);
            htabgroup=uitabgroup(obj.handle);
            htabgroup.Position(4)=.8;
            htabgroup.Position(2)=.2;
            for k=1:length(fn)
                phere=ap.(fn{k});
                category=phere.category;
                cs=['category_' category];
                if ~isfield(obj.guihandles,cs)
                    
                    obj.guihandles.(cs)=uitab(htabgroup,'Title',category);
                    positions.(cs)=0;
                end
                positions.(cs)=positions.(cs)+2;
                poshere=[0 1-positions.(cs)/numlines 1 .75*1/numlines];
                os=[category '_' fn{k}];
                if strcmp(phere.object.Style,'file')||strcmp(phere.object.Style,'dir')
                    poshere(3)=.75;
                    obj.guihandles.(os)=uicontrol('Parent',obj.guihandles.(cs),'Units','normalized','Style','edit',...
                        'FontSize',obj.guiPar.fontsize,'String',phere.object.String,'Position',poshere);
                    poshere(1)=.75; poshere(3)=.25;
                    b=uicontrol('Parent',obj.guihandles.(cs),'Units','normalized','Style','pushbutton',...
                        'FontSize',obj.guiPar.fontsize,'String','browse','Position',poshere);
                    b.Callback={@browsecallback,obj.guihandles.(os),phere.object.Style};
                elseif strcmp(phere.object.Style,'saveparameter')
                    if isfield(phere.object,'String')
                        str=phere.object.String;
                    else
                        str=save;
                    end
                    
                    b=uicontrol('Parent',obj.guihandles.(cs),'Units','normalized','Style','pushbutton',...
                        'FontSize',obj.guiPar.fontsize,'String',str,'Position',poshere,'Callback',{@obj.save_callback,fn{k}});
                    obj.guihandles.(os)=b;
                else
                    obj.guihandles.(os)=uicontrol('Parent',obj.guihandles.(cs),'Units','normalized','Style',phere.object.Style,...
                        'FontSize',obj.guiPar.fontsize,'String',phere.object.String,'Position',poshere,'Value',phere.object.Value);
                   
                    %any other controlo
                end
                poshere=[0 1-(positions.(cs)-1+.1)/numlines 1 .6*1/numlines];
                obj.guihandles.([os 't'])=uicontrol('Parent',obj.guihandles.(cs),'Units','normalized','Style','text',...
                        'FontSize',obj.guiPar.fontsize,'String',phere.description,'Position',poshere,'HorizontalAlignment','left');
            end
            poshere=[0,.11,1,.09];
            uicontrol('Parent',obj.handle,'Style','text','String','Some changes might take effect only after restart',...
                'FontSize',obj.guiPar.fontsize*.8,'Units','normalized','Position',poshere);
            poshere=[0.65,0,.25,.11];
            obj.guihandles.exitSettings=uicontrol('Parent',obj.handle,'Units','normalized','Style','pushbutton',...
                        'FontSize',obj.guiPar.fontsize,'String','Save and exit','Position',poshere);
            obj.guihandles.exitSettings.Callback=@obj.exitSave;
            poshere=[0.35,0,.25,.11];
            obj.guihandles.exitCancel=uicontrol('Parent',obj.handle,'Units','normalized','Style','pushbutton',...
                        'FontSize',obj.guiPar.fontsize,'String','Cancel','Position',poshere);
            obj.guihandles.exitCancel.Callback=@obj.exitCancel;
            
            
        end
        function exitSave(obj,a,b)
            ap=obj.P.globalSettings;
            fn=fieldnames(ap);
            for k=1:length(fn)
                phere=ap.(fn{k});
                hhere=obj.guihandles.([phere.category '_' fn{k}]);
                ap.(fn{k}).object=copyfields(ap.(fn{k}).object,hhere,{'String','Value'});
            end
            obj.P.globalSettings=ap;
            obj.saveGlobalSettings;
            
            delete(obj.handle)
            delete(obj)
        end

        function exitCancel(obj,a,b)
            delete(obj.handle)
            delete(obj)
        end  
        
        function save_callback(obj,a,b,field)
            par=obj.getPar(field);
            
            [file,path]=uiputfile([pwd '/settings/' field '.txt']);
            if file
                writestruct([path file],par);
            end
        end
    end
end

function browsecallback(a,b,hedit,style)
if strcmp(style,'file')
    strstart=hedit.String;
    if isempty(strstart)
        strstart='settings/*.txt';
    end
    [f,p]=uigetfile(strstart);
    if p
        p=makerelativetopwr(p);
        hedit.String=[p f];
    end
else
    strstart=hedit.String;
    p=uigetdir(strstart);
    if p
        p=makerelativetopwr(p);
        hedit.String=p;
    end
end

end