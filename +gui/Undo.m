classdef Undo< interfaces.GuiModuleInterface & interfaces.LocDataInterface
%     backs up localization data before executing plugin, restores it
    properties
        locDataOld
        undoModule
    end
    methods
        function obj=Undo(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})          
        end
        function makeGui(obj)
            obj.guihandles.undobutton=uicontrol(obj.handle,'Style','pushbutton','String','Undo','Units','normalized',...
                'Position',[0.8,0.002,.07,.03],'Callback',@obj.undo_callback);
             obj.guihandles.undobutton.Units='pixels';
             obj.guihandles.undobutton.Position(4)=28;
             
             obj.guihandles.undobutton.TooltipString='Many Plugins backup the localization data before modifying them. The old data can then be restored with Undo';
            addlistener(obj.P,'backup4undo',@obj.backup);
        end
        function backup(obj,event,data) 
            obj.undoModule=obj.getPar('undoModule');
            obj.locDataOld=obj.locData.copy;
        end
        function undo_callback(obj,a,b)
            if isempty(obj.locDataOld)
%                 disp('nothing stored for undo')
                obj.status(['nothing stored for undo ' ])
                return
            end        
            temp=obj.locData.copy;
            obj.locData.loc=obj.locDataOld.loc;
            obj.locData.grouploc=obj.locDataOld.grouploc;
            obj.locData.files=obj.locDataOld.files;
            obj.setPar('locFields',fieldnames(obj.locData.loc));
            
            fl=({obj.locData.files.file.name});
            fls=fl;
            for k=1:length(fl)
                flh=strrep(fls{k},'\','/');
                [~ ,fls{k}]=fileparts(flh);
            end

            obj.setPar('filelist_long',fl);
            obj.setPar('filelist_short',fls);
            fsx=[{'layer','all'} fls];
            obj.setPar('filelist_short_ext',fsx,'String');
            obj.locDataOld=temp;
            
            obj.status(['undo performed: ' obj.undoModule])
            obj.undoModule=['undo of ' obj.undoModule];
            
        end
    end
end

