classdef DialogProcessor<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    %extends GuiModuleInterface with a results window, a Process button adn
    %Info button. This is the default class for Analyzer and Processor
    %modules
    properties
        resultstabgroup;   %handle to results  
        processorgui=true; %switch. true if process button etc are to be rendered. false if called externally (for workflow)
        showresults=false; % defined state for results
        history=false;
         parent; %encapsulating object, e.g. used to call methods in plugin4workflow from plugin directly
         undo=false;
         showloadsavepar=false;
%         moduleinfo;
    end
    properties (SetAccess = private, GetAccess = private)
       resultshandle
    end
    methods
        function obj=DialogProcessor(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})  
        end
        
        function makeGui(obj,guidef)
            %calls makeGui@GuiModuleInterface, additionally provides info
            %and process button and checkbox to show results
            if nargin==1
                guidef=obj.guidef;
            end
            if obj.processorgui
                obj.guiPar.fontsize=obj.guiPar.fontsize-1;
                obj.guiPar.FieldHeight=obj.guiPar.FieldHeight-2;
            end
     
            anyoptional=makeGui@interfaces.GuiModuleInterface(obj,guidef);
            
            if obj.guiselector.show || anyoptional
                posh=obj.handle.Position;
                if isempty(obj.guiselector.position)
                    pos(1:2)=0*posh(1:2)+posh(3:4)-[17,21];
                    pos(3:4)=[15,20];
                else
                    pos=obj.guiselector.position;
                end
                obj.guihandles.simplegui=uicontrol(obj.handle,'Style','togglebutton','String','A','Position',pos,'Callback',@obj.simplegui_callback);
                obj.guihandles.simplegui.TooltipString='toggle between simple and advanced GUI';
                obj.fieldvisibility('guistate',obj.simplegui);
            end
            
            if obj.processorgui && ~isempty(obj.handle)
                vis='on';
            else
                vis='off';
            end
            if isvalid(obj.handle)
                hpos=obj.handle.Position;
                vrim=obj.guiPar.Vrim;
                obj.guihandles.showresults=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2+40, hpos(4)-vrim+20,120,20],...
                    'FontSize',obj.guiPar.fontsize,'Style','checkbox', 'String', 'Show results',...
                    'Value',obj.showresults,'Callback',{@showresults_callback,obj},'Visible',vis);
                obj.guihandles.processgo_b=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*4-100, hpos(4)-vrim+20,100,55],...
                    'Style','pushbutton','String','Run','FontSize',obj.guiPar.fontsize*1.5,'Callback',{@processgo_callback,obj},'Visible',vis);
                obj.guihandles.info=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2+40, hpos(4)-vrim+50,100,25],...
                    'Style','pushbutton','String','Info','FontSize',obj.guiPar.fontsize,'Callback',{@info_callback,obj},'Visible',vis,...
                    'BackgroundColor',[0.7,1,0.9]);
                if obj.showloadsavepar %show parameter load and save functions
                    obj.guihandles.saveguipar=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2-20, hpos(4)-vrim+52,40,20],...
                        'Style','pushbutton','String','save','FontSize',obj.guiPar.fontsize*1.,'Callback',{@saveguipar_callback,obj},'Visible',vis);
                     obj.guihandles.loadguipar=uicontrol(obj.handle,'Position',[obj.guiPar.FieldWidth*2-20, hpos(4)-vrim+22,40,20],...
                        'Style','pushbutton','String','load','FontSize',obj.guiPar.fontsize*1.,'Callback',{@loadguipar_callback,obj},'Visible',vis);            
                end
            end
            obj.initGuiFinal;
        end
        function initGuiFinal(obj)
        end
        function setvisibility(obj,name)
            %shows and hides the GUI. called from module selector
            %setvisibility(visible) visible='on'/'off'
            if isvalid(obj) && isvalid(obj.handle)&&isvalid(obj.handle.Parent)&&~isa(obj.handle.Parent,'matlab.ui.Figure')
            set(obj.handle,'Visible',name);
            end
        end
        
        function makeResultsWindow(obj)
            %creates a window with results tabs
            obj.resultshandle=figure;
            obj.resultshandle.Visible='off';
            obj.resultshandle.Renderer='painters';
            
            %make convert to grid
            hToolbar=findall(obj.resultshandle, 'type', 'uitoolbar');
            pushtool=uipushtool(hToolbar);
            img=ones(16,16,3);
            img(6,:,:)=0;img(11,:,:)=0;img(:,6,:)=0;img(:,11,:)=0;
            pushtool.CData=img;
            pushtool.ClickedCallback=@makegridoutput;
            pushtool.Tooltip='Create new figure with all tabs arranged next to each other';
            
            htab=uitabgroup(obj.resultshandle);
            obj.guihandles.resultstabgroup=htab;
            obj.resultstabgroup=obj.guihandles.resultstabgroup;
        end
        function ax=initaxis(obj,varargin)
            %initializes axis in results window
            if isempty(obj.resultstabgroup) || ~isvalid(obj.resultstabgroup)
                obj.makeResultsWindow;
                obj.resultshandle.Visible='on';
            end
            ax=initaxis(obj.resultstabgroup,varargin{:});
        end
        function results=processgo(obj)
            %provides external access to run module (usually via process
            %button)
            results=processgo_callback(0,0,obj);
        end
        function addhistory(obj)
            p.parameters=obj.getGuiParameters(true,true);
            p.name=class(obj);
            obj.locData.addhistory(p);
        end      
        function simplegui_callback(obj,~,~)
            simplegui=obj.getSingleGuiParameter('simplegui');
            obj.fieldvisibility('guistate',simplegui);
        end
        
%         function setguistate(obj,state)
%             if nargin>1
%                 obj.simplegui=state;
%             end
%             if ~obj.simplegui
%                 state='on';
%             else
%                 state='off';
%             end
%             pard=obj.guidef;
%             fn=fieldnames(pard);
%             for k=1:length(fn)
%                 if isfield(pard.(fn{k}),'Optional')&&pard.(fn{k}).Optional==true
%                     obj.guihandles.(fn{k}).Visible=state;
%                 end
%             end
%         end
% 
    end
    methods (Access=private)
%         function outhandle=addresultstab(obj,name)
%             %adds a tab to results figure. 
%             outhandle=uitab(obj.guihandles.resultstabgroup,'Title',name);
%         end
    end
end

function results=processgo_callback(~,~,obj)
% notify(obj.locData,'undo',recgui.simpleEvent('backup'));
results=[];
obj.status(['executing ' class(obj)])
drawnow;
if isempty(obj.resultshandle)||~isvalid(obj.resultshandle)
    obj.makeResultsWindow;
end


p=obj.getAllParameters;

 p.obj=obj;
p.resultstabgroup=obj.guihandles.resultstabgroup;
if obj.processorgui
obj.resultshandle.Visible=onoff(p.showresults);
end
if isempty(obj.locData.loc) && obj.getGlobalSetting('advancedoutput')
    warning('no localization data present')
end

if obj.undo
    obj.setPar('undoModule',obj.info.name);
    notify(obj.P,'backup4undo');
end
results=obj.run(p);

if ~isempty(results)
    obj.setAutoResults(obj.pluginpath,results);
    if isfield(results,'clipboard')
        cl=results.clipboard;
        if ~iscell(cl)
            cl={cl};
        end
        
        for k=1:length(cl)
            if isnumeric(cl{k})
                cl{k}=num2str(cl{k});
            end
        end
        ct=sprintf('%s\t',cl{:});
        clipboard('copy',ct); 
    end
end
if obj.history
    obj.addhistory;
end
obj.resultshandle.Visible=onoff(p.showresults);
if ~isfield(results,'error')||isempty(results.error)
    obj.status([class(obj) ' finished'])
else
    obj.status(['ERROR in ' class(obj) '. ' results.error])
    obj.setPar('errorindicator', ['ERROR in ' class(obj) '. ' results.error]);
end
end


function showresults_callback(object,~,obj)
if isempty(obj.resultshandle)||~isvalid(obj.resultshandle)
    obj.makeResultsWindow;
end
switch object.Value
    case 1
        state='on';
    case 0 
        state='off';
end
set(obj.resultshandle,'Visible',state)
end

function info_callback(a,b,obj)
obj.showinfo;
% warnid='MATLAB:strrep:InvalidInputType';
% warnstruct=warning('off',warnid);
% 
% obj.guihandles.showresults.Value=1;
% showresults_callback(obj.guihandles.showresults,0,obj)
% ax=obj.initaxis('Info');
% hp=ax.Parent;
% hp.BackgroundColor='w';
% delete(ax);
% 
%  htxt=annotation(hp,'textbox',[0.03,0.03,.9,.93],...
%      'FontSize',obj.guiPar.fontsize,'HorizontalAlignment','left','BackgroundColor','w','FitBoxToText','off','EdgeColor','w');
% %  htxt=uicontrol(hp,'Style','text','Units','normalized','Position',[0,0,.9,1],...
% %      'FontSize',obj.guiPar.fontsize,'HorizontalAlignment','left','Max',100); 
% td=obj.info.description;
%  if ~iscell(td)
%   txt=strrep(td,9,' ');
% % txt=strrep(td,10,13);
%  else
%      txt=td;
%  end
%  htxt.String=txt;
% %   htxt.Position=[0 0 1 1];
%   warning(warnstruct);
%   h=uicontrol(hp,'Style','pushbutton','Units','normalized','Position',[0.95,0.95,.05,.05],'String','Edit','Callback',{@edit_callback,obj});
end

function makegridoutput(a,b)
fig=a.Parent.Parent;
numtabs=length(fig.Children.Children);
numrows=ceil(sqrt(numtabs));
flattenTabs(fig,numrows);
end

function loadguipar_callback(a,b,obj)
   pfad=fileparts(obj.getPar('lastSMLFile'));
   outf=[pfad filesep obj.pluginpath{end} '*_parameters.mat'];   
   [f,pfad]=uigetfile(outf);
   if f
       l=load([pfad f]);
   end  
   obj.setGuiParameters(l.par,true);
end

function saveguipar_callback(a,b,obj)
   par=obj.getGuiParameters(true);
   pfad=fileparts(obj.getPar('lastSMLFile'));
   outf=[pfad filesep obj.pluginpath{end} '_parameters.mat'];
   [f,pfad]=uiputfile(outf);
   if f
       save([pfad f], 'par')
   end    
end
