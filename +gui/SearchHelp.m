classdef SearchHelp<interfaces.GuiModuleInterface
% 
    properties
       figure
%        functionlist
       alltxt
       helpdir
       annotation
       searchind
    end
    
    methods
        function obj=SearchHelp(varargin)
%             makeGui(obj)
        end
        
        function makeGui(obj)
            makeGui(obj)
        end
    end
    
end


function makeGui(obj)
if isempty(obj.figure) || ~isvalid(obj.figure)
    obj.figure=figure('MenuBar','none','Name','Search in Help','Units','normalized');
    obj.handle=obj.figure;
end
hf=obj.figure;
h.searchstring=uicontrol(hf, 'Style','edit','Units','normalized','Position',[0.02,0.93,0.75,0.05]);
h.searchbutton=uicontrol(hf, 'Style','pushbutton','String','Search','Units','normalized','Position',[0.8,0.93,0.18,0.05],...
    'Callback',{@search_callback,obj});
h.resultslist=uicontrol(hf,'Style','listbox','Units','normalized','Position',[0.02,0.02,0.3,0.9],'Callback',{@list_callback,obj});
h.resultspanel=uipanel(hf,'Units','normalized','Position',[0.32,0.22,0.66,0.7],'BackgroundColor','w');
h.resultstt=uicontrol(hf,'Style','edit','Units','normalized','Position',[0.32,0.02,0.66,0.2],'BackgroundColor','w','HorizontalAlignment','left');

obj.guihandles=h;
end

function search_callback(a,b,obj)
alltxt=obj.getPar('allhelpfiles');
if isempty(alltxt)
    alltxt=makehelptxt(obj);
    obj.setPar('allhelpfiles',alltxt);
    obj.alltxt=alltxt;
end
searchstring=obj.guihandles.searchstring.String;
ind=find(contains(alltxt.text,searchstring));
obj.searchind=ind;
obj.guihandles.resultslist.String=onlynames(alltxt.filenames(ind));
% obj.functionlist=alltxt.filenames(ind);
% obj.functiontxt=alltxt.text(ind);
list_callback(obj.guihandles.resultslist,0,obj);
end

function list_callback(a,b,obj)
select=a.Value;
allindh=obj.searchind(select);
delete(obj.annotation)
ss=obj.guihandles.searchstring.String;
[description,tooltips,interpreter]=parsehelpfile(obj.alltxt.text{allindh});
txtformat=strrep(description,'\n',newline);
ind=strfind(txtformat,ss);
txtformat=setbold(txtformat,ind,length(ss),interpreter);

fname=strrep(obj.alltxt.filenames{allindh},'_','\_');
txt=[setbold(fname, 1,length(fname),interpreter) newline txtformat];
%find searchstring and set bold.
fn=fieldnames(tooltips);
txttt='';
for k=1:length(fn)   
%     indf=strfind(fn{k},ss);
    fnh=fn{k};
%     indt=strfind(tooltips.(fn{k}) ,ss);
    tth=tooltips.(fn{k});
    txttt=[txttt fnh ': ' tth newline];
end


pos=[0 0 .65 1];
obj.annotation=annotation(obj.guihandles.resultspanel,'textbox',pos,...
             'HorizontalAlignment','left',...
             'BackgroundColor','w','FitBoxToText','off','EdgeColor','w',...
             'String',(txt),'Interpreter',interpreter,'FontSize',12);
          obj.annotation.Position=pos;
          obj.guihandles.resultstt.Max=100;
 obj.guihandles.resultstt.String=txttt;       
end

function out=onlynames(in)
for k=length(in):-1:1
    ind=find(in{k}=='.');
    out{k}=in{k}(ind(end-1)+1:ind(end)-1);
end
    
end

function txto=makehelptxt(obj)
settingsdir=obj.getPar('SettingsDirectory');
helpfiledir=[fileparts(settingsdir) filesep 'Documentation' filesep 'help' filesep ];
if strcmp(helpfiledir(1),filesep)
    helpfiledir(1)=[];
end
allfiles=dir([helpfiledir '*.txt']);
for k=length(allfiles):-1:1
    txt{k}=fileread([helpfiledir allfiles(k).name]);
    filenames{k}=allfiles(k).name;
end
txto.text=txt;
txto.filenames=filenames;
obj.helpdir=helpfiledir;
end

function txtformat=setbold(in,ind,lss,interpreter)
if isempty(ind)
    txtformat=in;
    return
end
switch interpreter
    case 'tex'
        openbf='\bf';
        closingbf='\rm';
    case 'latex'
        openbf='\textbf{';
        closingbf='}';
    otherwise
        openbf='';
        closingbf='';
end
for k=length(ind):-1:1
    txtformat=insertAfter(in, ind(k)-1+lss,closingbf);
    txtformat=insertAfter(txtformat, ind(k)-1,openbf);
end
end