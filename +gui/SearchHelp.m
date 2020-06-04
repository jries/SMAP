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
    obj.figure=figure('MenuBar','none','Name','Search in Help','Units','pixels');
    smappos=obj.getPar('mainGuihandle').Position;
    obj.figure.Position(1)=smappos(1)+smappos(3);
    obj.figure.Position(2)=smappos(2);
    obj.figure.Position(3:4)=[700 800];
    obj.figure.Units='normalized';
    obj.handle=obj.figure;
end
hf=obj.figure;
h.searchstring=uicontrol(hf, 'Style','edit','Units','normalized','Position',[0.02,0.93,0.75,0.05]);
h.searchbutton=uicontrol(hf, 'Style','pushbutton','String','Search','Units','normalized','Position',[0.8,0.93,0.18,0.05],...
    'Callback',{@search_callback,obj});
h.resultslist=uicontrol(hf,'Style','listbox','Units','normalized','Position',[0.02,0.02,0.3,0.9],'Callback',{@list_callback,obj});
h.resultspanel=uipanel(hf,'Units','normalized','Position',[0.32,0.22,0.66,0.7],'BackgroundColor','w');
h.resultstt=uicontrol(hf,'Style','edit','Units','normalized','Position',[0.32,0.02,0.66,0.2],'BackgroundColor','w','HorizontalAlignment','left','Max',100);

obj.guihandles=h;
end

function search_callback(a,b,obj)
alltxt=obj.getPar('allhelpfiles');
if isempty(alltxt)
    alltxt=makehelptxt(obj);
    obj.setPar('allhelpfiles',alltxt);
end
obj.alltxt=alltxt;
searchstring=obj.guihandles.searchstring.String;
ind=find(contains(alltxt.text,searchstring,'IgnoreCase',true));
if isempty(ind)
    obj.guihandles.resultslist.String={'no results'};
    delete(obj.annotation)
    obj.guihandles.resultstt.String='';
    return
end
obj.searchind=ind;
obj.guihandles.resultslist.String=onlynames(alltxt.filenames(ind));
% obj.functionlist=alltxt.filenames(ind);
% obj.functiontxt=alltxt.text(ind);
list_callback(obj.guihandles.resultslist,0,obj);
end

function list_callback(a,b,obj)
select=a.Value;
allindh=obj.searchind(select);
for  k=1:length(obj.annotation)
delete(obj.annotation{k})
end
ss=obj.guihandles.searchstring.String;
[description,tooltips]=parsehelpfile(obj.alltxt.text{allindh});
fn=fieldnames(description);
for k=1:length(fn)
    ind=strfind(description.(fn{k}),ss);
    description.(fn{k})=setbold(description.(fn{k}),ind,length(ss),fn{k});
end
% txtformat=strrep(description,'\n',newline);
% ind=strfind(txtformat,ss);
% txtformat=setbold(txtformat,ind,length(ss),interpreter);
if isempty(tooltips)
    txttt{1}='';
else
    fname=strrep(obj.alltxt.filenames{allindh},'_','\_');
    % txt=[setbold(fname, 1,length(fname),interpreter) newline txtformat];
    %find searchstring and set bold.
    fn=fieldnames(tooltips);
    txttt='';
    for k=1:length(fn)   
    %     indf=strfind(fn{k},ss);
        fnh=fn{k};
    %     indt=strfind(tooltips.(fn{k}) ,ss);
        tth=strrep(tooltips.(fn{k}), newline,' ');
        tth=regexprep(tth,' +',' ');
        txttt{k}=[fnh ': ' tth];
    %     txttt=[txttt fnh ': ' tth newline];
    end
end
pos=[0 0 .65 1];
obj.annotation=showpluginhelp(obj.guihandles.resultspanel,description,pos);
% obj.annotation=annotation(obj.guihandles.resultspanel,'textbox',pos,...
%              'HorizontalAlignment','left',...
%              'BackgroundColor','w','FitBoxToText','off','EdgeColor','w',...
%              'String',(txt),'Interpreter',interpreter,'FontSize',12);
%           obj.annotation.Position=pos;
obj.guihandles.resultstt.Max=100;
 obj.guihandles.resultstt.String=txttt;       
end

function out=onlynames(in)
out={};
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
        openbf='\color{red}';
        closingbf='\color{black}';
%         openbf='\bf';
%         closingbf='\rm';
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