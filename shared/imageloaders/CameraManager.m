classdef CameraManager<interfaces.GuiParameterInterface
% Stores all cameras with metadata structure and default setttings in a
% data base. This eliminates the need to specify camera or camera settings
% when loading raw image data.
    
    properties
        handle
        imloader
        cameras
        guihandles
        currentcam=1;
        currentstate=1;
        lastcamtableselected=[];
        defaultpath;
        cameraSettingsFile='settings/cameras.mat';
    end
    
    methods
        function obj=CameraManager
%             makeGui(obj)
        end
        function loadimages(obj,file) 
            if nargin > 1
                obj.imloader=imageloaderAll(file,[],obj.cameraSettingsFile);
            end
                [par,cam,state]=getCameraCalibration(obj.imloader,[],true,obj.cameraSettingsFile);
                if isempty(cam)
                    answ=questdlg('camera not found. Create new camera?');
                    if strcmp(answ,'Yes')
                        createnewcamera(obj);
                    end
                else
                    obj.currentcam=cam;
                    
                    prop2table(obj);
                end
                
                if isempty(state)&&~isempty(cam)
                    answ=questdlg('State not recognized. Create new state?');
                    if strcmp(answ,'Yes')
                        stateadd(0,0,obj,'add')
                    end
                else
                    obj.currentstate=state;
                    prop2table(obj);
                end
        end
        function makeGui(obj)
            makeGui(obj)
        end
    end
    
end

function makeGui(obj)
width=800;
height=700;
lineheight=25;
posbutton=width*.8;
buttonwidth=width*.15;
if isempty(obj.handle)||~isvalid(obj.handle)
     obj.handle=figure('Units','normalized','Units','pixels','Position',[150,200,width,height],'Name','CameraManager','NumberTitle','off');
     obj.handle.ToolBar='none';
     obj.handle.MenuBar='none';
     delete(obj.handle.Children);
end
%load images
hp=uicontrol('Style','text','String','camera Settings File','Position',[20 height-40,width*.2,lineheight]);
obj.guihandles.camerafile=uicontrol('Style','edit','String',obj.cameraSettingsFile,'Position',[width*.2 height-40,width*.6,lineheight]);
hp=uicontrol('Style','pushbutton','String','load','Position',[posbutton height-40,buttonwidth,lineheight],'Callback',{@loadcamerafile,obj});


hp=uicontrol('Style','pushbutton','String','Load images','Position',[posbutton height-90,buttonwidth,lineheight],'Callback',{@loadimages,obj});
hp=uicontrol('Style','pushbutton','String','test','Position',[posbutton height-115,buttonwidth,lineheight],'Callback',{@testcal,obj});
hp=uicontrol('Style','pushbutton','String','Add camera','Position',[posbutton height-140,buttonwidth,lineheight],'Callback',{@menu_callback,obj,'add'});


tcam=uitable(obj.handle,'Position',[10 height-lineheight*8-50 width-200,lineheight*7.5]);
tcam.ColumnName={'Camera Name','ID field','ID'};
tcam.Data={'Default','select Cam_ID','001'};
wh=tcam.Position(3);
tcam.ColumnWidth={'auto',wh*.5,wh*.27};
tcam.CellSelectionCallback={@cellselect,obj,'cam'};
tcam.ColumnEditable=[ true false true];

hc=uicontextmenu(obj.handle);
hui=uimenu('Parent',hc,'Label','add','Callback',{@menu_callback,obj});
hui=uimenu('Parent',hc,'Label','remove','Callback',{@menu_callback,obj});
hui=uimenu('Parent',hc,'Label','move up','Callback',{@menu_callback,obj});
hui=uimenu('Parent',hc,'Label','move down','Callback',{@menu_callback,obj});

tcam.UIContextMenu=hc;
tpar=uitable(obj.handle,'Position',[10 height-lineheight*17-80 width-40,lineheight*9.5]);
tpar.ColumnName={'Parameter','mode','fixvalue','metafield','Value','conversion','Converted'};

% parnames={'EMon','pixelsize','conversion','emgain','offset','roi','exposure','timediff','comment'};
% dat=cell(length(parnames),7);
% dat(:,1)=parnames;
% dat{1,2}='fix';
dat=intpartable;


wh=tpar.Position(3);
tpar.ColumnWidth={wh*.12,'auto','auto',wh*.25,'auto',wh*.18,'auto'};
tpar.CellSelectionCallback={@cellselect,obj,'par'};
tpar.CellEditCallback={@celledit,obj,'par'};
tpar.ColumnFormat={'char',{'fix','metadata','state dependent'},'char','char','char','char','char'};
tpar.ColumnEditable=[false true true false false true false];
tpar.Data=dat;


tstates=uicontrol('Style','listbox','String',{'State 1'},'Position',[10 50 width*.15,lineheight*4],'Callback',{@statecallback,obj});
tstatesadd=uicontrol('Style','pushbutton','String','add','Position',[10 lineheight*3+80 width*.06,lineheight],'Callback',{@stateadd,obj,'add'});
tstatesrem=uicontrol('Style','pushbutton','String','rem','Position',[width*.09+10 lineheight*3+80 width*0.06,lineheight],'Callback',{@stateadd,obj,'rem'});
uicontrol('Style','text','String','State defining parameters','Position',[width*.2 lineheight*3+80 width*.2,lineheight])
tstatesrem=uicontrol('Style','pushbutton','String','close','Position',[posbutton 15 buttonwidth,lineheight],'Callback',{@close_callback,obj});
hp=uicontrol('Style','pushbutton','String','Save','Position',[posbutton-buttonwidth-15, 15,buttonwidth,lineheight],'Callback',{@savecameras,obj});

hp=uicontrol('Style','pushbutton','String','Calibrate Camera','Position',[15, 15,buttonwidth*1.2,lineheight],'Callback',{@calibrate_cameras,obj});


tdef=uitable(obj.handle,'Position',[width*.18 50 width*.5,lineheight*4]);
tdef.ColumnName={'Meta Field','Value'};
tdef.Data={'select','';'select','';'select','';'select',''};
tdef.CellSelectionCallback={@cellselect,obj,'def'};
wh=tdef.Position(3);
tdef.ColumnWidth={wh*.55,wh*.35};
tdef.ColumnEditable=[false true];


tval=uitable(obj.handle,'Position',[width*.7 50 width-40-width*.7,lineheight*5]);
wh=tval.Position(3);
tval.ColumnName={'Parameter','Value'};
tval.ColumnEditable=[false true];
tval.Data=tpar.Data(:,[1 3]);
tval.ColumnWidth={wh*.45,wh*.25};


obj.guihandles.camtable=tcam;
obj.guihandles.partable=tpar;
obj.guihandles.statelist=tstates;
obj.guihandles.statedeftable=tdef;
obj.guihandles.statevaltable=tval;

tables2prop(obj);
showpartable(obj);

% obj.cameras(1).par=tpar.Data;
% obj.cameras(1).ID=struct('name',tcam.Data{1,1},'tag',tcam.Data{1,2},'value',tcam.Data{1,3});
% statestruct=struct('statelist',{tstates.String},'defpar',{tdef.Data},'par',{tval.Data});
% obj.cameras(1).state(1)=statestruct;
loadcameras(obj);
end

function celledit(table,data,obj,tname)
switch tname
    case 'cam'
%         obj.cameras(obj.currentcam).par=obj.guihandles.partable.Data;
%         obj.guihandles.partable.Data=obj.cameras(data.Indices(1)).par;
%         obj.currentcam=data.Indices(1);
    case 'par'
        if data.Indices(2)==2
            showpartable(obj);
        end
    case 'def'
%         indtag=1;
%         indval=2;
%         tab=obj.guihandles.statedeftable;
end
end

function showpartable(obj)
partable=obj.guihandles.partable.Data;
tval=obj.guihandles.statevaltable.Data;

dat=(partable(:,[1 2]));

col=ones(size(tval,1),3)*.7;
indg=strcmp(dat(:,2),'state dependent');
col(indg,:)=1;
obj.guihandles.statevaltable.BackgroundColor=col;
% obj.guihandles.statevaltable.Data=dat;
end

function loadcameras(obj)
% file='settings/cameras.mat';
file =obj.cameraSettingsFile;
if ~exist(file,'file')
    disp('camera file does not exist')
    return
end
l=load(file);
obj.cameras=l.cameras;
obj.guihandles.camtable.Data=l.camtab;
prop2table(obj);
showpartable(obj);
% obj.guihandles.partable.Data=l.cameras(1).par;
% 
% obj.guihandles.statedeftable.Data=l.cameras(1).state(1).defpar;
% obj.guihandles.statevaltable.Data=l.cameras(1).state(1).par;
% obj.guihandles.statelist.String=l.cameras(1).state(1).statelist;
end

function savecameras(a,b,obj)
% file='settings/cameras.mat';
file=obj.cameraSettingsFile;
tables2prop(obj);
camtab=obj.guihandles.camtable.Data;
cameras=obj.cameras;
save(file,'cameras','camtab')
end

function tables2prop(obj)
cameras=obj.cameras;
cameras(obj.currentcam).par=obj.guihandles.partable.Data;
cameras(obj.currentcam).state(obj.currentstate)=struct('statelist',{obj.guihandles.statelist.String},...
    'defpar',{obj.guihandles.statedeftable.Data},'par',{obj.guihandles.statevaltable.Data});
camtab=obj.guihandles.camtable.Data;
s=size(camtab);
for k=1:s(1)
    cameras(k).ID=struct('name',camtab{k,1},'tag',camtab{k,2},'value',camtab{k,3});
end
obj.cameras=cameras;
end

function prop2table(obj)

%hack to add property for roi mode
pardat=obj.cameras(obj.currentcam).par;
if size(pardat,1)<13
    pardat(13,:)={'roimode','fix','none','select',[],'',[]};
end
obj.guihandles.partable.Data=pardat;

obj.guihandles.statelist.Value=obj.currentstate;
obj.guihandles.statedeftable.Data=obj.cameras(obj.currentcam).state(obj.currentstate).defpar;
statepar=obj.cameras(obj.currentcam).state(obj.currentstate).par;
%hack to add property for roi mode
if size(statepar,1)<13
    statepar(13,:)={'roimode','none'};
    obj.cameras(obj.currentcam).state(obj.currentstate).par=statepar;
end
obj.guihandles.statevaltable.Data=statepar;

cams=obj.cameras;
for k=1:length(cams)
    data{k,1}=cams(k).ID.name;
    data{k,2}=cams(k).ID.tag;
    data{k,3}=cams(k).ID.value;
end

nstates=length(obj.cameras(obj.currentcam).state);
for k=1:nstates
    statestr{k}=['State ' num2str(k)];
end
obj.guihandles.statelist.String=statestr;

obj.guihandles.camtable.Data=data;
col=ones(size(data,1),3);
col(obj.currentcam,1)=.3;
obj.guihandles.camtable.BackgroundColor=col;
showpartable(obj)
end

function t=intpartable
t=cell(12,7);
parnames={'EMon','cam_pixelsize_um','conversion','emgain','offset','roi','exposure','timediff','comment','numberOfFrames','Width','Height','roimode'};
mode={'fix','fix','fix','fix','fix','fix','fix','fix','fix','metadata','metadata','metadata','fix'};
default={'1','0.1','1','100','100','','1','1','settings not initialized','0','0','0','none'};
conversion={'str2double(X)','str2double(X)','str2double(X)','str2double(X)','str2double(X)','str2num(X)','str2double(X)','str2double(X)','','str2double(X)','str2double(X)','str2double(X)',''};
metafield={'select','select','select','select','select','select','select','select','select','select','select','select','select','select'};

for k=1:size(t,1)
    t{k,1}=parnames{k};
    t{k,2}=mode{k};
    t{k,3}=default{k};
    t{k,4}=metafield{k};
    t{k,6}=conversion{k};
end
end

% function cellselect_cam(table,data,obj,indtag,indval)
% if isempty(data.Indices)
%     return
% end
% if data.Indices(2)==indtag
%     ma=obj.imloader.getmetadatatags;
%      tag = gettag(ma);
%      if ~isempty(tag)
%         table.Data(data.Indices(1),[indtag indval])=tag;
%      end
% end
% obj.cameras(obj.currentcam).par=obj.guihandles.partable.Data;
% obj.guihandles.partable.Data=obj.cameras(data.Indices(1)).par;
% obj.currentcam=data.Indices(1);
% end

function cellselect(table,data,obj,tname)
if isempty(data.Indices)
    return
end
switch tname
    case 'cam'
        indtag=2;
        indval=3;
    case 'par'
        indtag=4;
        indval=5;
    case 'def'
        indtag=1;
        indval=2;
end

if data.Indices(2)==indtag
    if isempty(obj.imloader)
        warndlg('please load images before assigning fields')
        return     
    end
    ma=obj.imloader.getmetadatatags;
     tag = gettag(ma);
     if ~isempty(tag)
        table.Data(data.Indices(1),[indtag indval])=tag;
        
        if strcmp(tname,'par')
            X=tag{2};
            if ~isempty(table.Data{data.Indices(1),6})&&~isempty(X)
                X=eval(table.Data{data.Indices(1),6});
            end
            if isnumeric(X) && length(X)>1
                X=num2str(X);
            end
            table.Data{data.Indices(1),7}=X;
        end
        tables2prop(obj);
     end
end


switch tname
    case 'cam'
        obj.cameras(obj.currentcam).par=obj.guihandles.partable.Data;
%         obj.guihandles.partable.Data=obj.cameras(data.Indices(1)).par;
        obj.currentcam=data.Indices(1);
        obj.lastcamtableselected=data.Indices;
        obj.currentstate=1;
        prop2table(obj);
        
    case 'par'
%         if data.Indices(2)==2
%             showpartable(obj);
%         end
    case 'def'
        indtag=1;
        indval=2;
        tab=obj.guihandles.statedeftable;
end

end

function tag = gettag(ma)
f=figure;
tc=uitable(f);
[~,ind]=sortrows(ma(:,1));


tc.Data=ma(ind,:);
tc.Position(3)=f.Position(3)-30;
tc.Position(2)=100;
tc.CellSelectionCallback=@cellselecth;
w=tc.Position(3);
tc.ColumnWidth={w*.75,.2*w};
uicontrol('Style','pushbutton','String','Ok','Position',[200 10 100 20],'Callback',@buttoncallback)
uicontrol('Style','pushbutton','String','Cancel','Position',[10 10 100 20],'Callback',@buttoncallback)
pos=[];
waitfor(f)

    function buttoncallback(a,b)
        if ~isempty(pos)&&strcmp(a.String,'Ok')
        tag=tc.Data(pos(1),:);
        else
            tag={};
        end
        close(f)
    end
    function cellselecth(t,d)
        pos=d.Indices;
    end
end


function loadimages(a,b,obj)
ph=obj.defaultpath;
if ~isempty(ph)
    ph=[fileparts(ph) filesep];
else
    ph='';
end
[file path]=uigetfile([ph '*.*']);
if file
    obj.loadimages([path file]);
    obj.defaultpath=path;
% obj.imloader=imageloaderAll([path file]);
end
end

function createnewcamera(obj)
l=length(obj.cameras)+1;
dat=obj.guihandles.camtable.Data;
dat(l,:)={'new','select',''};
obj.guihandles.camtable.Data=dat;
obj.cameras(l)=obj.cameras(1);
obj.cameras(l).par=intpartable;
obj.cameras(l).ID=struct('name','new','tag','select','value','');
cellselect(obj.guihandles.camtable,struct('Indices',[l,1]),obj,'cam');
end

function statecallback(object,data,obj)
newstate=object.Value;
tables2prop(obj);
obj.currentstate=newstate;
prop2table(obj);
end

function stateadd(a,b,obj,addrem)
switch addrem
    case 'add'
%         states=obj.guihandles.statelist.String;
        l=length(obj.cameras(obj.currentcam).state);
%         states{l+1}=['State ' num2str(l+1)];
        newstate=obj.cameras(obj.currentcam).state(obj.currentstate);
        obj.currentstate=l;
        try
            for k=1:size(newstate.defpar,1)
                fieldh=newstate.defpar{k,1};
                if ~strcmp(fieldh,'select')
                    v=obj.imloader.gettag(fieldh);
                    newstate.defpar{k,2}=v;
                end
            end
            
        end
        obj.cameras(obj.currentcam).state(l+1)=newstate;
        obj.currentstate=l+1;
        prop2table(obj);
%         obj.guihandles.statelist.String=states;
%         obj.guihandles.statelist.Value=l+1;
    case 'rem'
        obj.cameras(obj.currentcam).state(obj.currentstate)=[];
        obj.currentstate=1;
        prop2table(obj);
end
end

function menu_callback(object, data, obj,label)
if nargin<4
    label=object.Label;
end
switch label
    case 'remove'
        if isempty(obj.lastcamtableselected)
            return
        end
        if length(obj.cameras)==1
            warning('at least one camera required. You can overwrite the values')
            return
        end
        tables2prop(obj);
        obj.cameras(obj.lastcamtableselected(1))=[];
        obj.currentcam=1;
        prop2table(obj);
        return
    case 'add'
        createnewcamera(obj);
        return
    case 'move up'
       newpos=max(obj.currentcam-1,1);
    case 'move down'
        newpos=min(obj.currentcam+1,length(obj.cameras));
end
oldpos=obj.currentcam;
obj.cameras([oldpos newpos])=obj.cameras([newpos oldpos]);
prop2table(obj);
end

function testcal(a,b,obj)
tables2prop(obj);
p=getCameraCalibration(obj.imloader,obj,true,obj.cameraSettingsFile);
p
end

function close_callback(a,b,obj)
close(obj.handle)
end

function calibrate_cameras(a,b,obj)
%display info on how to calibrate cameras.
end

function loadcamerafile(a,b,obj)
file=obj.cameraSettingsFile;
[f,p]=uigetfile(file);
if f
    obj.cameraSettingsFile=[p f];
    obj.guihandles.camerafile=[p f];
    loadcameras(obj);
end


end