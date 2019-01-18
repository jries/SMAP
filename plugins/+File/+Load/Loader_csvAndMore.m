classdef Loader_csvAndMore<interfaces.DialogProcessor
    properties
        loaderpath='settings/csvloaderconversion/';
    end
    methods
        function obj=Loader_csvAndMore(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            if nargin<4
                mode=getfilemode(file);
            end
            loadfile(obj,p,file,mode);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function run(obj,p)
            [f,path]=uigetfile(obj.info.extensions);
            obj.load(p,[path f]);
            initGuiAfterLoad(obj);
        end
        function clear(obj,file,isadd)
            if isadd 
                obj.locData.clear('filter');
            else
                obj.locData.clear;
            end
        end
        function initGui(obj)
            obj.loaderpath='settings/csvloaderconversion/';
            files=dir([obj.loaderpath '*.txt']);
            string={'New format', files(:).name};
            obj.guihandles.importdef.String=string;
%             obj.guihandles.importdef.Callback=@obj.csvconvert_callback;
        end
%         function csvconvert_callback(obj,object,b)
%            if object.Value==1
%                newformat(obj);
%            end
%         end
    end
end

% function newformat(obj)
% [f,p]=uigetfile(obj.info.extensions);
% end

function table_callback(object,value)
if value.Indices(2)==4 %popupmenu
    if strcmp(value.NewData,'edit')
        object.Data(value.Indices(1),3)=object.Data(value.Indices(1),1);
    else
    object.Data{value.Indices(1),3}=value.NewData;
    end

end
end

function pfileo=importdialog(tab,pfile,obj)
 defnames={'edit','ID','frame','channel','xnm','ynm','znm','phot','bg','locprecnm','locprecznm','PSFxnm','PSFynm'};
%         varnames=tab.Properties.VariableNames;
        varnames=fieldnames(tab);
        columnformat={'char','logical','char',defnames};
%         f=dialog;
f=figure;
        f.Position(3:4)=[400,500];
        ht=uitable('units','normalized','Parent',f,'Position',[0 0.2 1 .5]);
        ht.ColumnFormat=columnformat;
        ht.ColumnEditable=[false true true true];
        ht.ColumnName={'Column name','import','field name','select field'};
        ds(:,1)=varnames;ds(:,3)=varnames;
        ds(:,4)=repmat({'edit'},length(varnames),1);
        
        ht2=uitable('units','normalized','Parent',f,'Position',[0 0.7 1 .3]);
        ht2.ColumnName=varnames;
        a=table2array(struct2table(tab));
%         a=struct2cell(tab);
        ht2.Data=a(1:100,:);
        
        fn=fieldnames(pfile);
        for k=1:length(fn)
            vh=pfile.(fn{k});
            if isnumeric(vh)
                ds{vh,3}=fn{k};
                try
                    ds{vh,4}=fn{k};
                    ds{vh,2}=true;
                catch
                end
            end
        end
        
        ht.Data=ds;
        ht.CellEditCallback=@table_callback;
        
        uicontrol('Parent',f,'Style','text','String','Cam pixel size (nm)','Position',[5,50,150,20]);
        hcam=uicontrol('Parent',f,'Style','edit','String','100','Position',[155,50,100,20]);
        
        uicontrol('Parent',f,'Style','text','String','Factor (xy or [xy z])','Position',[5,70,150,20]);
        hfac=uicontrol('Parent',f,'Style','edit','String','1','Position',[155,70,100,20]);
        
        b1=uicontrol('Parent',f,'Style','pushbutton','String','Cancel','Callback',{@button_callback,0},'Position',[5,5,50,20]);
        b2=uicontrol('Parent',f,'Style','pushbutton','String','Save conversion structure','Callback',{@button_callback,2},'Position',[70,5,190,20]);
        b3=uicontrol('Parent',f,'Style','pushbutton','String','Ok','Callback',{@button_callback,1},'Position',[275,5,50,20]);
        
         uiwait(f);
        
    function button_callback(a,b,number)
       dsh=ht.Data; 
       
       if number>0
           sd=size(dsh);
       for k2=1:sd(1)
             if dsh{k2,2}
                 pfileo.(dsh{k2,3})=dsh{k2,1};
             end
       end
       else
           pfileo=[];
       end
       pfileo.cam_pixelsize_um=str2double(hcam.String);
       pfileo.factor=str2num(hfac.String);
       if number==2 %save
%            path='settings/csvloaderconversion/';
           [ file, path]=uiputfile([obj.loaderpath 'csv_.txt']);
           writestruct([path file],pfileo);
           disp('converstion structure saved');
           obj.initGui;
       end
       delete(f);
    end
% pfileo=pfile;
        
end

function loadfile(obj,p,file,mode)

[~,~,ext]=fileparts(file);
switch ext
    case {'.csv','.txt'}
        tab=readtable(file);
    case '.mat'
        tab=load(file);
    case '.hdf5'
        info=h5info(file);
        tab=h5read(file,['/' info.Datasets(1).Name]);

end

if isstruct(tab)
fn=fieldnames(tab);
while length(fn)==1 
    tab=tab.(fn{1});
    if isstruct(tab)
        fn=fieldnames(tab);
    else
        break
    end
end
% remove other entries
fn=fieldnames(tab);
tabo=[];
for k=1:length(fn)
    vh=tab.(fn{k});
    if isnumeric(vh) && length(vh)>10 && sum(size(vh)>1)==1
        tabo.(fn{k})=tab.(fn{k});
    end
end
tab=tabo;
end

if isnumeric(tab)
     tab=array2table(tab);
%     tab=table2struct(tab);
% elseif isstruct(tab)
%    tab=struct2table(tab);
end

if istable(tab)
    tab=table2struct(tab,'ToScalar',true);
    fns=fieldnames(tab);
    for k=1:length(fns)
    if iscell(tab.(fns{k}))
        tab=rmfield(tab,fns{k});
    end
    end
end

if ~isfield(p,'importdef') %if called from loader
    p.importdef.Value=1;
end
if p.importdef.Value==1 %new file
%look for file description
path=fileparts(file);
descfile=dir([path filesep 'file-description.xml']);
    if ~isempty(descfile)
        pfile=pxml2p(xml2struct([path filesep descfile(1).name]));
    else
        
        pfile.frame=1;pfile.xnm=2;pfile.ynm=3;pfile.znm=4;pfile.phot=5;
        fnh=lower(fieldnames(tab));
        ind=find(~cellfun(@isempty,strfind(fnh,'frame'))); if ~isempty(ind), pfile.frame=ind(1); end
        ind=find(~cellfun(@isempty,strfind(fnh,'x'))); if ~isempty(ind), pfile.xnm=ind(1); end
        ind=find(~cellfun(@isempty,strfind(fnh,'y'))); if ~isempty(ind), pfile.ynm=ind(1); end
        ind=find(~cellfun(@isempty,strfind(fnh,'z'))); if ~isempty(ind), pfile.znm=ind(1); end
        ind=find(~cellfun(@isempty,strfind(fnh,'intens'))); if ~isempty(ind), pfile.phot=ind(1); end
        ind=find(~cellfun(@isempty,strfind(fnh,'phot'))); if ~isempty(ind), pfile.phot=ind(1); end
        ind=find(~cellfun(@isempty,strfind(fnh,'bg'))); if ~isempty(ind), pfile.bg=ind(1); end
        ind=find(~cellfun(@isempty,strfind(fnh,'backg'))); if ~isempty(ind), pfile.bg=ind(1); end

    end
        %dialog: rename
     pfile=importdialog(tab,pfile,obj);     
     if isempty(pfile)
         return
     end
else
    cfile=[obj.loaderpath p.importdef.selection];
    pfile=readstruct(cfile);
    %load pfile structure
end



% header=csvread(file,[0 pfile.firstrow],0);
% dat=csvread(file,pfile.firstrow-1,0);
% sdat=size(dat);
% filedat=load(file);
% filedat.filename=file;
% try
filenumber=obj.locData.files.filenumberEnd+1;

locData=interfaces.LocalizationData;
fn=setdiff(fieldnames(pfile),{'cam_pixelsize_um','factor'});


for k=1:length(fn)
    fnsmap=fn{k};
    fnimport=pfile.(fnsmap);
    if strcmp(fnsmap,'frame')
        locData.setloc(fnsmap,double(tab.(fnimport)));
    elseif isfield(tab,fnimport)
        locData.setloc(fnsmap,single(tab.(fnimport)));
    end
end

fieldsxy=intersect(fieldnames(locData.loc),{'xnm','ynm','locprecnm','PSFxnm','PSFynm'});
fieldsz=intersect(fieldnames(locData.loc),{'znm','locprecznm'});
if isfield(pfile,'factor')&&any(pfile.factor~=1)
    facxy=pfile.factor(1);
    if length(pfile.factor)>1
        facz=pfile.factor(end);
    else
        facz=1;
    end
    for k=1:length(fieldsxy)
        locData.loc.(fieldsxy{k})=locData.loc.(fieldsxy{k})*facxy;
    end
    for k=1:length(fieldsz)
        locData.loc.(fieldsz{k})=locData.loc.(fieldsz{k})*facz;
    end
end

if ~isfield(locData.loc,'xnm')||~isfield(locData.loc,'ynm')
    warning('Define at least the fields xnm and ynm')
end
zd=zeros(size(tab.(pfile.(fn{1}))),'single');
if ~isfield(locData.loc,'frame')
    locData.setloc('frame',double(zd+1));
end
if ~isfield(locData.loc,'phot')
    locData.setloc('phot',0*single(zd+1)+1000);
end
% locData.addloc('frame',dat(:,pfile.frame));
% locData.addloc('xnm',single(dat(:,pfile.xnm)));
% locData.addloc('ynm',single(dat(:,pfile.ynm)));
% if isfield(pfile,'znano')
% locData.addloc('znm',single(dat(:,pfile.znm)));
% end
% locData.addloc('phot',single(dat(:,pfile.phot)));


if isfield(locData.loc,'phot')
    phot=locData.loc.phot;
else
    phot=zd+1000; %default: 1000 photons
end
psfnm=150;psfznm=500;

if ~isfield(locData.loc,'locprecznm')&&isfield(locData.loc,'znm')
    locData.setloc('locprecznm',psfznm./sqrt(phot));
end

if ~isfield(locData.loc,'locprecnm')
    locData.setloc('locprecnm',psfnm./sqrt(phot));
end

if ~isfield(locData.loc,'bg')
    locData.setloc('bg',zd);
end

if ~isfield(locData.loc,'PSFxnm')
    locData.setloc('PSFxnm',zd+psfnm);
end

if ~isfield(locData.loc,'channel')
    locData.setloc('channel',zd);
end

locData.setloc('filenumber',zd+filenumber);
% catch err
%     err
%     warndlg('import did not work. Select a different import format file or define new file')
%     return
% %     rethrow(err);
% end



if isfield(pfile,'cam_pixelsize_um')
    pixnm=pfile.cam_pixelsize_um;
else
    pixnm=100;
end

obj.locData.addLocData(locData);

filestruc=locData.files.file;
filestruc.name=file;
mx=ceil(max(locData.loc.xnm)/pixnm);
my=ceil(max(locData.loc.ynm)/pixnm);

filestruc.info=struct('Width',mx,'Height',my,'roi',[0 0 mx my],'cam_pixelsize_um',pixnm/1000);
if obj.locData.files.filenumberEnd==0
    obj.locData.files.file=filestruc;
    
else
    obj.locData.files.file(filenumber)=copyfields(obj.locData.files.file(1),filestruc,fieldnames(obj.locData.files.file(1)));
end
obj.locData.files.filenumberEnd=length(obj.locData.files.file);
obj.locData.files.file(filenumber).number=filenumber;



end

function [po,p]=pxml2p(pxml)
fna=fieldnames(pxml.description);
fn=setdiff(fna,{'separator'});
pn=copyfields([],pxml.description,fn);
for k=1:length(fn)
    p.(fn{k})=str2double(pn.(fn{k}).Text)+1;
end
replace={'xnano','xnm','ynano','ynm','znano','znm','intensity','phot','frame','frame'};
for k=1:2:length(replace)
    if isfield(p,replace{k})
    po.(replace{k+1})=p.(replace{k});
    end
end
if isfield(p,'firstrow')
    p.firstrow=p.firstrow+1;
else
    p.firstrow=1;
end
end


function pard=guidef
info.name='Import CSV/MAT/HDF5';
info.extensions={'*.csv;*.xls;*.mat;*.hdf5;*.txt','*.*'};
info.dialogtitle='select any .csv .mat or .hdf5 file';
pard.plugininfo=info;  
pard.plugininfo.type='LoaderPlugin';

pard.importdef.object=struct('Style','popupmenu','String',{{'select import'}});
pard.importdef.position=[1,1];
pard.importdef.Width=2;
pard.importdef.TooltipString='Select definition file for import';
end