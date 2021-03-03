classdef Loader_decode<interfaces.DialogProcessor
    properties
        
    end
    methods
        function obj=Loader_decode(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            out=[];
            [~,~,ext]=fileparts(file);
            switch ext
                case '.h5'
                    [locs,info]=loadh5(file);
                case '.csv'
                    [locs,info]=loadcsv(file);
            end
            locData=interfaces.LocalizationData;
            if strcmp(info.unit,'px')
                pix2nm=info.px_size;
            else
                pix2nm=[1 1];
            end
            filenumber=obj.locData.files.filenumberEnd+1;
            zd=zeros(size(locs.x),'single');
            
            locData.setloc('ynm',single(locs.x*pix2nm(1)));
            locData.setloc('xnm',single(locs.y*pix2nm(2)));
            locData.setloc('znm',single(locs.z));
            locData.setloc('phot',single(locs.phot));
            locData.setloc('frame',double(locs.frame_ix+1));
            locData.setloc('prob',single(locs.prob));
            locData.setloc('LLrel',single(locs.prob));
            
            locData.setloc('filenumber',zd+filenumber);
            locData.setloc('channel',zd);
            
            if ~info.thin
                locData.setloc('ynmerr',single(locs.x_sig*pix2nm(1)));
                locData.setloc('xnmerr',single(locs.y_sig*pix2nm(2)));
                locData.setloc('locprecznm',single(locs.z_sig));
                locData.setloc('bg',single(locs.bg));
                locData.setloc('phot_err',single(locs.phot_sig));
                locData.setloc('locprecnm',(locData.loc.xnmerr+locData.loc.ynmerr)/2);
            else
                locData.setloc('bg',zd);
                locData.setloc('locprecznm',mean(pix2nm)./sqrt(locs.phot)*3);
                locData.setloc('locprecnm',mean(pix2nm)./sqrt(locs.phot));
            end

            obj.locData.addLocData(locData);

            filestruc=locData.files.file;
            filestruc.name=file;
            mx=ceil(max(locData.loc.xnm)/pix2nm(2));
            my=ceil(max(locData.loc.ynm)/pix2nm(1));

            filestruc.info=struct('Width',mx,'Height',my,'roi',[0 0 mx my],'cam_pixelsize_um',pix2nm([2 1])/1000);
            filestruc.info=copyfields(filestruc.info,info);
            if obj.locData.files.filenumberEnd==0
                obj.locData.files.file=filestruc;

            else
                obj.locData.files.file(filenumber)=copyfields(obj.locData.files.file(1),filestruc,fieldnames(obj.locData.files.file(1)));
            end
            obj.locData.files.filenumberEnd=length(obj.locData.files.file);
            obj.locData.files.file(filenumber).number=filenumber;

        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function out=run(obj,p)
            [f,path]=uigetfile(obj.info.extensions);
            if exist([path f],'file')
                obj.load(p,[path f]);
                initGuiAfterLoad(obj);
                out.file=[f,path];
            else
                out.error='file not found. Cannot be loaded.';
            end
        end
        function clear(obj,file,isadd)
            if isadd 
                obj.locData.clear('filter');
            else
                obj.locData.clear;
            end
        end
    end
end


function [locs,io]=loadh5(file)
info=h5info(file);
io.version=info.Groups(2).Attributes.Value;
io.unit=info.Groups(3).Attributes(1).Value;
io.px_size=info.Groups(3).Attributes(2).Value;
for k=1:length(info.Groups(1).Datasets) 
    locs.(info.Groups(1).Datasets(k).Name)=h5read(file,['/data/' info.Groups(1).Datasets(k).Name]);
end
locs.x(:,1)=locs.xyz(1,:);
locs.y(:,1)=locs.xyz(2,:);
locs.z(:,1)=locs.xyz(3,:);

io.thin=true;
if ~isempty(locs.xyz_sig)
    io.thin=false;
    locs.x_sig(:,1)=locs.xyz_sig(1,:);
    locs.y_sig(:,1)=locs.xyz_sig(2,:);
    locs.z_sig(:,1)=locs.xyz_sig(3,:);
end

if ~isempty(locs.xyz_cr)
    locs.x_cr(:,1)=locs.xyz_cr(1,:);
    locs.y_cr(:,1)=locs.xyz_cr(2,:);
    locs.z_cr(:,1)=locs.xyz_cr(3,:);
end

end

function [locs,io]=loadcsv(file)
locs=readtable(file,'NumHeaderLines',3);
fid=fopen(file);
fgetl(fid);
l2=fgetl(fid);
l3=fgetl(fid);
fclose(fid);
io.version=sscanf(l2,'# {%*s "%s}');
io.version(end-2:end)='';
io.unit=sscanf(l3,'# {%*s "%s}');
io.unit(end-2:end)=[];
io.px_size=sscanf(l3,'# {%*s %*s %*s [%f, %f }');
locs=table2struct(locs);
io.thin=false;
if all(isnan(locs.x_sig))
    io.thin=true;
end

end

function loadfile(obj,p,file,mode)
% loades localzation data from a variety of files including text (.csv,
% .txt), hdf5 or MATLAB files. Localization data properties can be
% converted to those used in SMAP, and conversions can be saved for
% repeated use.
[~,~,ext]=fileparts(file);
switch ext
    case {'.csv','.txt'}
        tab=readtable(file);
    case '.mat'
        tab=load(file);
    case {'.hdf5','.h5'}
        info=h5info(file);
        if length(info.Datasets)>1
            answ=listdlg('ListString',{info.Datasets(:).Name});
        else
            answ=1;
        end
        tab=h5read(file,['/' info.Datasets(answ).Name]);
        file=[file '/' info.Datasets(answ).Name];
        if isnumeric(tab)
             if size(tab,1)<size(tab,2)
                 tab=tab';
             end
        end
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
if isstruct(tab)

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
end

if isnumeric(tab)
     
%      if size(tab,1)<size(tab,2)
%          tab=tab';
%      end
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
mx=ceil(max(locData.loc.xnm)/pixnm(1));
my=ceil(max(locData.loc.ynm)/pixnm(end));

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
info.name='Import DECODE .csv/.h5';
info.extensions={'*.csv;*.h5','*.*'};
info.dialogtitle='select a DECODE .csv .mat or .h5 file';
pard.plugininfo=info;  
pard.plugininfo.type='LoaderPlugin';
pard.plugininfo.description='DECODE Loader';

% pard.importdef.object=struct('Style','popupmenu','String',{{'select import'}});
% pard.importdef.position=[1,1];
% pard.importdef.Width=2;
% pard.importdef.TooltipString='Select definition file for import';
end