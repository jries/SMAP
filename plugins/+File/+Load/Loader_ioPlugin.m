classdef Loader_ioPlugin<interfaces.DialogProcessor
    methods
        function obj=Loader_ioPlugin(varargin)        
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
            [f,p]=uigetfile(obj.info.extensions);
            obj.load(p,[p f]);
            initGuiAfterLoad(obj);
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




function pard=guidef
info.name='imagin optics loader';
info.extensions={'*.csv';'*.*'};
info.dialogtitle='select any CSV file';
pard.plugininfo=info;     
pard.plugininfo.type='LoaderPlugin';
end

function loadfile(obj,p,file,mode)
%look for file description
path=fileparts(file);
descfile=dir([path filesep 'file-description.xml']);
if ~isempty(descfile)
    pfile=pxml2p(xml2struct([path filesep descfile(1).name]));
else
    pfile.frame=8;pfile.xnano=1;pfile.ynano=2;pfile.znano=3;pfile.intensity=7;
    pfile.sx=4;pfile.sy=5;pfile.sz=6;
end
delimiter='\t';

dat=dlmread(file,delimiter,1,0);
sdat=size(dat);
% filedat=load(file);
% filedat.filename=file;
filenumber=obj.locData.files.filenumberEnd+1;

locData=interfaces.LocalizationData;
locData.addloc('frame',dat(:,pfile.frame));
locData.addloc('xnm',single(dat(:,pfile.xnano)));
locData.addloc('ynm',single(dat(:,pfile.ynano)));
if isfield(pfile,'znano')
locData.addloc('znm',single(dat(:,pfile.znano)));
end
locData.addloc('phot',single(dat(:,pfile.intensity)));

phot=single(dat(:,pfile.intensity));
zd=zeros(sdat(1),1,'single');
locData.addloc('bg',zd);

if isfield(pfile,'sx')
    locData.addloc('locprecnm',single(dat(:,pfile.sx)));
end
if isfield(pfile,'sz')
    locData.addloc('locprecznm',single(dat(:,pfile.sz)));
end
 psfnm=150;psfznm=500;
% locData.addloc('locprecznm',psfznm./sqrt(phot));
% locData.addloc('locprecnm',psfnm./sqrt(phot));
% locData.addloc('PSFxnm',zd+psfnm);
locData.addloc('channel',zd);
locData.addloc('filenumber',zeros(sdat(1),1,'uint8')+filenumber);
locData.removelocs(locData.loc.phot==0);

obj.locData.addLocData(locData);

filestruc=locData.files.file;
filestruc.name=file;
filestruc.info=struct('Width',256,'Height',256,'roi',[0 0 256 256],'pixsize',psfnm/1000);
if obj.locData.files.filenumberEnd==0
    obj.locData.files.file=filestruc;
else
    obj.locData.files.file(filenumber)=copyfields(obj.locData.files.file(1),filestruc,fieldnames(obj.locData.files.file(1)));
end
obj.locData.files.filenumberEnd=length(obj.locData.files.file);

end

function p=pxml2p(pxml)
fna=fieldnames(pxml.description);
fn=setdiff(fna,{'separator'});
pn=copyfields([],pxml.description,fn);
for k=1:length(fn)
    p.(fn{k})=str2double(pn.(fn{k}).Text)+1;
end
end
