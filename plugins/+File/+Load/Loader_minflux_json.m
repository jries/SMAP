classdef Loader_minflux_json<interfaces.DialogProcessor
    properties
        
    end
    methods
        function obj=Loader_minflux_json(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file)
            loadfile(obj,p,file);
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function run(obj,p)
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




function loadfile(obj,p,file)
[~,~,ext]=fileparts(file);
switch ext
    case '.mat'
        jt=load(file);
        loc=minfluxmat2loc(jt,p.onlyvalid,~p.simple);
    case '.json'
        txt=fileread(file);
        jt=jsondecode(txt);
        loc=minfluxjson2loc(jt);
end





filenumber=obj.locData.files.filenumberEnd+1;

locData=interfaces.LocalizationData;

locData.loc=loc;
zd=0*loc.xnm;
phot=zd+1000;
if ~isfield(locData.loc,'phot')
    locData.setloc('phot',phot);
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
pixnm=100;


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


function loc=minfluxjson2loc(jt)
valid=[jt.vld];
jtv=jt(valid);
for k=length(jtv):-1:1
    locs=single(jtv(k).itr(end).loc*1e9);
    loc.xnm(k,1)=locs(1);
    loc.ynm(k,1)=locs(2);
    if length(locs)>2
        loc.znm(k,1)=locs(3);
    end
    loc.time(k,1)=single(jtv(k).tim)*1e3;  %from seconds to milliseconds
    loc.iter(k,1)=single(jtv(k).itr(end).itr);
    loc.cfr(k,1)=single(jtv(k).itr(end).cfr);
    loc.dcr(k,1)=single(jtv(k).itr(end).dcr);
    loc.frame(k,1)=k;
end

loc.xnm=loc.xnm-min(loc.xnm);
loc.ynm=loc.ynm-min(loc.ynm);
end

function loc=minfluxmat2loc(jt, onlyvalid,loadall)
numiter=size(jt.cfr,2);
locs=single(jt.loc*1e9);
if onlyvalid
    valid=find(jt.vld);
    indx=sub2ind(size(locs),valid,numiter*ones(size(valid)),ones(size(valid)));
    indy=sub2ind(size(locs),valid,numiter*ones(size(valid)),2*ones(size(valid)));
    indz=sub2ind(size(locs),valid,numiter*ones(size(valid)),3*ones(size(valid)));
    ind2=sub2ind(size(jt.cfr),valid,numiter*ones(size(valid)));
    ind1=jt.vld; 
else  
    lx=jt.loc(:,:,1);
    lxg=~isnan(lx);
    cols=1:numiter;
    lxgc=cols.*lxg;
    lm=max(lxgc,[],2);
    goodind=lm>0;
    g1=find(goodind);g2=lm(goodind);
    indx=sub2ind(size(locs),g1,g2,ones(size(g1)));
    indy=sub2ind(size(locs),g1,g2,2*ones(size(g1)));
    indz=sub2ind(size(locs),g1,g2,3*ones(size(g1)));
    ind2=sub2ind(size(jt.cfr),g1,g2);
    ind1=goodind;
    loc.vld(:,1)=jt.vld(ind1);
    loc.iterations(:,1)=g2;
end

loc.xnm(:,1)=locs(indx);
loc.ynm(:,1)=locs(indy);
znm=locs(indz);
if any(znm>0)
    loc.znm(:,1)=locs(indz);
end
loc.time(:,1)=single(jt.tim(ind1))*1e3;  %from seconds to milliseconds
loc.frame(:,1)=1:length(loc.xnm);
loc.dcr(:,1)=single(jt.eco(ind2));
loc.cfr(:,1)=single(jt.cfr(ind2));

if loadall
    loc.eco(:,1)=single(jt.eco(ind2));
    loc.ecc(:,1)=single(jt.ecc(ind2));
    loc.efo(:,1)=single(jt.efo(ind2));
    loc.efc(:,1)=single(jt.efc(ind2));
    loc.sta(:,1)=single(jt.sta(ind2));
    if isfield(jt,'fbg')
        loc.fbg(:,1)=single(jt.fbg(ind2));
    end
    loc.tid(:,1)=single(jt.tid(ind1));
    loc.act(:,1)=single(jt.act(ind1));
    loc.sky(:,1)=single(jt.sky(ind1));
end

loc.xnm=loc.xnm-min(loc.xnm(loc.vld));
loc.ynm=loc.ynm-min(loc.ynm(loc.vld));
end

function pard=guidef(obj)
info.name='Import MINFLUX Json';
info.extensions={'*.json;*.mat'};
info.dialogtitle='select any .json or .mat  file';
pard.plugininfo=info;  
pard.plugininfo.type='LoaderPlugin';
pard.plugininfo.description='loades localzation data from a variety of files including text (.csv, .txt), hdf5 or MATLAB files. Localization data properties can be converted to those used in SMAP, and conversions can be saved for repeated use.';

pard.simple.object=struct('Style','checkbox','String','load only main fields','Value',1);
pard.simple.position=[1,1];
pard.simple.Width=2;
pard.simple.TooltipString='Only load main localization attributes needed for rendering.';

pard.onlyvalid.object=struct('Style','checkbox','String','load only valid','Value',1);
pard.onlyvalid.position=[2,1];
pard.onlyvalid.Width=2;
pard.onlyvalid.TooltipString='Load only localizations with the value tag == true.';

end