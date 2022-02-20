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
        if ~isfield(jt,'loc') && isfield(jt,'itr')
            
            jt.itr.tim=jt.tim;
            jt.itr.vld=jt.vld;
            jt.itr.tid=jt.tid;
            jt.itr.act=jt.act;
            jt.itr.sky=jt.sky;
            jt=jt.itr;
        end
        loc=minfluxmat2loc(jt,p.onlyvalid,~p.simple);
        xoff=min(loc.xnm(loc.vld));yoff=min(loc.ynm(loc.vld));
        loc.xnm=loc.xnm-xoff;
        loc.ynm=loc.ynm-yoff;
        if isfield(loc,'xncnm')
            loc.xncnm=loc.xncnm-xoff;
            loc.yncnm=loc.yncnm-yoff;
        end
    case '.json'
        txt=fileread(file);
        jt=jsondecode(txt);
        loc=minfluxjson2loc(jt);
end
filenumber=obj.locData.files.filenumberEnd+1;

locData=interfaces.LocalizationData;

locData.loc=loc;
zd=0*loc.xnm;
phot=loc.phot;
% phot=zd+1000;
% if ~isfield(locData.loc,'phot')
%     locData.setloc('phot',phot);
% end

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



function pard=guidef(obj)
info.name='Import MINFLUX Abberior';
info.extensions={'*.json;*.mat'};
info.dialogtitle='select any .json or .mat  file';
pard.plugininfo=info;  
pard.plugininfo.type='LoaderPlugin';
pard.plugininfo.description='loades localzation data from a variety of files including text (.csv, .txt), hdf5 or MATLAB files. Localization data properties can be converted to those used in SMAP, and conversions can be saved for repeated use.';

pard.simple.object=struct('Style','checkbox','String','load only main fields','Value',1);
pard.simple.position=[1,1];
pard.simple.Width=2;
pard.simple.TooltipString='Only load main localization attributes needed for rendering.';

pard.onlyvalid.object=struct('Style','checkbox','String','load only valid','Value',0);
pard.onlyvalid.position=[2,1];
pard.onlyvalid.Width=2;
pard.onlyvalid.TooltipString='Load only localizations with the value tag == true.';

end