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

function [locso,io]=loadcsv(file)
locs=readtable(file,'NumHeaderLines',3);
fid=fopen(file);
fgetl(fid);
l2=fgetl(fid);
l3=fgetl(fid);
fclose(fid);
io.version=sscanf(l2,'# {%*s "%s}');
io.version(end-2:end)='';
io.unit=sscanf(l3,'# {%*s "%s}');
io.unit(end-1:end)=[];
io.px_size=sscanf(l3,'# {%*s %*s %*s [%f, %f }');

varnames=locs.Properties.VariableNames;
for k=1:length(varnames)
    if isnan(locs.(varnames{k})(1))
        if all(isnan(locs.(varnames{k})))
            continue
        end
    end
    locso.(varnames{k})=locs.(varnames{k});
end

% locs=table2struct(locs);
io.thin=false;
if ~isfield(locso, 'x_sig') || all(isnan(locs.x_sig))
    io.thin=true;
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