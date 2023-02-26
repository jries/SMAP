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
            [loc,info]=decodeh5ToLoc(file);
            loc.filenumber=0*loc.xnm+obj.locData.files.filenumberEnd+1;
%             [~,~,ext]=fileparts(file);
%             switch ext
%                 case '.h5'
%                     [locs,info]=loadh5(file);
%                 case '.csv'
%                     [locs,info]=loadcsv(file);
%             end
%             locData=interfaces.LocalizationData;
%             if strcmp(info.unit,'px')
%                 pix2nm=info.px_size;
%             else
%                 pix2nm=[1 1];
%             end
              filenumber=obj.locData.files.filenumberEnd+1;
%             zd=zeros(size(locs.x),'single');
%             
%             locData.setloc('ynm',single(locs.x*pix2nm(1)));
%             locData.setloc('xnm',single(locs.y*pix2nm(2)));
%             locData.setloc('znm',single(locs.z));
%             locData.setloc('phot',single(locs.phot));
%             locData.setloc('frame',double(locs.frame_ix+1));
%             locData.setloc('prob',single(locs.prob));
%             locData.setloc('LLrel',single(locs.prob));
%             
%             locData.setloc('filenumber',zd+filenumber);
%             locData.setloc('channel',zd);
%             
%             if ~info.thin
%                 locData.setloc('ynmerr',single(locs.x_sig*pix2nm(1)));
%                 locData.setloc('xnmerr',single(locs.y_sig*pix2nm(2)));
%                 locData.setloc('locprecznm',single(locs.z_sig));
%                 locData.setloc('bg',single(locs.bg));
%                 locData.setloc('phot_err',single(locs.phot_sig));
%                 locData.setloc('locprecnm',(locData.loc.xnmerr+locData.loc.ynmerr)/2);
%             else
%                 locData.setloc('bg',zd);
%                 locData.setloc('locprecznm',mean(pix2nm)./sqrt(locs.phot)*3);
%                 locData.setloc('locprecnm',mean(pix2nm)./sqrt(locs.phot));
%             end

            obj.locData.addLocData(loc);
%             filestruc=obj.locData.files.file(end);
%             filestruc=locData.files.file;
            filestruc.name=file;
            mx=ceil(max(loc.xnm)/info.pix2nm(2));
            my=ceil(max(loc.ynm)/info.pix2nm(1));

            filestruc.info=struct('Width',mx,'Height',my,'roi',[0 0 mx my],'cam_pixelsize_um',info.pix2nm([2 1])/1000);
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