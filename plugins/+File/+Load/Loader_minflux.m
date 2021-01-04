classdef Loader_minflux<interfaces.DialogProcessor
    properties
        loaderpath
    end
    methods
        function obj=Loader_minflux(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file)
            loadfile(obj,p,file)
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
%         function initGui(obj)
%             obj.loaderpath='settings/csvloaderconversion/';
%             files=dir([obj.loaderpath '*.txt']);
%             string={'New format', files(:).name};
%             obj.guihandles.importdef.String=string;
% %             obj.guihandles.importdef.Callback=@obj.csvconvert_callback;
%         end
%         function csvconvert_callback(obj,object,b)
%            if object.Value==1
%                newformat(obj);
%            end
%         end
    end
end


function loadfile(obj,p,file)

l=load(file);
l2=l.position_lists;

ind=0;

for k=1:length(l2)
    locmh=l2{k};
    if isempty(locmh)
        continue
    end
    numl=length(locmh.t);
    loco.frame(ind+1:ind+numl)=locmh.t*100;
    loco.xnm(ind+1:ind+numl)=locmh.x_est_absolute*1000;
    loco.ynm(ind+1:ind+numl)=locmh.y_est_absolute*1000;
    loco.phot(ind+1:ind+numl)=locmh.N;
    loco.validLocalization(ind+1:ind+numl)=single(locmh.validLocalization);
    loco.indFiltered(ind+1:ind+numl)=locmh.indFiltered;
    loco.p0(ind+1:ind+numl)=locmh.p0;
    loco.locprecnm(ind+1:ind+numl)=locmh.crb1D;
    loco.moleculeID1(ind+1:ind+numl)=locmh.moleculeID1;
    loco.moleculeID2(ind+1:ind+numl)=locmh.moleculeID2;
    loco.moleculeID3(ind+1:ind+numl)=locmh.moleculeID3;
    loco.roinumber(ind+1:ind+numl)=locmh.moleculeID3*0+k;
end

loco.xnm=loco.xnm-min(loco.xnm);
loco.ynm=loco.ynm-min(loco.ynm);


indbad=isnan(loco.xnm) | isnan(loco.ynm);

filenumber=obj.locData.files.filenumberEnd+1;

locData=interfaces.LocalizationData;

fn=fieldnames(loco);

for k=1:length(fn)
    fnsmap=fn{k};
    fnimport=loco.(fnsmap)';
    locData.setloc(fnsmap,fnimport(~indbad));
%     if strcmp(fnsmap,'frame')
%         locData.setloc(fnsmap,double(tab.(fnimport)));
%     elseif isfield(tab,fnimport)
%         locData.setloc(fnsmap,single(tab.(fnimport)));
%     end
end

psfnm=100;
zd=zeros(length(locData.loc.frame),1,'single');
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




    pixnm=100;


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


function pard=guidef
info.name='Import MINFLUX';
info.extensions={'*.mat','*.*'};
info.dialogtitle='select any .mat file';
pard.plugininfo=info;  
pard.plugininfo.type='LoaderPlugin';

% pard.importdef.object=struct('Style','popupmenu','String',{{'select import'}});
% pard.importdef.position=[1,1];
% pard.importdef.Width=2;
% pard.importdef.TooltipString='Select definition file for import';
end