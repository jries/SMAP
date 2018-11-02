classdef Loader_settings<interfaces.DialogProcessor
    methods
        function obj=Loader_settings(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'mainGui'};
        end
        
        function out=load(obj,p,file,mode)
            filedat=load(file);
            if ~isfield(filedat,'parameters')
                errordlg(['no parameters found in file ' file]);
            else
                p.mainGui.setGuiParameters(filedat.parameters,true)
            end
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
            obj.locData.clear('filter');
        end        
    end
end

function loadfile(obj,p,file,mode)
% fobj=obj.locData.files;
% obj.locData.files.filenumberEnd=obj.locData.files.filenumberEnd+1; %write back to .files
filedat=load(file);
filedat.filename=file;

filenumber=obj.locData.files.filenumberEnd;
switch mode
    case 'sml'
        [templocData,GUIsettings,siteexplorer]=load_smlV3(filedat);
    case 'fitpos'
        templocData=loadfitposV2(filedat);
        GUIsettings=[];
        siteexplorer=[];
    case 'sites'
        [templocData,siteexplorer]=load_sites(filedat);
        GUIsettings=[];
    otherwise 
        disp('file format not recognized');
        return;
end
indout=templocData.loc.locprecnm>250|imag(templocData.loc.locprecnm)~=0|isnan(templocData.loc.locprecnm);
templocData.removelocs(indout);

%correct filenumber: .loc, files.filenumber, add files.
nfiles=length(templocData.files.file);
templocData.loc.filenumber=templocData.loc.filenumber+filenumber;
obj.locData.addLocData(templocData);

for k=1:nfiles
%     templocData.files.file(k).number=templocData.files.file(k).number+filenumber;
    templocData.files.file(k).number=k+filenumber;
    if ~isfield(templocData.files.file(k).info,'roi')||isempty(templocData.files.file(k).info.roi)||~isnumeric(templocData.files.file(k).info.roi)
        templocData.files.file(k).info.roi=[0 0 templocData.files.file(k).info.Width templocData.files.file(k).info.Height];
    end
end

newfilenumbers=filenumber+1:filenumber+nfiles;
filestruc=templocData.files.file;
if obj.locData.files.filenumberEnd==0
    obj.locData.files.file=filestruc;
else


obj.locData.files.file(filenumber+1:filenumber+nfiles)=copyfields(obj.locData.files.file(1),filestruc,fieldnames(obj.locData.files.file(1)));
end
obj.locData.files.filenumberEnd=length(obj.locData.files.file);

if ~isempty(GUIsettings) %write back parameters
%     button=questdlg('Restore saved GUI Parameters?','GUI parameters','Yes','No','No');
    if p.updateGuiPar %strcmpi(button,'Yes')
        if isfield(GUIsettings,'par')
            GUIsettings=convertparameters(GUIsettings);
        end
        p.mainGui.setGuiParameters(GUIsettings,true)
    end
end

if ~isempty(siteexplorer)
    se=obj.locData.SE;
    se.empty;
    if siteexplorer.numberOfSites>0
    se.addSites(siteexplorer,newfilenumbers)
    end   
end
end

function pard=guidef
info.name='Gui settings';
info.extensions={'*.mat';'*.*'};
info.dialogtitle='select any GUI settings file or _sml file';
pard.plugininfo=info;
pard.plugininfo.type='LoaderPlugin';
end
