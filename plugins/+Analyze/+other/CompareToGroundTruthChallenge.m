classdef CompareToGroundTruthChallenge<interfaces.DialogProcessor
    properties
        induse
    end
    methods
        function obj=CompareToGroundTruthChallenge(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'numberOfLayers','sr_layerson','mainfile','cam_pixelsize_nm'};
%             end   
        end
        
        function out=run(obj,p)
            out=[];
            %get localizations
            [path,file,ext]=fileparts(p.mainfile);
            if ~exist(path,'dir')
                path='settings';
            end
            
            % fit time
            fn=obj.locData.files.file.name;
            l=load(fn);
            fitt=l.saveloc.fitparameters.processfittime;
%             disp(['fit time without loading: ' num2str(fitt,3) ' s']);

            filenew=fullfile(path,['F_' file '_T' num2str(fitt,'%3.0f') '_temp.csv']);
            
            lochere=obj.locData.copy;
            if p.shiftpix
            shiftx=-0.5*p.cam_pixelsize_nm(1);
            shifty=-0.5*p.cam_pixelsize_nm(2);
            else
                shiftx=0;
                shifty=0;
            end
            if p.shiftframe
                lochere.loc.frame=lochere.loc.frame+1;
            end
            lochere.loc.xnm=lochere.loc.xnm+shiftx;
            lochere.loc.ynm=lochere.loc.ynm+shifty;
            lochere.loc.xnm=lochere.loc.xnm+p.offsetxyz(1);
            lochere.loc.ynm=lochere.loc.ynm+p.offsetxyz(2);
            if isfield(lochere.loc,'znm')
            lochere.loc.znm=lochere.loc.znm+p.offsetxyz(3);
            end
             lochere.loc.phot=lochere.loc.phot*p.photonfactor;
            [descfile]=saveLocalizationsCSV(lochere,filenew,p.onlyfiltered,p.numberOfLayers,p.sr_layerson,obj.induse);
            
            %modify challenge data
            challengepath=['External' filesep 'SMLMChallenge' filesep];
            switch p.comparer.Value
                case 1
                    javapath=['"' pwd filesep challengepath 'challenge.jar"'];
                case 2   
                    javapath=['"' pwd filesep challengepath 'CompareLocalization3D.jar"'];
            end
%             
            settingsfile=[challengepath 'CompareLocalizationSettings.txt'];
            
            replacements={'firstRow1','0','shiftY1','0','txtFile1',strrep(filenew,'\','/'),'colY1','3','colX1','2','colZ1','4','shiftX1','0','colF1','1','colI1','5','shiftUnit1','nm','txtDesc1',strrep(descfile,'\','/')};
            modifySettingsFile(settingsfile,replacements{:});
            oldp=pwd;
            cd(challengepath);
            disp('to contintue with Matlab, close SMLMChallenge application');
            if ispc
                system(javapath)
            else
            system(['java -jar ' javapath]) 
            end
            %later fix jave program and call via
            %smlm.assessment.application.Application
            %after adding javaclasspath(javapath)
            cd(oldp)
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function wobble_callback(a,b,obj)
filename=(obj.locData.files.file(1).name);
path=fileparts(filename);
gtfile=[path filesep 'activations.csv'];
if ~exist(gtfile,'file')
    ind=strfind(path,filesep); 
    gtfile=[path(1:ind(end))  'activations.csv'];
end
 
if ~exist(gtfile,'file')
    [fi,pa]=uigetfile(gtfile);
    if fi
    gtfile=[pa fi];
    else
        return
    end
end

gtData=csvread(gtfile,1,0);
% fullnameLoc = get(handles.text5,'String');
%hasHeader = fgetl(fopen(fullnameLoc));
%hasHeader = 1*(sum(isstrprop(hasHeader,'digit'))/length(hasHeader) < .6);
%localData = csvread(fullnameLoc, hasHeader, 0);
%SH: switched to importdata tool and defined columns to make more general
% localData =importdata(fullnameLoc);
% if isstruct(localData)
%     %strip the header
%     localData = localData.data;
% end
% xCol = str2num(get(handles.edit_x,'String'));
% yCol = str2num(get(handles.editY,'String'));
% frCol = str2num(get(handles.editFr,'String'));
% 
% fullnameGT = get(handles.text6,'String');
% %assumes GT file is as defined in competition
% %CSV file. X col 3, y col 4.
% gtData = importdata(fullnameGT);
XCOLGT =3;
YCOLGT =4;
gtAll = gtData(:,[XCOLGT,YCOLGT]);
gt = unique(gtAll,'rows');

% frameIsOneIndexed = get(handles.radiobutton_is1indexed,'Value');
% 
% [pathstr,~,~] = fileparts(fullnameLoc); 
% output_path = pathstr;
% xnm = localData(:,xCol);
% ynm = localData(:,yCol);
% frame = localData(:,frCol);
cam_pixelsize_nm=obj.getPar('cam_pixelsize_nm');
p=obj.getGuiParameters;
if p.shiftpix
shiftx=-0.5*cam_pixelsize_nm(1);
shifty=-0.5*cam_pixelsize_nm(end);
else
    shiftx=0;
    shifty=0;
end
%might be set by the users in future updates
zmin = -750;zmax = 750;zstep = 10;%nm
roiRadius = 500;%nm
frameIsOneIndexed=true;
output_path=path;
wobbleCorrectSimBead(double(obj.locData.loc.xnm+shiftx),double(obj.locData.loc.ynm+shifty),double(obj.locData.loc.frame), gt,zmin,zstep,zmax,roiRadius,frameIsOneIndexed,filename)

% addpath('External/SMLMChallenge')
% wobble_correct;
end

function pard=guidef(obj)
pard.onlyfiltered.object=struct('Style','checkbox','String','Export filtered (displayed) localizations.','Value',1);
pard.onlyfiltered.position=[2,1];
pard.onlyfiltered.Width=2;

pard.forceungrouped.object=struct('Style','checkbox','String','Force ungrouped localiyations','Value',1);
pard.forceungrouped.position=[3,3];
pard.forceungrouped.Width=2;

pard.shiftpix.object=struct('Style','checkbox','String','Shift by 0.5 camera pixels','Value',1);
pard.shiftpix.position=[3,1];
pard.shiftpix.Width=2;




pard.shiftframe.object=struct('Style','checkbox','String','Shift frame by +1','Value',1);
pard.shiftframe.position=[4,1];
pard.shiftframe.Width=4;

pard.comparer.object=struct('Style','popupmenu','String',{{'2D, 2013','3D, 2016'}},'Value',2);
pard.comparer.position=[5,1];
pard.comparer.Width=4;

pard.wobblebutton.object=struct('Style','pushbutton','String','Wobble.m','Callback',{{@wobble_callback,obj}});
pard.wobblebutton.position=[4,4];
pard.wobblebutton.Width=1;


pard.offsetxyzt.object=struct('Style','text','String','Shift by x,y,z nm');
pard.offsetxyzt.position=[6,1];
pard.offsetxyzt.Width=1;
pard.offsetxyz.object=struct('Style','edit','String','0 0 0');
pard.offsetxyz.position=[6,2];
pard.offsetxyz.Width=1;

pard.photonfactort.object=struct('Style','text','String','factor for photon count');
pard.photonfactort.position=[7,1];
pard.photonfactort.Width=1;
pard.photonfactor.object=struct('Style','edit','String','1');
pard.photonfactor.position=[7,2];
pard.photonfactor.Width=1;

pard.plugininfo.type='ProcessorPlugin';
end