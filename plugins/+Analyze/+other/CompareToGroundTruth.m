classdef CompareToGroundTruth<interfaces.DialogProcessor
%     Compares fitted data of simulated images to ground truth and
%     calculates several quality metrices
    methods
        function obj=CompareToGroundTruth(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'numberOfLayers','sr_layerson','mainfile','cam_pixelsize_nm'};
%             end   
        end
        
        function out=run(obj,p)
            if p.checklayer
                % do something with grouping? at least warn?
                obj.locData.sort('frame');
                obj.locData.filter;
                lr=obj.locData.getloc({'xnm','ynm','znm','phot','frame'},'layer',p.reflayer.Value);
                lt=obj.locData.getloc({'xnm','ynm','znm','phot','frame'},'layer',p.targetlayer.Value);
%                 lt.frame=lt.frame+1;
                [outlayer2D, outlayer3D,mr,mt]=getmatch(lr,lt,p);
                outlayer2D.name='layer2D';
                outlayer3D.name='layer3D';
            end
            if p.checkdata
            end
            tab=maketable(outlayer2D,outlayer3D);
            ax=obj.initaxis('results');
            f=ax.Parent;
            delete(ax);
%             pos=f.Position;pos(1:2)=20;pos(3:4)=pos(3:4)*.9;
            uit=uitable(f,'Units','normalized','Position',[0 0 1 0.9]);
            uit.Units='pixels';
            uit.Data=table2cell(tab);
            uit.ColumnName=tab.Properties.VariableNames;
            ncol=length(uit.ColumnName);
            w=uit.Position(3)/(ncol+1);
            for k=1:ncol
              cw{k}=w;
            end
            uit.ColumnWidth=cw;
            
            ax=obj.initaxis('z compare gt');
            plot(ax,lr.znm(mr),lt.znm(mt),'.');
            zrh=lr.znm(mr);
            zth=lt.znm(mt);
            inrange=abs(zrh)<300;
            zrh=zrh(inrange);
            zth=zth(inrange);
            B=[zrh ones(length(zrh),1)];
            fit=B\zth;
            slope=fit(1);off=fit(2);
            hold(ax,'on');
            plot(ax, zrh,slope*zrh+off);
            hold(ax,'off');
            title(ax,['slope: ' num2str(slope) ', off: ' num2str(off)])
            xlabel(ax,'z reference (nm)')
            ylabel(ax,'z target(nm)')
            
            
            out=[];
            ax=obj.initaxis(' phot');
            photr=lr.phot(mr);photr=photr(inrange);
            phott=lt.phot(mt);phott=phott(inrange);
            plot(ax,lr.phot(mr),lt.phot(mt),'.');
            plot(ax,photr,phott,'r.');
            xlabel('reference')
            ylabel('target')
            title(['reference/target' num2str(nanmedian(photr./phott))])
            
            ax=obj.initaxis('phot (z)');
            plot(ax,lr.znm(mr),lr.phot(mr)./lt.phot(mt),'.')
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function [out,out3D,mr,mt]=getmatch(lr,lt,p)
if length(p.searchradius)==1
    p.searchradius(2)=p.searchradius(1);
end
    [mr,mt,ur,ut]=matchlocs(lr.xnm,lr.ynm,lt.xnm,lt.ynm,[],p.searchradius(1));
    lr.x=lr.xnm;lr.y=lr.ynm;lt.x=lt.xnm;lt.y=lt.ynm;
    [mr,mt,ur,ut,nseen]=matchlocsall(lr,lt,0,0,p.searchradius(1));
     out=getmatchstat(lr,lt,mr,mt,ur,ut);  
    
    %look at z
    if isfield(lr,'znm')&&~isempty(lr.znm)&&isfield(lt,'znm')&&~isempty(lt.znm)
        dz=lr.znm(mr)-lt.znm(mt);
        inz=abs(dz)<=p.searchradius(2);
        
        ur2=vertcat(ur',mr(~inz));
        ut2=vertcat(ut',mt(~inz));
        mr2=mr(inz);
        mt2=mt(inz);
    end
    out3D=getmatchstat(lr,lt,mr2,mt2,ur2,ut2);  
end

function out=getmatchstat(lr,lt,mr,mt,ur,ut)
    nrtot=length(mr)+length(ur);
    nttot=length(mt)+length(ut);
    out.falsePositives=length(ut)/nrtot;
    out.falseNegatives=length(ur)/nrtot;
    
    dx=lr.xnm(mr)-lt.xnm(mt);
    dy=lr.ynm(mr)-lt.ynm(mt);
    if isfield(lr,'znm')&&~isempty(lr.znm)&&isfield(lt,'znm')&&~isempty(lt.znm)
        dz=lr.znm(mr)-lt.znm(mt);
    else
        dz=0*dx;
    end
    out.shift=[mean(dx),mean(dy),mean(dz)];
    dx=dx-mean(dx);dy=dy-mean(dy);dz=dz-mean(dz);
    out.stderr=[std(dx),std(dy),std(dz)];
    out.err=sqrt([mean(dx.^2),mean(dy.^2),mean(dz.^2)]);
    out.meanabs=([mean(abs(dx)),mean(abs(dy)),mean(abs(dz))]);

end

function tab=maketable(varargin)
for k=1:length(varargin)
    inh=varargin{k};
    ss.name{k,1}=inh.name;
    ss.fPos(k,1)=inh.falsePositives;
    ss.fNeg(k,1)=inh.falseNegatives;
    ss.Jaccard(k,1)=(1-inh.falseNegatives)*(1-inh.falsePositives);
    ss.stdxy(k,1)=sqrt(sum(inh.stderr(1:2).^2));
    ss.meanabsxy(k,1)=mean(inh.meanabs);
    ss.errxy(k,1)=sqrt(sum(inh.err(1:2).^2));
    ss.stdz(k,1)=(inh.stderr(3));
    ss.meanabsz(k,1)=mean(inh.meanabs(3));
    ss.errz(k,1)=(inh.err(3));
end
tab=struct2table(ss);

end
function wobble_callback(a,b,obj)
filename=(obj.locData.files.file(1).name);
path=fileparts(filename);
gtfile=[path filesep 'activations.csv','load ground truth for data'];
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

gtData=csvread(gtfile);
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
shiftx=-0.5*cam_pixelsize_nm;
shifty=-0.5*cam_pixelsize_nm;
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

pard.t1.object=struct('Style','text','String','Reference');
pard.t1.position=[1,1.5];
pard.t2.object=struct('Style','text','String','Target');
pard.t2.position=[1,2.5];

pard.t3.object=struct('Style','text','String','search radius x,y (nm)');
pard.t3.position=[1,3.5];
pard.t3.Width=1;
pard.searchradius.object=struct('Style','edit','String','100 300');
pard.searchradius.position=[1,4.5];
pard.searchradius.Width=.5;


pard.checklayer.object=struct('Style','checkbox','String','','Value',1);
pard.checklayer.position=[2,1];
pard.checklayer.Width=.5;

pard.reflayer.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',1);
pard.reflayer.position=[2,1.5];
pard.reflayer.Width=1;

pard.targetlayer.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',2);
pard.targetlayer.position=[2,2.5];
pard.targetlayer.Width=1;

pard.checkdata.object=struct('Style','checkbox','String','','Value',0);
pard.checkdata.position=[3,1];
pard.checkdata.Width=.5;

pard.refdata.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',1);
pard.refdata.position=[3,1.5];
pard.refdata.Width=1;

pard.targetdata.object=struct('Style','popupmenu','String',{{'layer1','layer2','layer3','layer4','layer5'}},'Value',2);
pard.targetdata.position=[3,2.5];
pard.targetdata.Width=1;


pard.syncParameters={{'filelist_short','refdata',{'String'}},{'filelist_short','targetdata',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Compares fitted data of simulated images to ground truth and calculates several quality metrices';
end