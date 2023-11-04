classdef driftcorrection_dme<interfaces.DialogProcessor
%     1.Cnossen, J., Cui, T. J., Joo, C. & Smith, C. S. Drift correction in
%     localization microscopy using entropy minimization. bioRxiv
%     2021.03.30.437682 (2021) doi:10.1101/2021.03.30.437682.

    methods
        function obj=driftcorrection_dme(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'layer1_','cam_pixelsize_nm'};
                obj.history=true;
                obj.showresults=true;
                obj.guiselector.show=true;
        end
        
        function out=run(obj,p)
            out=[];
            if isempty(obj.locData.loc)
                out.error='no localizations loaded';
                return
            end

            obj.setPar('undoModule','driftDME');
            notify(obj.P,'backup4undo');
            numberOfFiles=obj.locData.files.filenumberEnd;
            layers=find(obj.getPar('sr_layerson'));
            
            %determine which files to correct
            if contains(p.drift_whatfiles.selection,'all')
                files=1:numberOfFiles;
                region='all';
                rmfilter={'filenumber'};
                driftall=true;
            else %visible
                files=p.layer1_.ch_filelist.Value;
                rmfilter={};
                region='roi';
                driftall=false;
            end
            
                for k=files
                    lochere=obj.locData.copy;
                    lochere.files.fileNumberEnd=1;
                    lochere.files.file=lochere.files.file(k);
                    if driftall
                        badind=lochere.loc.filenumber~=k;
                        lochere.removelocs(badind);
                    else
                        badind=false(size(lochere.loc.filenumber));
                    end
                    lochere.regroup;
                    lochere.loc.filenumber=lochere.loc.filenumber*0+1;
                    
%                     locs=lochere.getloc({'frame','xnm','ynm','znm'},'position','all','grouping',groupcheck);
%                      locs=lochere.getloc({'frame','xnm','ynm','znm'},'position','all','grouping',groupcheck,'layer',1,'removeFilter',{'filenumber'});
                    locs=lochere.getloc({'frame','xnm','ynm','znm','locprecnm','locprecznm'},'position',region,'layer',layers,'removeFilter',rmfilter);
                    if length(locs.xnm)<200
                        out.error='Too few localizations. Remove ROI?';
%                         answ=questdlg(['Only ' num2str(length(locs.xnm)/p.drift_timepoints) ' localizations per time window. Abort drift correction?']); %htis is modal: no way to see output.
%                         if ~contains(answ,'No') %not use this
                            return
%                         end
                    end
                        
                    p.maxframeall=max(lochere.loc.frame);
                    p.framestart=min(locs.frame);
                    p.framestop=max(locs.frame);
                    p.roi=obj.locData.files.file(k).info.roi;
                    p.ax=obj.initaxis('driftxyz');
                    [drift,driftinfo,fieldc]=getxyzdrift(locs,p);
                    
                    locsall=copyfields([],lochere.loc,{fieldc{:},'frame','filenumber'});
                    
                    if p.drift_ask %check if to apply drift correction
%                         answ=questdlg('apply drift correction?'); %htis is modal: no way to see output.
                          answ=nonmodaldialog('apply drift correction?','YES, apply', 'NO');
                        if ~contains(answ,'YES') %not use this
                            continue
                        end
                           
                    end
                    
                    locsnew=applydriftcorrection(drift,locsall);
                    lochere.loc=copyfields(lochere.loc,locsnew,fieldc);
                    
%                     lochere.loc=copyfields(lochere.loc,locsnew,{'xnm','ynm'});
                    if isfield(lochere.files.file(1),'driftinfo')
                        driftinfoh=copyfields(lochere.files.file(1).driftinfo,driftinfo);
                    else
                        driftinfoh=driftinfo;
                    end
                    lochere.files.file(1).driftinfo=driftinfoh;
                    if isfield(obj.locData.files.file(k),'driftinfo')
                        driftinfoh=copyfields(obj.locData.files.file(k).driftinfo,driftinfo);
                    else
                        driftinfoh=driftinfo;
                    end
                    obj.locData.files.file(k).driftinfo=driftinfoh;
                    fn=lochere.files.file(1).name;
                    if contains(fn,'_sml.mat')
                        fnn=strrep(fn,'_sml.mat','_dme_sml.mat');
                    end
                    if p.save_dc
                        obj.addhistory;
                        lochere.savelocs(fnn); 
                    end
                    obj.locData.loc.xnm(~badind)=lochere.loc.xnm;
                    obj.locData.loc.ynm(~badind)=lochere.loc.ynm;
                    if isfield(lochere.loc,'znm')
                        obj.locData.loc.znm(~badind)=lochere.loc.znm;
                    end
                end
                obj.locData.regroup;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end



function [drifto,driftinfo,fieldc]=getxyzdrift(locs,p)
isconstc=p.crlbconst;
drifto=[];driftinfo=[];
norm(1)=p.cam_pixelsize_nm(1);
norm(2)=p.cam_pixelsize_nm(end);
norm(3)=1000;
if ~isempty(locs.znm)&&p.correctz
    coords=horzcat(locs.xnm/norm(1),locs.ynm/norm(2),locs.znm/norm(3));
%     crlb=horzcat(locs.locprecnm.^2/norm(1)^2,locs.locprecnm.^2/norm(2)^2,locs.locprecznm.^2/norm(3)^2);
    crlb=horzcat(locs.locprecnm/norm(1),locs.locprecnm/norm(2),locs.locprecznm/norm(3));
%     crlb0=[0.2 0.2 0.2]';
    crlb0=median(crlb)';
    fieldc={'xnm','ynm','znm'};
    isz=1;
else
    coords=horzcat(locs.xnm/norm(1),locs.ynm/norm(2));
    crlb=horzcat(locs.locprecnm/norm(1),locs.locprecnm/norm(2));
%     crlb=horzcat(locs.locprecnm.^2/norm(1)^2,locs.locprecnm.^2/norm(2)^2);
%     crlb0=[0.2 0.2]';
    crlb0=median(crlb)';
    fieldc={'xnm','ynm'};
    isz=0;
    norm=norm(1:2);
end
framenum=int32(locs.frame);
maxframe=max(p.maxframeall);
numspots=length(framenum);
maxit=10000;
drift=zeros(isz+2,maxframe,'single');
framesperbin=p.dme_fbin;
gradientStep=1e-6;
maxdrift=0;
scores=zeros(1,maxit,'single');


maxneighbors=p.maxneighbours;

nIterations =int32([
0;
0;
]);



if isconstc
    crlb=crlb0;
end

isgpu=p.usegpu;

flags=4*isconstc+2*isgpu+isz;

% Flags
% Flags               cuda    cpu
% 2D & constant CRLB  6       4
% 2D & variable CRLB  2       0
% 3D & constant CRLB  7       5
% 3D & variable CRLB  3       1
tic
if isgpu
    dme_cuda(single(coords'), single(crlb'), int32(framenum),...
        numspots, maxit, drift, framesperbin, gradientStep, maxdrift, scores,...
        flags, maxneighbors, nIterations);
else

    dme_cpu(single(coords'), single(crlb'), int32(framenum),...
        numspots, maxit, drift, framesperbin, gradientStep, maxdrift, scores,...
        flags, maxneighbors, nIterations);
end
t=toc;
disp(['Time for DME drift correction: ' num2str(t) 's'])


shift_dme = double(drift');
fcc = [1:maxframe]';
fitrange = round(linspace(1,numel(fcc),min([1000,numel(fcc)])));
shiftavg = movmean(shift_dme,200);
%drift_sm = shiftavg;
drift_sm = [];
for ii = 1:3
   fz = fit(fcc(fitrange),shiftavg(fitrange,ii),'smoothingspline','SmoothingParam',1e-7);
   drift_sm = cat(2,drift_sm,fz(fcc));
end


%drifto.xy.x=drift(1,:)'*norm(1);
%drifto.xy.y=drift(2,:)'*norm(2);

drifto.xy.x=drift_sm(:,1)*norm(1);
drifto.xy.y=drift_sm(:,2)*norm(2);

frames=1:maxframe;
hold(p.ax,'off')
plot(p.ax,frames,drift(1,:)'*norm(1),'color',[1,1,1]*0.7)
hold(p.ax,'on')
plot(p.ax,frames,drifto.xy.x)
plot(p.ax,frames,drift(2,:)'*norm(2)+10,'color',[1,1,1]*0.7)
plot(p.ax,frames,drifto.xy.y+10)


if isz
    %drifto.z=drift(3,:)'*norm(3);
    drifto.z=drift_sm(:,3)*norm(3);
    plot(p.ax,frames,drift(3,:)'*norm(3)-10,'color',[1,1,1]*0.7)
    plot(p.ax,frames,drifto.z-10)
end

xlabel(p.ax,'frames');
ylabel(p.ax,'drift (nm)');
legend(p.ax,'','x','','y','','z','Location','northwest')

% driftinfo=p;
% flags=7;
% dme_cuda(single(coords'), single(crlb'), int32(framenum),...
%     numspots, maxit, drift, framesperbin, gradientStep, maxdrift, scores,...
%  flags, maxneighbors, nIterations);

% if 1
%     finddriftfeatureM(locs,p);
% end
                

   
%     switch p.drift_mirror2c.Value
%         case 1 %all
%             ind=true(length(locs.xnm),1);
%             rep=false;
%             mirror='none';
%             midpoint=0;
%              p.repetitionname='';
%         case 2 %horizontal
%             midpoint=p.cam_pixelsize_nm(1)*(p.roi(1)+(p.roi(1)+p.roi(3))/2);
%             ind=locs.xnm<=midpoint;
%             rep=true;
%             mirror='horizontal';
%              p.repetitionname='1';
%         case 3 %vertical
%             midpoint=p.cam_pixelsize_nm(2)*(p.roi(2)+(p.roi(2)+p.roi(4))/2);
%             ind=locs.ynm<=midpoint;
%             rep=true;
%             mirror='vertical';
%             p.repetitionname='1';
%     end
%     
    
    
%     
%     [driftxy,driftinfoxy]=finddriftfeature(copystructReduce(locs,ind),p);
%     driftinfoxy.mirror=mirror; driftinfoxy.midpoint=midpoint;
%     driftinfo.xy=driftinfoxy;
%     drift.xy=copyfields([],driftxy,{'x','y'});
%     drift.xy(1).mirror=mirror;drift(1).xy(1).midpoint=midpoint;

%     if rep
%          p.repetitionname='2';
%         [driftxy,driftinfoxy]=finddriftfeature(copystructReduce(locs,~ind),p);
%         driftinfoxy.mirror=mirror; driftinfoxy.midpoint=midpoint;
%         driftinfo.xy(2)=(driftinfoxy);
%         drift.xy(2)=copyfields(drift.xy(1),driftxy,{'x','y'});
% %         drift.xy(2).mirror=mirror;drift(2).midpoint=midpoint;
%     end

    
% if ~isempty(locs.znm)&&p.correctz
%     if p.correctxy
%         locsnew=copyfields(locs,applydriftcorrection(drift,locs),{'xnm','ynm'});
%     else
%         locsnew=locs;
%         drift=[];
%     end
%     [driftz,driftinfoz]=finddriftfeatureZ(locsnew,p);
%     
%     drift.z=driftz.z;%copyfields(drift,driftz,'z');
%      driftinfo.z=driftinfoz;%copyfields(driftinfo,driftinfoz);
%      fieldc={'xnm','ynm','znm'};
% else
%     fieldc={'xnm','ynm'};
% end

end

function pard=guidef(obj)


pard.dme_fbint.object=struct('String','frames/bin','Style','text');
pard.dme_fbint.position=[2,1];
pard.dme_fbint.Width=1;

pard.dme_fbin.object=struct('String','100','Style','edit');
pard.dme_fbin.position=[2,2];
pard.dme_fbin.isnumeric=1;
pard.dme_fbin.Width=0.5;

pard.maxneighbourst.object=struct('String','max neighbours','Style','text');
pard.maxneighbourst.position=[3,1];
pard.maxneighbourst.Width=1;

pard.maxneighbours.object=struct('String','10000','Style','edit');
pard.maxneighbours.position=[3,2];
pard.maxneighbours.isnumeric=1;
pard.maxneighbours.Width=0.5;

p(1).value=0; p(1).on={}; p(1).off={'zranget','zrange'};
p(2).value=1; p(2).on=p(1).off; p(2).off={};
pard.correctz.object=struct('String','Correct z-drift','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
pard.correctz.position=[1,3];


pard.zranget.object=struct('String','zrange nm','Style','text','Visible','on');
pard.zranget.position=[3,3];
pard.zranget.Optional=true;

pard.zrange.object=struct('String','-400 400','Style','edit','Visible','on');
pard.zrange.position=[3,4];
pard.zrange.Optional=true;
pard.zrange.Width=0.9;

pard.crlbconst.object=struct('String','use median CRLB','Style','checkbox','Visible','on');
pard.crlbconst.position=[4,1];
pard.crlbconst.Width=1.5;
% pard.crlbconst.Optional=true;

pard.usegpu.object=struct('String','use GPU','Style','checkbox','Visible','on','Value',1);
pard.usegpu.position=[5,1];
% pard.usegpu.Optional=true;

pard.drift_reference.object=struct('String','reference is last frame','Style','checkbox');
pard.drift_reference.position=[7,3];
pard.drift_reference.Optional=true;
pard.drift_reference.object.TooltipString=sprintf('If checked, drift at end of data set is set to zero. \n Useful for sequential acquisition, use this for first data set.');
pard.drift_reference.Width=2;
pard.drift_reference.Optional=true;

% pard.drift_mirror2c.object=struct('String',{{'no mirror', '2 Channels, mirrored, vertical split', '2 Channels, mirrored, horizontal split'}},'Style','popupmenu');
% pard.drift_mirror2c.position=[7,1];
% pard.drift_mirror2c.Optional=true;
% pard.drift_mirror2c.object.TooltipString=sprintf('Vertical split: next to each other, horizontal split: below each other.');
% pard.drift_mirror2c.Width=2;
% pard.drift_mirror2c.Optional=true;


% pard.drift_individual.object=struct('String','correct every file individually','Style','checkbox','Value',1);
% pard.drift_individual.position=[8,1];
% pard.drift_individual.Width=2;
% pard.drift_individual.Optional=true;
pard.drift_whatfiles.object=struct('String',{{'visible','all files'}},'Style','popupmenu','Value',1);
pard.drift_whatfiles.position=[8,1];
pard.drift_whatfiles.Width=1.5;
pard.drift_whatfiles.Optional=true;

pard.drift_ask.object=struct('String','?','Style','checkbox','Value',0);
pard.drift_ask.position=[8,2.6];
pard.drift_ask.Width=.4;
pard.drift_ask.Optional=true;

pard.save_dc.object=struct('String','Save driftcorrected SML','Style','checkbox','Value',1);
pard.save_dc.position=[8,3];
pard.save_dc.Width=2;
pard.save_dc.Optional=true;

pard.plugininfo.name='drift correction DME';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description={'Entropy based drift correction'};
end