classdef driftcorrectionXYZFiducials<interfaces.DialogProcessor
    methods
        function obj=driftcorrectionXYZFiducials(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'layer1_','cam_pixelsize_nm'};
                obj.history=true;
                obj.showresults=true;
                obj.guiselector.show=true;
        end
        
        function out=run(obj,p)
            out=[];
            if ~p.correctxy && ~p.correctz
                return
            end

            
%             if p.beadsource.Value==1 %ROImanager
%                 beads=sites2beads(obj,p);
%             else % segment new
%                 beads=segmentb(obj,p);
%             end
%             %
            
            %batch: many sml files loaded:
                %do it per file, save only file.
                %locData.copy, remove other filenames from .loc, .grouploc,
                %. files.file(k)
                %save copy
                obj.setPar('undoModule','driftfeature');
            notify(obj.P,'backup4undo');
            groupcheck=obj.locData.isgrouped(1);
            numberOfFiles=obj.locData.files.filenumberEnd;
            if p.drift_individual&&numberOfFiles>1
                for k=1:numberOfFiles
                    lochere=obj.locData.copy;
                    lochere.files.fileNumberEnd=1;
                    lochere.files.file=lochere.files.file(k);
                    badind=lochere.loc.filenumber~=k;
                    lochere.removelocs(badind);
                    lochere.regroup;
                    lochere.loc.filenumber=lochere.loc.filenumber*0+1;
                    
%                     locs=lochere.getloc({'frame','xnm','ynm','znm'},'position','all','grouping',groupcheck);
                     locs=lochere.getloc({'frame','xnm','ynm','znm','numberInGroup','groupindex'},'position','all','grouping',groupcheck,'layer',1,'removeFilter',{'filenumber'});
                    p.maxframeall=max(lochere.loc.frame);
                    p.framestart=min(lochere.loc.frame);
                    p.framestop=p.maxframeall;
                    p.roi=obj.locData.files.file(k).info.roi;
                    [drift,driftinfo,fieldc]=getxyzdriftF(lochere,p);
                    
%                     [drift,driftinfo]=finddriftfeature(locs,p);
%                     locsall=copyfields([],lochere.loc,{'xnm','ynm','frame','filenumber'});
                    locsall=copyfields([],lochere.loc,{fieldc{:},'frame','filenumber'});
                    
                    locsnew=applydriftcorrection(drift,locsall);
                    lochere.loc=copyfields(lochere.loc,locsnew,fieldc);
                    
%                     lochere.loc=copyfields(lochere.loc,locsnew,{'xnm','ynm'});
                    lochere.files.file(1).driftinfo=driftinfo;
                    obj.locData.files.file(k).driftinfo=driftinfo;
                    fn=lochere.files.file(1).name;
                    if strfind(fn,'_sml')
                        fnn=strrep(fn,'_sml','_driftc_sml');
                    else
                        fnn=strrep(fn,'fitpos','driftc_sml');
                    end
                    if p.save_dc
                        lochere.savelocs(fnn); 
                    end
                    obj.locData.loc.xnm(~badind)=lochere.loc.xnm;
                    obj.locData.loc.ynm(~badind)=lochere.loc.ynm;
                end
                obj.locData.regroup;
            else
                locs=obj.locData.getloc({'frame','xnm','ynm','znm','numberInGroup','groupindex'},'position','roi','layer',1);
%                  locs=obj.locData.getloc({'frame','xnm','ynm','znm'},'position','fov','grouping',groupcheck);
                if length(locs.xnm)<100
                    locs=obj.locData.getloc({'frame','xnm','ynm','znm','numberInGroup','groupindex'},'position','all','grouping',groupcheck);
                end

            
                p.maxframeall=max(obj.locData.loc.frame);
%                 p.framestart=p.layer1_.frame_min;
                p.framestart=(min(locs.frame));
                p.framestop=max(locs.frame);
                p.roi=obj.locData.files.file(1).info.roi;  
                [drift,driftinfo,fieldc]=getxyzdriftF(obj.locData,p);
                locsall=copyfields([],obj.locData.loc,{fieldc{:},'frame','filenumber'});
                locsnew=applydriftcorrection(drift,locsall);
                obj.locData.loc=copyfields(obj.locData.loc,locsnew,fieldc);
                if length(unique(locsall.filenumber))>1
                    locsnew.channel=locsall.filenumber;
                    obj.locData.loc=copyfields(obj.locData.loc,locsnew,{'channel'});
                end

               obj.locData.files(obj.locData.loc.filenumber(1)).file.driftinfo=driftinfo;
                fn=obj.locData.files(obj.locData.loc.filenumber(1)).file.name;
                            
                if strfind(fn,'_sml')
                    fnn=strrep(fn,'_sml','_driftc_sml');
                else
                    fnn=strrep(fn,'fitpos','driftc_sml');
                end
%                 fnn=strrep(fn,'_sml','_driftc_sml');
                if p.save_dc
                    obj.locData.savelocs(fnn); 
                end
                obj.locData.regroup;
                
            end
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end


function [drift,driftinfo,fieldc]=getxyzdriftF(locData,p)
minframes=1000;

locdc=locData.copy;
dX=50; dT=2;
locdc.regroup(dX,dT);
locg=locdc.getloc({'xnm','ynm','znm','frame','numberInGroup','groupindex'},'grouping','grouped');
loc=locdc.getloc({'xnm','ynm','znm','frame','numberInGroup','groupindex'},'grouping','ungrouped');
ind=locg.numberInGroup>minframes;

%group beads
indf=find(ind);
xb=locg.xnm(ind);
yb=locg.ynm(ind);
groupindex=locg.groupindex(ind);
used=false(length(xb),1);
indbead=1;
for k=1:length(indf)
    if ~used(k)
        nn=find((xb(k)-xb).^2+(yb(k)-yb).^2<dX^2*.001);
        used(nn)=true;
        bead(indbead).xnm=[];
        bead(indbead).znm=[];
        bead(indbead).ynm=[];
        bead(indbead).frame=[];
        for l=1:length(nn)
            indh=loc.groupindex==groupindex(nn(l));
            bead(indbead).xnm=vertcat(bead(indbead).xnm,loc.xnm(indh));
            bead(indbead).ynm=vertcat(bead(indbead).ynm,loc.ynm(indh));
            bead(indbead).znm=vertcat(bead(indbead).znm,loc.znm(indh));
            bead(indbead).frame=vertcat(bead(indbead).frame,loc.frame(indh));
            
        end
        bead(indbead).numpoints=length(bead(indbead).frame);
        indbead=indbead+1;
    end
end
indbad=[bead(:).numpoints]<max(loc.frame)*.95;
bead(indbad)=[];
allt=zeros(max(loc.frame),length(bead),3)+NaN;

for k=1:length(bead)
    allt(bead(k).frame,k,1)=bead(k).xnm-0*nanmedian(bead(k).xnm);
    allt(bead(k).frame,k,2)=bead(k).ynm-0*nanmedian(bead(k).ynm);
    allt(bead(k).frame,k,3)=bead(k).znm-0*nanmedian(bead(k).znm);
end


drift=[];driftinfo=[];
% if 1
%     finddriftfeatureM(locs,p);
% end
                
if p.correctxy
    
    switch p.drift_mirror2c.Value
        case 1 %all
            ind=true(length(locs.xnm),1);
            rep=false;
            mirror='none';
            midpoint=0;
        case 2 %horizontal
            midpoint=p.cam_pixelsize_nm*(p.roi(1)+(p.roi(1)+p.roi(3))/2);
            ind=locs.xnm<=midpoint;
            rep=true;
            mirror='horizontal';
        case 3 %vertical
            midpoint=p.cam_pixelsize_nm*(p.roi(2)+(p.roi(2)+p.roi(4))/2);
            ind=locs.ynm<=midpoint;
            rep=true;
            mirror='vertical';
    end
    [driftxy,driftinfoxy]=finddriftfeature(copystructReduce(locs,ind),p);
    driftinfoxy.mirror=mirror; driftinfoxy.midpoint=midpoint;
    driftinfo.xy=driftinfoxy;
    drift.xy=copyfields([],driftxy,{'x','y'});
    drift.xy(1).mirror=mirror;drift(1).xy(1).midpoint=midpoint;

    if rep
        [driftxy,driftinfoxy]=finddriftfeature(copystructReduce(locs,~ind),p);
        driftinfoxy.mirror=mirror; driftinfoxy.midpoint=midpoint;
        driftinfo.xy(2)=(driftinfoxy);
        drift.xy(2)=copyfields(drift.xy(1),driftxy,{'x','y'});
%         drift.xy(2).mirror=mirror;drift(2).midpoint=midpoint;
    end
end
if ~isempty(locs.znm)&&p.correctz
    if p.correctxy
        locsnew=copyfields(locs,applydriftcorrection(drift,locs),{'xnm','ynm'});
    else
        locsnew=locs;
        drift=[];
    end
    [driftz,driftinfoz]=finddriftfeatureZ(locsnew,p);
    
    drift.z=driftz.z;%copyfields(drift,driftz,'z');
     driftinfo.z=driftinfoz;%copyfields(driftinfo,driftinfoz);
     fieldc={'xnm','ynm','znm'};
else
    fieldc={'xnm','ynm'};
end
% locsall=copyfields([],obj.locData.loc,{fieldc{:},'frame','filenumber'});
% locsnew=applydriftcorrection(drift,locsall);
end

function pard=guidef

% pard.beadsource.object=struct('String',{{'RoiManager','Segment'}},'Style','popupmenu','Value',2);
% pard.beadsource.position=[1,2];
% pard.beadsource.Width=1;

pard.correctxy.object=struct('String','Correct xy-drift','Style','checkbox','Value',1);
pard.correctxy.position=[1,1];

pard.texta.object=struct('String','timepoints','Style','text');
pard.texta.position=[2,1];


pard.drift_timepoints.object=struct('String','10','Style','edit');
pard.drift_timepoints.object.TooltipString=sprintf('whole data is divided into timepoints individual \n blocks. Range: 7-20');
pard.drift_timepoints.position=[2,2];
pard.drift_timepoints.isnumeric=1;

pard.text1.object=struct('String','pixrec nm','Style','text');
pard.text1.position=[3,1];
pard.text1.Optional=true;

pard.drift_pixrec.object=struct('String','10','Style','edit');
pard.drift_pixrec.position=[3,2];
% pard.drift_pixrec.isnumeric=1;
pard.drift_pixrec.object.TooltipString=sprintf('pixel size (nm) for reconstruction. \n Smaller for well defined peak. But slower \n Range: 10-25');
pard.drift_pixrec.Optional=true;

pard.text2.object=struct('String','window pix','Style','text');
pard.text2.position=[4,1];
pard.text2.Optional=true;

pard.drift_window.object=struct('String','7','Style','edit');
pard.drift_window.position=[4,2];
% pard.drift_window.isnumeric=1;
pard.drift_window.object.TooltipString=sprintf('size of region for peakfinding (ellipt. Gaussian). \n should be small, but cover clear maximum. \n Range: 7-15');
pard.drift_window.Optional=true;

pard.text3.object=struct('String','maxdrift nm','Style','text');
pard.text3.position=[5,1];
pard.text3.Optional=true;

pard.drift_maxdrift.object=struct('String','1000','Style','edit');
pard.drift_maxdrift.position=[5,2];

pard.drift_maxdrift.object.TooltipString=sprintf('Maximum drift expected. \n Smaller if data is sparse and wrong peak found. \n larger if no clear peak found. \n Range 250-2000');
pard.drift_maxdrift.Optional=true;

pard.drift_maxpixelst.object=struct('String','max size (pix)','Style','text');
pard.drift_maxpixelst.position=[6,1];
pard.drift_maxpixelst.Optional=true;

pard.drift_maxpixels.object=struct('String','2048','Style','edit');
pard.drift_maxpixels.position=[6,2];
pard.drift_maxpixels.object.TooltipString=sprintf('Maximum size of the reconstructed images. Smaller for speed and lower memory consumption, larger for noisy signal. 128-4096');
pard.drift_maxpixels.Optional=true;


pard.correctz.object=struct('String','Correct z-drift','Style','checkbox','Value',0);
pard.correctz.position=[1,3];

pard.textaz.object=struct('String','timepoints z','Style','text');
pard.textaz.position=[2,3];
pard.textaz.Optional=true;

pard.drift_timepointsz.object=struct('String','10','Style','edit');
pard.drift_timepointsz.object.TooltipString=sprintf('whole data is divided into timepoints individual \n blocks. Range: 10-40');
pard.drift_timepointsz.position=[2,4];
pard.drift_timepointsz.Optional=true;


pard.drift_pixreczt.object=struct('String','z binwidth nm','Style','text');
pard.drift_pixreczt.position=[3,3];
pard.drift_pixreczt.Optional=true;

pard.drift_pixrecz.object=struct('String','5','Style','edit');
pard.drift_pixrecz.position=[3,4];
pard.drift_pixrecz.isnumeric=1;
pard.drift_pixrecz.object.TooltipString=sprintf('pixel size (nm) for reconstruction. \n Smaller for well defined peak. But slower \n Range: 10-25');
pard.drift_pixrecz.Optional=true;

pard.drift_windowzt.object=struct('String','z fit window pix','Style','text');
pard.drift_windowzt.position=[4,3];
pard.drift_windowzt.Optional=true;

pard.drift_windowz.object=struct('String','15','Style','edit');
pard.drift_windowz.position=[4,4];
pard.drift_windowz.isnumeric=1;
pard.drift_windowz.object.TooltipString=sprintf('size of region for peakfinding (ellipt. Gaussian). \n should be small, but cover clear maximum. \n Range: 7-15');
pard.drift_windowz.Optional=true;

pard.zranget.object=struct('String','zrange nm','Style','text');
pard.zranget.position=[5,3];
pard.zranget.Optional=true;

pard.zrange.object=struct('String','-400 400','Style','edit');
pard.zrange.position=[5,4];
pard.zrange.Optional=true;

pard.slicewidtht.object=struct('String','slice width nm','Style','text');
pard.slicewidtht.position=[6,3];
pard.slicewidtht.Optional=true;

pard.slicewidth.object=struct('String','200','Style','edit');
pard.slicewidth.position=[6,4];
pard.slicewidth.Optional=true;



pard.drift_reference.object=struct('String','reference is last frame','Style','checkbox');
pard.drift_reference.position=[7,3];
pard.drift_reference.Optional=true;
pard.drift_reference.object.TooltipString=sprintf('If checked, drift at end of data set is set to zero. \n Useful for sequential acquisition, use this for first data set.');
pard.drift_reference.Width=2;
pard.drift_reference.Optional=true;

pard.drift_mirror2c.object=struct('String',{{'no mirror', '2 Channels, mirrored, horizontal', '2 Channels, mirrored, vertical'}},'Style','popupmenu');
pard.drift_mirror2c.position=[7,1];
pard.drift_mirror2c.Optional=true;
pard.drift_mirror2c.object.TooltipString=sprintf('If checked, drift at end of data set is set to zero. \n Useful for sequential acquisition, use this for first data set.');
pard.drift_mirror2c.Width=2;
pard.drift_mirror2c.Optional=true;


pard.drift_individual.object=struct('String','correct every file individually','Style','checkbox','Value',1);
pard.drift_individual.position=[8,1];
pard.drift_individual.Width=2;
pard.drift_individual.Optional=true;

pard.save_dc.object=struct('String','Save driftcorrected SML','Style','checkbox','Value',1);
pard.save_dc.position=[8,3];
pard.save_dc.Width=2;
pard.save_dc.Optional=true;

pard.plugininfo.name='Fiducial drift correctiom X,Y,Z ';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description={'Drift correction based on cross-correlation.','Algorithm: the data set is divided into [timepoints] blocks, for which superresolution images are calculated. The displacement between all images is calcualted with a FFT-based cross-correlation algorithm. The position of the maxima of the cross-correlation curve are fitted with sub-pixel accuracy with a free elliptical Gaussian.',...
    'A robust estimator is used to calculate the drift vs frame from all pairwise displacements.','All localiaztions visible in the superresolution image are used to infer the drift. Use [Render]...[Layer] to control this.',...
    'If two files are loaded, their drift is calculated together and they are saved as one file with their filenumbers copied to the channel field.',' ','(c) Jonas Ries, EMBL, 2015'};
end