classdef driftcorrection_beads<interfaces.DialogProcessor
    methods
        function obj=driftcorrection_beads(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                % obj.inputParameters={'layer1_','cam_pixelsize_nm'};
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
            if numberOfFiles>1
                warndlg('please use only on a single file at a time')
                return
            end
            [dx,dy,dz,frames]=getbeads(obj,p,'');
            
            
            maxframe=max(obj.locData.loc.frame);
            framerange=(1:maxframe)';
            drift.xy.x=runningWindowAnalysis(vertcat(frames{:}),vertcat(dx{:}),framerange,p.filterwin,'mean');
            drift.xy.y=runningWindowAnalysis(vertcat(frames{:}),vertcat(dy{:}),framerange,p.filterwin,'mean');
            
            axd=obj.initaxis('drift'); 
            plot(axd, framerange, drift.xy.x,framerange,drift.xy.y)

            if ~isempty(dz)
                drift.xy.z=runningWindowAnalysis([frames{:}],[dz{:}],framerange,p.filterwin,'mean');
                hold(axd,'on')
                plot(axd, framerange, drift.xy.z)
                hold(axd,'off')
            end

                 
            
                        
                    % p.maxframeall=max(lochere.loc.frame);
                    % p.framestart=min(locs.frame);
                    % p.framestop=max(locs.frame);
                    % p.roi=obj.locData.files.file(k).info.roi;
                    % p.ax=obj.initaxis('driftxyz');
                    % [drift,driftinfo,]=getxyzdrift(locs,p);
                    
                    % locsall=copyfields([],obj.locData.loc,{fieldc{:},'frame','filenumber'});

                    
                    locsnew=applydriftcorrection(drift,obj.locData.loc);
                    % locall=copyfields(obj.locData.loc,locsnew);
                    obj.locData.loc.xnm=locsnew.xnm;
                    obj.locData.loc.ynm=locsnew.ynm;
                     if ~isempty(dz)
                         obj.locData.loc.znm=locsnew.znm;
                     end

                    [dx,dy,dz,frames]=getbeads(obj,p,'_corr'); %plot again
                    
% %                     lochere.loc=copyfields(lochere.loc,locsnew,{'xnm','ynm'});
%                     if isfield(lochere.files.file(1),'driftinfo')
%                         driftinfoh=copyfields(lochere.files.file(1).driftinfo,driftinfo);
%                     else
%                         driftinfoh=driftinfo;
%                     end
%                     lochere.files.file(1).driftinfo=driftinfoh;
%                     if isfield(obj.locData.files.file(k),'driftinfo')
%                         driftinfoh=copyfields(obj.locData.files.file(k).driftinfo,driftinfo);
%                     else
%                         driftinfoh=driftinfo;
%                     end
                    % obj.locData.files.file(k).driftinfo=driftinfoh;
                    fn=obj.locData.files.file(1).name;
                    if contains(fn,'_sml.mat')
                        fnn=strrep(fn,'_sml.mat','_driftbeads_sml.mat');
                    end
                    if p.save_dc
                        obj.addhistory;
                        obj.locData.savelocs(fnn); 
                    end
                    % obj.locData.loc.xnm(~badind)=lochere.loc.xnm;
                    % obj.locData.loc.ynm(~badind)=lochere.loc.ynm;
                    % if isfield(lochere.loc,'znm')
                    %     obj.locData.loc.znm(~badind)=lochere.loc.znm;
                    % end
                
                obj.locData.regroup;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end



function [drifto,driftinfo,fieldc]=getxyzdrift(locs,p)


end

function [dx,dy,dz,frames,isz]=getbeads(obj,p,suff)
if contains(p.beadsource.selection,'View')
    lochere=obj.locData.getloc({'xnm','ynm','znm','frame','numberInGroup','groupindex'},'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','ungrouped');
    if p.beadsource_minlocsuse
        beadid=unique(lochere.groupindex(lochere.numberInGroup>p.beadsource_minlocs));
    else
        beadid=unique(lochere.groupindex);
    end
else
    sites=obj.locData.SE.sites;
    for k=length(sites):-1:1
        lochere=obj.locData.getloc({'xnm','ynm','znm','frame','numberInGroup','groupindex'},'layer',find(obj.getPar('sr_layerson')),'position',sites(k),'grouping','ungrouped');
        beadid(k)=mode(lochere.groupindex);
    end
    
end
lochere=obj.locData.loc;
    axx1=obj.initaxis(['x' suff]); hold(axx1,'off')
    axy1=obj.initaxis(['y' suff]); hold(axy1,'off')
    if isfield(lochere,'znm') && ~isempty(lochere.znm)
        axz1=obj.initaxis(['z' suff]); hold(axz1,'off')
    end
    for k=1:length(beadid)
        indh=lochere.groupindex==beadid(k);
        if sum(indh)<p.beadsource_minlocs
            continue
        end
        plot(axx1,lochere.frame(indh),lochere.xnm(indh)-mean(lochere.xnm(indh)));
        hold(axx1,'on')
        plot(axy1,lochere.frame(indh),lochere.ynm(indh)-mean(lochere.ynm(indh)));
        hold(axy1,'on')
        dx{k}=lochere.xnm(indh)-mean(lochere.xnm(indh));
        dy{k}=lochere.ynm(indh)-mean(lochere.ynm(indh));
        if isfield(lochere,'znm') && ~isempty(lochere.znm)
            dz{k}=lochere.znm(indh)-mean(lochere.znm(indh));
            plot(axz1,lochere.frame(indh),lochere.znm(indh)-mean(lochere.znm(indh)));
            hold(axz1,'on')
            xlabel(axz1,'frames'); ylabel(axz1,'y (nm)')
        else
            dz=[];
        end
        frames{k}=lochere.frame(indh);
    end

    xlabel(axx1,'frames'); ylabel(axx1,'x (nm)')
    xlabel(axy1,'frames'); ylabel(axy1,'y (nm)')
end

function pard=guidef(obj)
pard.beadsourcet.object=struct('String','Beads from ','Style','text');
pard.beadsourcet.position=[2,1];
pard.beadsource.object=struct('String',{{'View/Roi','Roi manager'}},'Style','popupmenu');
pard.beadsource.position=[2,2];
pard.beadsource.Width=1.5;
pard.beadsource_minlocsuse.object=struct('String','#locs >','Style','checkbox');
pard.beadsource_minlocsuse.position=[3,1];
pard.beadsource_minlocs.object=struct('String','100','Style','edit');
pard.beadsource_minlocs.position=[3,2];
pard.beadsource_minlocs.Width=0.5;


pard.filtert.object=struct('String','Filter window (frames)','Style','text');
pard.filtert.position=[4,1];
pard.filterwin.object=struct('String','5','Style','edit');
pard.filterwin.position=[4,2];
pard.filterwin.Width=0.5;

% pard.drift_reference.object=struct('String','reference is last frame','Style','checkbox');
% pard.drift_reference.position=[7,3];
% pard.drift_reference.Optional=true;
% pard.drift_reference.object.TooltipString=sprintf('If checked, drift at end of data set is set to zero. \n Useful for sequential acquisition, use this for first data set.');
% pard.drift_reference.Width=2;
% pard.drift_reference.Optional=true;

% pard.drift_whatfiles.object=struct('String',{{'visible','all files'}},'Style','popupmenu','Value',1);
% pard.drift_whatfiles.position=[8,1];
% pard.drift_whatfiles.Width=1.5;
% pard.drift_whatfiles.Optional=true;

% pard.drift_ask.object=struct('String','?','Style','checkbox','Value',0);
% pard.drift_ask.position=[8,2.6];
% pard.drift_ask.Width=.4;
% pard.drift_ask.Optional=true;

pard.save_dc.object=struct('String','Save driftcorrected SML','Style','checkbox','Value',1);
pard.save_dc.position=[8,3];
pard.save_dc.Width=2;
pard.save_dc.Optional=true;

pard.plugininfo.name='drift correction Beads';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description={'Fiducial based drift correction'};
end