classdef driftcorrectionM<interfaces.DialogProcessor
%     As DriftcorrectionXYZ, but additionally changes the magnification
%     over time. This can sometimes locally correct the drift in live-cell
%     imaging or in case the sample deforms.
    methods
        function obj=driftcorrectionM(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'layer1_'};
                obj.history=true;
                obj.showresults=true;
                obj.guiselector.show=true;
        end
        
        function out=run(obj,p)
            out=[];
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
                     locs=lochere.getloc({'frame','xnm','ynm','znm'},'position','all','grouping',groupcheck,'layer',1,'removeFilter',{'filenumber'});
                    p.maxframeall=max(lochere.loc.frame);
                    p.framestart=min(lochere.loc.frame);
                    p.framestop=p.maxframeall;
                    [drift,driftinfo,fieldc]=getMdrift(locs,p);
                    
%                     [drift,driftinfo]=finddriftfeature(locs,p);
%                     locsall=copyfields([],lochere.loc,{'xnm','ynm','frame','filenumber'});
                    locsall=copyfields([],lochere.loc,{fieldc{:},'frame','filenumber'});
                    
                    locsnew=applydriftcorrectionM(drift,locsall);
                    lochere.loc=copyfields(lochere.loc,locsnew,fieldc);
                    
%                     lochere.loc=copyfields(lochere.loc,locsnew,{'xnm','ynm'});
                    lochere.files.file(1).driftinfo=driftinfo;
                    obj.locData.files.file(k).driftinfo=driftinfo;
                    fn=lochere.files.file(1).name;
                    if strfind(fn,'_sml')
                        fnn=strrep(fn,'_sml','_driftcM_sml');
                    else
                        fnn=strrep(fn,'fitpos','driftcM_sml');
                    end
                    if p.save_dc
                        lochere.savelocs(fnn); 
                    end
                    obj.locData.loc.xnm(~badind)=lochere.loc.xnm;
                    obj.locData.loc.ynm(~badind)=lochere.loc.ynm;
                end
                obj.locData.regroup;
            else
                locs=obj.locData.getloc({'frame','xnm','ynm','znm'},'position','roi','layer',1);
%                  locs=obj.locData.getloc({'frame','xnm','ynm','znm'},'position','fov','grouping',groupcheck);
                if length(locs.xnm)<100
                    locs=obj.locData.getloc({'frame','xnm','ynm','znm'},'position','all','grouping',groupcheck);
                end

            
                p.maxframeall=max(obj.locData.loc.frame);
                p.framestart=p.layer1_.frame_min;
                p.framestart=max(min(locs.frame),(p.layer1_.frame_min));
                p.framestop=min(p.layer1_.frame_max,p.maxframeall);
                [drift,driftinfo,fieldc]=getMdrift(locs,p);
                locsall=copyfields([],obj.locData.loc,{fieldc{:},'frame','filenumber'});
                midpoint=[mean([max(locs.xnm),min(locs.xnm)]) mean([max(locs.ynm),min(locs.ynm)])];
                locsnew=applydriftcorrectionM(drift,locsall,midpoint);
                obj.locData.loc=copyfields(obj.locData.loc,locsnew,fieldc);
                if length(unique(locsall.filenumber))>1
                    locsnew.channel=locsall.filenumber;
                    obj.locData.loc=copyfields(obj.locData.loc,locsnew,{'channel'});
                end

               obj.locData.files(obj.locData.loc.filenumber(1)).file.driftinfo=driftinfo;
                fn=obj.locData.files(obj.locData.loc.filenumber(1)).file.name;
                              
                if strfind(fn,'_sml')
                    fnn=strrep(fn,'_sml','_driftcM_sml');
                else
                    fnn=strrep(fn,'fitpos','driftcM_sml');
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


function [drift,driftinfo,fieldc]=getMdrift(locs,p)

[drift,driftinfo]=finddriftfeatureM(locs,p);

fieldc={'xnm','ynm'};

end

function poso=applydriftcorrectionM(drift,pos,midpoint)
if nargin<3
    midpoint=[mean(max(pos.xnm),min(pos.xnm)) mean(max(pos.xnm),min(pos.xnm))];
end
indg=pos.frame>0;
if isfield(drift,'M') 
    poso.xnm=pos.xnm;
    poso.ynm=pos.xnm;
    poso.xnm(indg)=(pos.xnm(indg)-midpoint(1)).*drift.M(pos.frame(indg))+midpoint(1);
    poso.ynm(indg)=(pos.ynm(indg)-midpoint(2)).*drift.M(pos.frame(indg))+midpoint(2);
end
end

function pard=guidef


pard.texta.object=struct('String','timepoints','Style','text');
pard.texta.position=[2,1];


pard.drift_timepoints.object=struct('String','10','Style','edit');
pard.drift_timepoints.object.TooltipString=sprintf('whole data is divided into timepoints individual \n blocks. Range: 7-20');
pard.drift_timepoints.position=[2,2];
pard.drift_timepoints.isnumeric=1;

pard.text1.object=struct('String','pixrec nm','Style','text');
pard.text1.position=[3,1];
pard.text1.Optional=true;
pard.text1.Width=1.5;

pard.drift_pixrec.object=struct('String','150','Style','edit');
pard.drift_pixrec.position=[3,2.5];
pard.drift_pixrec.Width=0.5;
% pard.drift_pixrec.isnumeric=1;
pard.drift_pixrec.object.TooltipString=sprintf('pixel size (nm) for reconstruction. \n Smaller for well defined peak. But slower \n Range: 10-25');
pard.drift_pixrec.Optional=true;

pard.text2.object=struct('String','Magnification search range','Style','text');
pard.text2.position=[4,1];
pard.text2.Optional=true;
pard.text2.Width=1.5;
pard.drift_window.object=struct('String','.02','Style','edit');
pard.drift_window.position=[4,2.5];
% pard.drift_window.isnumeric=1;
pard.drift_window.object.TooltipString=sprintf('maximum change in magnification from one block to next');
pard.drift_window.Optional=true;
pard.drift_window.Width=0.5;



pard.drift_reference.object=struct('String','reference is last frame','Style','checkbox');
pard.drift_reference.position=[7,1];
pard.drift_reference.Optional=true;
pard.drift_reference.object.TooltipString=sprintf('If checked, drift at end of data set is set to zero. \n Useful for sequential acquisition, use this for first data set.');
pard.drift_reference.Width=2;
pard.drift_reference.Optional=true;

pard.drift_individual.object=struct('String','correct every file individually','Style','checkbox','Value',1);
pard.drift_individual.position=[8,1];
pard.drift_individual.Width=2;
pard.drift_individual.Optional=true;

pard.save_dc.object=struct('String','Save driftcorrected SML','Style','checkbox','Value',1);
pard.save_dc.position=[8,3];
pard.save_dc.Width=2;
pard.save_dc.Optional=true;

pard.plugininfo.name='drift correction Magnification';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description={' As DriftcorrectionXYZ, but additionally changes the magnification over time. This can sometimes locally correct the drift in live-cell imaging or in case the sample deforms. Needs modification for M. Drift correction based on cross-correlation.','Algorithm: the data set is divided into [timepoints] blocks, for which superresolution images are calculated. The displacement between all images is calcualted with a FFT-based cross-correlation algorithm. The position of the maxima of the cross-correlation curve are fitted with sub-pixel accuracy with a free elliptical Gaussian.',...
    'A robust estimator is used to calculate the drift vs frame from all pairwise displacements.','All localiaztions visible in the superresolution image are used to infer the drift. Use [Render]...[Layer] to control this.',...
    'If two files are loaded, their drift is calculated together and they are saved as one file with their filenumbers copied to the channel field.',' ','(c) Jonas Ries, EMBL, 2015'};
end