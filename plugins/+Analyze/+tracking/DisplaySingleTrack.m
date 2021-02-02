classdef DisplaySingleTrack<interfaces.DialogProcessor
    %Interactive analysis of SPT data
    methods
        function obj=DisplaySingleTrack(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters={'sr_pixrec','numberOfLayers','sr_pos','sr_size','layers','sr_layerson','cam_pixelsize_nm'};
            obj.showresults=true;
        end
        function out=run(obj,p)
            out=[];
            locs=obj.locData.getloc({'xnm','ynm','znm','frame','time','xnmline','ynmline'},'layer',find(obj.getPar('sr_layerson')), 'grouping','ungrouped','Position','roi');
            if ~isempty(locs.time)
                time=locs.time;
                xax='time (xx)';
                ds=quantile(diff(locs.time),0.05);
                
            else
                time=locs.frame;
                xax='frame (xx)';
                ds=quantile(diff(locs.frame),0.05);
            end
            if ~isempty(locs.xnmline)
                x=locs.xnmline;
            else
                x=locs.xnm;
            end
            ax=obj.initaxis('x');
            plot(ax,time,x,'r');
            xlabel(ax,xax)
            ylabel(ax,'position (nm)')
            windowsize=p.filtersize*ds;
            xf=runningWindowAnalysis(time,x,time,windowsize,p.filtermode.selection);
            hold(ax,'on')
            plot(ax,time,xf,'b');
            hold(ax,'off')
            
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end


function pard=guidef

pard.filtert.object=struct('String','filter size','Style','text');
pard.filtert.position=[1,1];

pard.filtersize.object=struct('String','10','Style','edit');
pard.filtersize.position=[1,2];
pard.filtermode.object=struct('String',{{'median','mean'}},'Style','popupmenu');
pard.filtermode.position=[1,3];
% pard.analysismode.Width=3;
% % ,'grid based diffusion coefficients'
% 
% pard.lentrackst.object=struct('String','Minimum length of tracks','Style','text');
% pard.lentrackst.position=[2,1];
% pard.lentrackst.Width=2.5;
% pard.lentrackst.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
%             ' fewer than param.good valid positions.  This is useful \n'...
%            'due to blinking noise particles in the data stream.']);
% pard.lentracks.object=struct('String','2','Style','edit');
% pard.lentracks.position=[2,3.5];
% pard.lentracks.Width=.5;
% pard.lentracks.TooltipString=pard.lentrackst.TooltipString;
% 
% pard.cutoffimmobilet.object=struct('String','Cutoff immobile log10(D(um^2/s))','Style','text');
% pard.cutoffimmobilet.position=[3,1];
% pard.cutoffimmobilet.Width=2.5;
% pard.cutoffimmobilet.TooltipString=sprintf(['set this keyword to eliminate all trajectories with \n'...
%             ' fewer than param.good valid positions.  This is useful \n'...
%            'due to blinking noise particles in the data stream.']);
% pard.cutoffimmobile.object=struct('String','-3','Style','edit');
% pard.cutoffimmobile.position=[3,3.5];
% pard.cutoffimmobile.Width=.5;
% pard.cutoffimmobile.TooltipString=pard.cutoffimmobilet.TooltipString;
% 
% pard.saveD.object=struct('String','Write diffusion coefficient to localization data','Style','checkbox');
% pard.saveD.position=[4,1];
% pard.saveD.Width=3;


pard.plugininfo.description=sprintf('Interactive analysis of SPT data');
pard.plugininfo.type='ProcessorPlugin';
% pard.plugininfo.Name='AnalyzeSPT';
end

