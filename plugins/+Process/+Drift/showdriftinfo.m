classdef showdriftinfo<interfaces.DialogProcessor
    methods
        function obj=showdriftinfo(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
               
        end
        
        function out=run(obj,p)
            out=[];
            driftinfo=obj.locData.files.file(p.dataselect.Value).driftinfo;
            if isfield(driftinfo,'xy')
                for k=1:length(driftinfo.xy)
                    dx=driftinfo.xy(k).dx;
                    dy=driftinfo.xy(k).dy;
                    dxplot=driftinfo.xy(k).dxplot;
                    dyplot=driftinfo.xy(k).dyplot;
                    dxt=driftinfo.xy(k).dxt;
                    dyt=driftinfo.xy(k).dyt;
                    binframes=driftinfo.xy(k).binframes;
                    framesall=1:length(dxt);
                    results_ax2=initaxis(p.resultstabgroup,['di(frame)' num2str(k)]);
                    subplot(1,2,1)
                    hold off
                    plot(dxplot)
                    hold on
                    plot(dx,'k','LineWidth',1.5);
                    sx=(max(dx)-min(dx));
                    ylim([min(dx)-sx/2 max(dx)+sx/2])
                    axis tight

                    subplot(1,2,2)
                    hold off
                    plot(dyplot)
                    hold on
                    plot(dy,'k','LineWidth',1.5);

                    sy=(max(dx)-min(dx));
                    ylim([min(dy)-sy/2 max(dy)+sy/2])
                    axis tight

                    initaxis(p.resultstabgroup,['d(frame) ' num2str(k)]);

                    hold off
                    plot(binframes,dx,'x',framesall,dxt,'k')
                    hold on
                    plot(binframes,dy,'o',framesall,dyt,'r')
                    xlabel('frame')
                    ylabel('dx, dy (nm)')
                    drawnow

                    initaxis(p.resultstabgroup,['dx vs dy ' num2str(k)]);
                    hold off
                    plot(dxt,dyt,'k')
                    hold on
                    plot(dx,dy,'ro')
                    plot(dx(1),dy(1),'gx')
                    xlabel('dx')
                    ylabel('dy')
                    drawnow
                    axis equal
                end
            end
            if isfield(driftinfo,'z')
                for k=1:length(driftinfo.z)
                    dz=driftinfo.z(k).dz;
                    dzplot=driftinfo.z(k).dzplot;
                    dzt=driftinfo.z(k).dzt;
                    
                    binframes=driftinfo.z(k).binframesz;
                    framesall=1:length(dzt);
                    results_ax3=initaxis(p.resultstabgroup,['dzi(frame)' num2str(k)]);
%                     subplot(1,2,1)
                    hold off
                    plot(dzplot)
                    hold on
                    plot(dz,'k','LineWidth',1.5);
                    sz=(max(dz)-min(dz));
                    ylim([min(dz)-sz/2 max(dz)+sz/2])
                    axis tight

%                     subplot(1,2,2)
%                     hold off
%                     plot(dyplot)
%                     hold on
%                     plot(dy,'k','LineWidth',1.5);

%                     sy=(max(dx)-min(dx));
%                     ylim([min(dy)-sy/2 max(dy)+sy/2])
%                     axis tight

                    initaxis(p.resultstabgroup,['dz(frame) ' num2str(k)]);

%                     hold off
                    plot(binframes,dz,'x',framesall,dzt,'k')
%                     hold on
%                     plot(binframes,dy,'o',framesall,dyt,'r')
                    xlabel('frame')
                    ylabel('dz(nm)')
                    drawnow

%                     initaxis(p.resultstabgroup,['dx vs dy ' num2str(k)]);
%                     hold off
%                     plot(dxt,dyt,'k')
%                     hold on
%                     plot(dx,dy,'ro')
%                     plot(dx(1),dy(1),'gx')
%                     xlabel('dx')
%                     ylabel('dy')
%                     drawnow
%                     axis equal
                end
            end
            
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[2,1];
pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.type='ProcessorPlugin';
end