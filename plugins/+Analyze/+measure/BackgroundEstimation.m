classdef BackgroundEstimation<interfaces.DialogProcessor
    % Calculates average background per pixel. For the SlowStorm paper
    properties
    end
    methods
        function obj=BackgroundEstimation(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=false; %set true, if results are shown by default
            obj.history=false; %if set true, every time the plugin is called, its parameters are added to the history. Makes sense only if plugin changes the data
        end       
        function initGui(obj)
%             obj.guihandles.locfield=obj.getPar('locFields');
        end
        function out=run(obj,p)
            layers=find(obj.getPar('sr_layerson'));
            for k=1:length(layers)
                locsu=obj.locData.getloc({'bg','numberInGroup'},'position','roi','layer',layers(k),'grouping','ungrouped');
                locsg=obj.locData.getloc({'bg','numberInGroup','phot'},'position','roi','layer',layers(k),'grouping','grouped');
                locsgall=obj.locData.getloc({'phot'},'position','roi','grouping','grouped');
                out.bg_ungrouped(k)=mean(locsu.bg);
                out.bg_grouped_tot(k)=mean(locsg.bg);
                out.bg_grouped_frame(k)=mean(locsg.bg./locsg.numberInGroup);
                dphot=10;maxphot=1e5;
                photrange=0:dphot:maxphot;
                hc=histcounts(locsg.phot,photrange);
                hcall=histcounts(locsgall.phot,photrange);
                photrange(end)=[]; photrange=photrange+dphot/2;
                out.photonrange=photrange;
                out.phothistfiltered=hc;
                out.phothistall=hcall;
            end

        end
        function pard=guidef(obj)
            pard.text.object=struct('String','calculate background', 'Style','text');
            pard.text.position=[1,1];
            pard.text.Width=4;
            pard.plugininfo.type='ProcessorPlugin';
            pard.plugininfo.description='Calculates average background per pixel. For the SlowStorm paper';
         
        end
    end
end


function callbackfunction(uiobject,data,obj,extradata)
disp('callback')
end

function aftersync_callback(obj)
disp('aftersync')
end