classdef annotatePeaks<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
        peakline
        roihandle
        axis
    end
    methods
        function obj=annotatePeaks(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
        %draw 
        modality=p.modality.selection;
        obj.axis=obj.setoutput('profile');
        fs=obj.site.evaluation.fibrilStatistics.measurement;
        switch modality
            case 'deviation'
                dev=fs.deviation.value;
                dsmooth=fs.deviation.value_smo;
            case 'polarization'
                dev=fs.P.value;
                dsmooth=fs.P.value_smo;
        end
        
        posx=(1:length(dev))'*10;
        hold(obj.axis,'off');
        plot(obj.axis,posx,dev,'-')
        hold(obj.axis,'on');
        plot(obj.axis,posx,dsmooth,'r-','LineWidth',3);
        if isfield(obj.site.evaluation,obj.name)
            out=obj.site.evaluation.(obj.name);
        else
            out.(modality).Position=[];
        end
        if isfield(obj.site.evaluation,obj.name) && isfield(obj.site.evaluation.(obj.name),modality) && ~isempty(obj.site.evaluation.(obj.name).(modality).Position)
            obj.roihandle=images.roi.Polyline(obj.axis,'Position',obj.site.evaluation.(obj.name).(modality).Position);
            out=obj.site.evaluation.(obj.name);
            addlistener(obj.roihandle,'ROIMoved',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexAdded',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexDeleted',@(src,evt) obj.updateposition(src,evt));
            obj.plotdistances;
        end

        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function select_callback(obj,a,b)
            if ~isempty(obj.roihandle)&&isvalid(obj.roihandle)
                delete(obj.roihandle);
            end
            obj.roihandle=drawpolyline;
%             obj.site.evaluation.(obj.name).Position=obj.roihandle.Position;
            addlistener(obj.roihandle,'ROIMoved',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexAdded',@(src,evt) obj.updateposition(src,evt));
            addlistener(obj.roihandle,'VertexDeleted',@(src,evt) obj.updateposition(src,evt));
            updateposition(obj,a,b)
        end
        function updateposition(obj,a,b)
            if obj.getSingleGuiParameter('equaldistance')
                pos=obj.roihandle.Position;
                pos(:,1)=linspace(pos(1,1),pos(end,1),size(pos,1));
                obj.roihandle.Position=pos;
            end
            modality=obj.getSingleGuiParameter('modality').selection;
            obj.site.evaluation.(obj.name).(modality).Position=obj.roihandle.Position;
            obj.plotdistances;
        end
        function plotdistances(obj)
            modality=obj.getSingleGuiParameter('modality').selection;
            pos=obj.site.evaluation.(obj.name).(modality).Position;
            
            period2=mean(diff(pos(:,1)));
            form='%2.2f';
            period=(pos(end,1)-pos(1,1))/(size(pos,1)-1);
            
            title(['Period: ' num2str(period,form)])
            
            obj.site.evaluation.(obj.name).(modality).Period=period;
            
        end
%         function fit_callback(obj,a,b)
%             pos=obj.site.evaluation.(obj.name).Position;
%             period=mean(diff(pos(:,1)));
%             
%             
%         end
    end

end




function pard=guidef(obj)
pard.drawline.object=struct('Style','pushbutton','String','select peaks','Callback',@obj.select_callback);
pard.drawline.position=[1,1];
pard.drawline.Width=2;

pard.modality.object=struct('Style','popupmenu','String',{{'deviation','polarization'}});
pard.modality.position=[1,3];
pard.modality.Width=2;

pard.equaldistance.object=struct('Style','checkbox','String','equal distances');
pard.equaldistance.position=[2,1];
pard.equaldistance.Width=2;

% pard.fitpeaks.object=struct('Style','pushbutton','String','fit peaks','Callback',@obj.fit_callback);
% pard.fitpeaks.position=[1,1];
% pard.fitpeaks.Width=2;



% pard.dxt.Width=3;
% pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end
