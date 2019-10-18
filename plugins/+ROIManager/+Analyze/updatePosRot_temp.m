classdef updatePosRot_temp<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=updatePosRot_temp(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'se_viewer'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            out =[];
            %% Int.
%             p.relativeAmpLimit = 3;         % max acceptable relative amplitude
%             p.rotLimit = 44;                % max angle of rotation
%             p.applyFilter = 0;              % filter out bad (extreme) site or not
%             p.viz = 0;                      % show plotSElink
            
            %% Sort sites according to 'use'
            fdcal=figure(519);
            dcal=plugin('ROIManager','Analyze','SortROIs',fdcal,obj.P);
            dcal.attachLocData(obj.SE.locData);
            dcal.makeGui;
            pp=dcal.getGuiParameters;
            pp.sortprop1.Value=2;                               % select use
            pp.sortedit1 = 'annotation.use';                    % select use
            pp.direction1.Value=2;                              % decend
            dcal.setGuiParameters(pp);
            dcal.processgo;
            pause(10e-1000);
            
            %% to be added
            p.oldversion = 0;
            p.rmBound = 1;
            
            %% Main
            %% Filtering out sites with extreme parameters
            lastUsed = find(getFieldAsVector(obj.locData.SE.sites,'annotation.use'),1, 'last');
            
            % overwite parameters for unused sites by 0
            for k=lastUsed+1:obj.locData.SE.numberOfSites
                obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters=zeros([1 13]);
            end
            
            ax=newAxes([], 'Name','Log10(amp) to height(filtering displayed)','NumberTitle','off');
            ax2=newAxes([], 'Name','Rotation to height','NumberTitle','off');
            distance = zeros([lastUsed 1]);
            ratio = zeros([lastUsed 1]);
            rot = zeros([lastUsed 1]);
            altRatio = zeros([lastUsed 1]);
            for k=1:lastUsed
                distance(k) = obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(3);
                ratio(k) = obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(5);
                rot(k) = obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(6);
                xpos(k) = obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(1);
                ypos(k) = obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(2);
                ID(k) = obj.locData.SE.sites(k).ID;
                if p.oldversion
                    altRatio(k) = obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(5);
                else
                    altRatio(k) = obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(14);
                end
            end
            indZero = ratio == 0;
            logRatio = log10(ratio);
            logRatio(indZero) = log10(altRatio(indZero));
            rmLogInf = ratio == 0 & altRatio==0;
            rmRot = abs(rot) >= p.rotLimit;
            rmLogRatio = logRatio<=log10(p.relativeAmpLimitLower)|logRatio>=log10(p.relativeAmpLimit);
            %% not yet in GUI
            if p.rmBound
                rmXpos = abs(xpos) >= 60;
                rmYpos = ypos >= 120|ypos <= 30;
                rmDis = distance <= -200;
                rm = rmRot|rmLogRatio|rmLogInf|rmXpos'|rmYpos'|rmDis;
            else
                rm = rmRot|rmLogRatio|rmLogInf;
            end
            site2rm = find(rm);
            
            if p.viz
                plotSElink(ax,distance,logRatio,ID,obj.locData.SE,' ob')
                hold(ax, 'on')
                plot(ax,distance(rmRot),logRatio(rmRot),' or')
                plot(ax,distance(rmLogRatio),logRatio(rmLogRatio),' xg')
                hold(ax, 'off')
                ax.XLabel.String = 'Distance (nm)';
                ax.YLabel.String = 'Log10(Ratio)';

                plotSElink(ax2,distance,rot,ID,obj.locData.SE,' ob')
                ax2.XLabel.String = 'Distance (nm)';
                ax2.YLabel.String = 'Rotation (degree)';
            end

            %% Disable bad sites
            if p.applyFilter
                for k=site2rm'
                    obj.locData.SE.sites(k).annotation.use = 0;
                end

                %% Sort sites according to 'use'
                fdcal=figure(519);
                dcal=plugin('ROIManager','Analyze','SortROIs',fdcal,obj.P);
                dcal.attachLocData(obj.SE.locData);
                dcal.makeGui;
                pp=dcal.getGuiParameters;
                pp.sortprop1.Value=2;                               % select use
                pp.sortedit1 = 'annotation.use';                    % select use
                pp.direction1.Value=2;                              % decend
                dcal.setGuiParameters(pp);
                dcal.processgo;
                pause(10e-1000);
            end
            
            %% Update position
            if p.updatePos
                lastUsed = find(getFieldAsVector(obj.locData.SE.sites,'annotation.use'),1, 'last');
                for k = 1:lastUsed
                    siteRotPos = obj.locData.SE.sites(k).annotation.rotationpos.pos;
                    [obj.locData.SE.sites(k).annotation.rotationpos.pos(:,1),obj.locData.SE.sites(k).annotation.rotationpos.pos(:,2)] = rotcoord(siteRotPos(:,1), siteRotPos(:,2), obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(6)*pi/180);
                end

                for k=1:lastUsed
                    dx = 0 - obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(1);
                    dy = 120 - obj.locData.SE.sites(k).evaluation.CME2Dfitting2.parameters(2);
                    [dxrot,dyrot] = rotcoord(dx,dy,-obj.locData.SE.sites(k).annotation.rotationpos.angle*pi/180);
                    obj.locData.SE.sites(k).pos = obj.locData.SE.sites(k).pos - [dxrot dyrot 0];
                end
            end
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function pard=guidef(obj)

pard.applyFilter.object=struct('String','Apply filter','Style','checkbox', 'Value', 0);
pard.applyFilter.position=[2,1];
pard.applyFilter.Width=1;

pard.updatePos.object=struct('String','Update pos','Style','checkbox', 'Value', 0);
pard.updatePos.position=[2,2];
pard.updatePos.Width=1;

pard.viz.object = struct('String','Visualization', 'Style', 'checkbox', 'Value',1);
pard.viz.position=[3,1];
pard.viz.Width=1;

pard.t1.object = struct('String','Max relative amp.:', 'Style', 'text');
pard.t1.position=[4,1];
pard.t1.Width=1;

pard.relativeAmpLimit.object = struct('String','3', 'Style', 'edit');
pard.relativeAmpLimit.position=[4,2];
pard.relativeAmpLimit.Width=0.5;

pard.relativeAmpLimitLower.object = struct('String','0.1', 'Style', 'edit');
pard.relativeAmpLimitLower.position=[4,2.5];
pard.relativeAmpLimitLower.Width=0.5;

pard.t2.object = struct('String','Max rotating angle:', 'Style', 'text');
pard.t2.position=[5,1];
pard.t2.Width=1;

pard.rotLimit.object = struct('String','40', 'Style', 'edit');
pard.rotLimit.position=[5,2];
pard.rotLimit.Width=0.5;

end
