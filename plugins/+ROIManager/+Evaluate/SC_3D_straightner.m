classdef SC_3D_straightner<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
        %roicoordinates
        roihandles={[]};
        roisize
        roisizey
%         currentroi=0;
        %axis
        %images
%         roihandle
    end
    methods
        function obj=SC_3D_straightner(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            out=[];
            display=obj.display;

            polylinesource=obj.getSingleGuiParameter('polylinesource').selection;
            expandPL=obj.getSingleGuiParameter('polylineexpandT');
            expandby=obj.getSingleGuiParameter('polylineexpand');
            
            saveLocsforsiteT=p.saveLocsforsiteT;

            axeslayer=str2num(obj.getSingleGuiParameter('axesLayer').selection);
            calangle=obj.getSingleGuiParameter('axesangle');
            angleplot=obj.getSingleGuiParameter('angleplot');
            
            cutoffdistance=200;
            markerlayer=str2num(obj.getSingleGuiParameter('COLayer').selection);
            binwidth=5;
            centre=~obj.getSingleGuiParameter('centerCO');

            createsubseg=obj.getSingleGuiParameter('createsubsegT');
            subseglength=obj.getSingleGuiParameter('createsubseg');
            out.createsubseg=createsubseg;
            out.subseglength=subseglength;

            if ~isfield(obj.site.evaluation,polylinesource)
                error('Please choose a polyline source that is available. %s has not been evaulated',polylinesource)
           
            end
            layerson = find(obj.locData.getPar('sr_layerson'));             % check the layers used
            visitFlag = false;                                      % a flag for the conditional loop for the first round
            fieldQ = {'locprecnm','locprecznm','PSFxnm','xnm','znm','ynm','frame','xnmrot','ynmrot','phot1','phot2','filenumber','channel'};    % fields will be used %%%%%%%%%%%%%%%%%%% need to add more for displayer
            
            %layerson=layerson(layerson~=cosa1layer);
            fieldQNLayer = [fieldQ 'layer'];
            for k = layerson                                        % go through layers, collect all filtered locs -----from Yu-Le's LocMoFitGUI.m
                if ~visitFlag
                    locs=obj.getLocs(fieldQ,'layer',k,'size','freeroi');
                    locs.layer = ones(size(locs.(fieldQ{1})))*k;
                    fieldExact = fieldnames(locs);
                    lRm = ~ismember(fieldExact, fieldQNLayer);
                    locs = rmfield(locs,fieldExact(lRm));
                    visitFlag=true;
                else
                    locsNewLayer = obj.getLocs(fieldQ,'layer',k,'size','freeroi');
                    for l = length(fieldQNLayer):-1:1
                        if l == length(fieldQNLayer)
                            locsNewLayer.layer = ones(size(locsNewLayer.(fieldQNLayer{l-1})))*k;
                        end
                        if isfield(locs, fieldQNLayer{l})
                            locs.(fieldQNLayer{l})=[locs.(fieldQNLayer{l});locsNewLayer.(fieldQNLayer{l})];
                        end
                    end
                end
            end
            

            
            %%Missing an option to just straighten based on the
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%polyline

            %derivedparameters=obj.site.evaluation.(polylinesource).GuiParameters.fitter.getDerivedPars;
            %descriptors=derivedparameters{1,1};%derivePars(obj,obj.site.evaluation.(polylinesource).allParsArg,polylinesource);%
            if contains(polylinesource,'LocMoFitGUI')
                derivedparameters=obj.site.evaluation.(polylinesource).fitInfo.derivedPars;
                descriptors=derivedparameters{1,1};
                if createsubseg
                    npts=round(arclength(descriptors.ctrlpoints.ctrlX,descriptors.ctrlpoints.ctrlY,descriptors.ctrlpoints.ctrlZ,'linear')/subseglength);
                    [descriptors.pt, descriptors.dudt, descriptors.funthd]=interparc(npts,descriptors.ctrlpoints.ctrlX,descriptors.ctrlpoints.ctrlY,descriptors.ctrlpoints.ctrlZ,'pchip');% to lin 230404
                    [L,R,K] = curvature(descriptors.pt);
                    descriptors.curvature=1./R;
                end
            elseif contains(polylinesource,'plot3D_annotate')
                descriptors=struct();
                rotated=PolylineToInit(obj.site.evaluation.(polylinesource).roicoordinates.Position,obj.site.pos,obj.site.annotation.rotationpos.angle/180*pi,size(obj.site.evaluation.(polylinesource).roicoordinates.Position,1));
                rotated=[rotated.xnmR rotated.ynmR rotated.znmR];
                descriptors.ctrlpoints.ctrlX=rotated(:,1);
                descriptors.ctrlpoints.ctrlY=rotated(:,2);
                descriptors.ctrlpoints.ctrlZ=rotated(:,3);
                if createsubseg
                    npts=round(arclength(descriptors.ctrlpoints.ctrlX,descriptors.ctrlpoints.ctrlY,descriptors.ctrlpoints.ctrlZ,'linear')/subseglength);
                    [descriptors.pt, descriptors.dudt, descriptors.funthd]=interparc(npts,descriptors.ctrlpoints.ctrlX,descriptors.ctrlpoints.ctrlY,descriptors.ctrlpoints.ctrlZ,'pchip'); %pchip
                    [L,R,K] = curvature(descriptors.pt);
                    descriptors.curvature=1./R;
                else
                    descriptors.pt=rotated;
                    [L,R,K] = curvature(descriptors.pt);
                    descriptors.curvature=1./R;
                end
            else
                error('It is not possible to use %s as source of polyline',polylinesource)
            end
%             [ctrlX,ctrlY,ctrlZ]=getCtrlPointsXYZloc(obj.site.evaluation.(polylinesource).allParsArg);
            
%             minD = 50;%*0.75;%
%             arcLen=arclength(ctrlX, ctrlY, ctrlZ,'linear');
%             samplingFactor = round(arcLen/minD);
%             [pt,dudt]= interparc(round(samplingFactor/2), double(ctrlX), double(ctrlY), double(ctrlZ),'pchip');%100
%             descriptors.pt=pt;
%             [L,R,K] = curvature(pt);
%             descriptors.curvature = 1./R;
%             descriptors.ctrlpoints.ctrlX = ctrlX;
%             descriptors.ctrlpoints.ctrlY = ctrlY;
%             descriptors.ctrlpoints.ctrlZ = ctrlZ;
            
            %If there is no user specified polygon mask in annotation
            %tab, create a polygon mask from a polyline
            if (~isvalid(obj.locData.SE.processors.preview.hlines.line3)||isempty(obj.locData.SE.processors.preview.hlines.line3)) || expandPL
                polygon=polybuffer(descriptors.pt(:,1:2),'lines',expandby); 
                [in,on] = inpolygon(locs.xnmrot,locs.ynmrot, polygon.Vertices(:,1),polygon.Vertices(:,2));
                locsOri=locs;
                locs=getFields(locs,in);
            end

            %locs=obj.site.evaluation.LocMoFitGUI_2.GuiParameters.fitter.locsHandler(locs,obj.site.evaluation.LocMoFitGUI_2.GuiParameters.fitter.exportPars(1,'lPar'))
            locs.xnmA=locs.xnmrot;
            locs.ynmA=locs.ynmrot;
            locs.znmA=locs.znm;

            warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale'); %descriptors.pt 
            %Fix version 2 and use it for better straigtening - shoud go as
            %-interparc and then use every other point as transformation
            %centre
            %Also need to see if 40 gausian spread is enough if fittd with
            %axes[descriptors.ctrlpoints.ctrlX descriptors.ctrlpoints.ctrlY descriptors.ctrlpoints.ctrlZ]
%             straightner=straigthen(locs,descriptors.pt,createsubseg,subseglength, axeslayer, calangle);
            straightner=straigthen_v2(locs,descriptors.pt,axeslayer, calangle);
            warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale'); 

            out.straightner=straightner;
            out.cspline=descriptors.pt;
            out.curvature=descriptors.curvature;
            out.ctrlpoints=descriptors.ctrlpoints;
            out.polylinesource=polylinesource;
            
            %Assosiate CO marker to the curvature
            %Add clustering for MSH-5 and other SMLM, and then cal mean of
            %cluster or assosiate all and identify peaks-v1
%             cosa1file=2;
            
            warning('off', 'stats:pdist2:DataConversion');
            if centre==1
                if sum(locs.layer==markerlayer)>0 && sum(locs.layer==markerlayer)<10 
                    markerpoints=[locs.xnmA(locs.layer==markerlayer) locs.ynmA(locs.layer==markerlayer) locs.znmA(locs.layer==markerlayer)];
                    distances=pdist2(markerpoints,descriptors.pt,'euclidean');
                    nmarkerpt=sum(locs.layer==markerlayer);
                    [minimumd,indices] = min(distances,[],2);
                    [mofm,idcentre]=min(minimumd); %centre adound closest, in future for SMLM I need to add cluster size
                    idcentre=indices(idcentre);
                elseif sum(locs.layer==markerlayer)==0
                    %no centering but replicate output
                    mofm=0;
                    idcentre=0;
                    nmarkerpt=0;
                else
                    markerpoints=[locs.xnmA(locs.layer==markerlayer) locs.ynmA(locs.layer==markerlayer) locs.znmA(locs.layer==markerlayer)];
                    idx=dbscan(markerpoints,100,50);
                    [GC,GR] = groupcounts(idx);
                    [maxV,pickcluster]=max(GC);
                    meanmarkerpoints=mean(markerpoints(idx==GR(pickcluster),:),1);
                    distances=pdist2(meanmarkerpoints,descriptors.pt,'euclidean');
                    nmarkerpt=size(markerpoints,2);
                    [minimumd,indices] = min(distances,[],2);
                    [mofm,idcentre]=min(minimumd);
                    idcentre=indices(idcentre);

                end
                out.center.dist=mofm;
                out.center.centerpoint=idcentre;
                out.center.nmarkerpt=nmarkerpt;
                
            end
            warning('on', 'stats:pdist2:DataConversion');
            %Calculate angle numerical gradient
            if calangle==1
                dAXangle=gradient(straightner.axesangles);
                fordistance=straightner.locs.znmS(straightner.straight & locs.layer==axeslayer);
                bins=linspace(min(fordistance),max(fordistance),round((max(fordistance)-min(fordistance))/binwidth));
                histD=histogram(fordistance,bins);
                binsC=arrayfun(@(i)(bins(i)+bins(i+1))/2,1:(size(bins,2)-1));
                out.DGcoefX=myDoubleGaussian_simetric(binsC,histD.Values);

                fordistance=straightner.locs.ynmS(straightner.straight & locs.layer==axeslayer);
                bins=linspace(min(fordistance),max(fordistance),round((max(fordistance)-min(fordistance))/binwidth));
                histD=histogram(fordistance,bins);
                binsC=arrayfun(@(i)(bins(i)+bins(i+1))/2,1:(size(bins,2)-1));
                out.DGcoefZ=myDoubleGaussian_simetric(binsC,histD.Values);
                out.axesanglesQ1=restricttoquadrant1(straightner.axesangles);
                out.daxesanglesQ1=gradient(out.axesanglesQ1);
            end
            

            if display
                figure('Name','3dpreview - localizations and polyline');
                hold on
%                 scatter3(descriptors.ctrlpoints.ctrlX,descriptors.ctrlpoints.ctrlY,descriptors.ctrlpoints.ctrlZ,'green','filled')
                scatter3(descriptors.pt(:,1),descriptors.pt(:,2),descriptors.pt(:,3),[30],'green','filled')
                scatter3(straightner.correctedPt(:,1),straightner.correctedPt(:,2),straightner.correctedPt(:,3),[30],'magenta','filled')
                scatter3(locs.xnmA, locs.ynmA, locs.znmA,[5],locs.layer,'filled'); % 
                colormap("lines")
                daspect([1 1 1])
                hold off

                figure('Name','3dpreview - localizations and polyline');
                hold on
%                 scatter3(descriptors.ctrlpoints.ctrlX,descriptors.ctrlpoints.ctrlY,descriptors.ctrlpoints.ctrlZ,'green','filled')
                scatter3(descriptors.pt(:,1),descriptors.pt(:,2),descriptors.pt(:,3),[30],straightner.curvature_2,'filled')
                colormap("hot")
                colorbar
                daspect([1 1 1])
                hold off

                figure('Name','3dpreview - straightened localizations');
                hold on
                scatter3(straightner.locs.xnmS(straightner.straight), straightner.locs.ynmS(straightner.straight),straightner.locs.znmS(straightner.straight),[5],straightner.locs.layer(straightner.straight),'filled'); % 
                colormap("lines")
                daspect([1 1 1])
                hold off

%                 figure('Name','3dpreview - straightened localizations');
%                 hold on
%                 scatter3(straightner.locs.xnmS, straightner.locs.ynmS,straightner.locs.znmS,[5],straightner.locs.layer,'filled'); % 
%                 colormap("lines")
%                 daspect([1 1 1])
%                 hold off

%                 figure('Name','curvature-not centered yet');
%                 hold on
%                 plot(straightner.curvature); % 
%                 hold off
            end  
            
           if saveLocsforsiteT
               
               locsS=getFields(straightner.locs,straightner.straight);
               dateC=string(datetime('now', 'Format','yyMMdd'));
               %filepath=obj.locData.files.file(obj.site.info.filenumber).name;
               path=regexp(obj.locData.files.file(obj.site.info.filenumber).name,'\\','split');
               newpath=[strjoin([path(1:end-1)],"\\") '\\straightenedSCs'];
               if ~exist(newpath, 'dir')
                   mkdir(newpath);
               end

%                
%                modelPts.xnm=descriptors.pt(:,1)+abs(min(descriptors.pt(:,1)));
%                modelPts.ynm=descriptors.pt(:,2)+abs(min(descriptors.pt(:,2)));
%                modelPts.znm=descriptors.pt(:,3)+abs(min(descriptors.pt(:,3)));
               modelPts.curvature=straightner.curvature_2;
               modelPts.torsion=straightner.torsion_2;
%                locsS.xnmS=locsS.xnmS+abs(min(locS.xnmS));
%                locsS.ynmS=locsS.xnmS+abs(min(locS.ynmS));
%                locsS.znmS=locsS.xnmS+abs(min(locS.znmS));
               if calangle==1
                   dangles=straightner.axesangles;
                   modelPts.twist=dangles;
                   CSV_2=[newpath '\\' path{end}(1:end-4) obj.site.name(3:end-1) '_' dateC{1} '_LocMoFitModel.csv'];
                   CSV=[newpath '\\' path{end}(1:end-4) obj.site.name(3:end-1) '_' dateC{1} '_straightenedSC.csv'];
                   else
                   
                   CSV_2=[newpath '\\' path{end}(1:end-4) obj.site.name(3:end-1) '_' dateC{1} '_notwistcorr_LocMoFitModel.csv'];
                   CSV=[newpath '\\' path{end}(1:end-4) obj.site.name(3:end-1) '_' dateC{1} '_notwistcorr_straightenedSC.csv'];
               
               end
               final_2=struct2table(modelPts);
               
               final=struct2table(locsS);
               writetable(final,CSV,'WriteRowNames',true);
               writetable(final_2,CSV_2,'WriteRowNames',true);
           end

           if angleplot==1 && calangle==1 && display
%                forplot=straightner.axesangles;
%                changes=abs(diff(straightner.axesangles));
%                for s=1:size(changes,1)
%                    if changes(s)>2.7
%                        forplot((s+1):end)=forplot((s+1):end)
%                    
%                    end end
               figure('Name','Angle between axes');
               hold on
               plot(rad2deg(restricttoquadrant1(straightner.axesangles))); % 
%                ylim([0 180])
               hold off

           end

%             make3Daxis(obj) %make axis,
%             plot3Dviews(obj,p);%plot image,
%             if ~isfield(obj.site.annotation,'polarangle')
%                 obj.site.annotation.polarangle=0;
%             end
%             obj.currentroi=0;
%             if isfield(obj.site.evaluation,obj.name) && isfield(obj.site.evaluation.(obj.name),'roicoordinates')
%                 obj.roicoordinates=obj.site.evaluation.(obj.name).roicoordinates;
%                 obj.plotrois;
%                 out.GuiParameters=obj.site.evaluation.(obj.name).GuiParameters; %ICfix_210515 --- CHECK 
%                 out.roicoordinates=obj.site.evaluation.(obj.name).roicoordinates; %ICfix_210515 --- CHECK 
%                 %out=obj.site.evaluation.(obj.name);  - use this above
%                 %instead of two lines
%             else
%                 obj.roicoordinates=[];
%                
%             end
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
%         
    end

end

% function ind=findextra(A,B)
% xA=A(:,1);
% xB=B(:,1);
% mind=min((xA-xB').^2,[],2);
% [~,ind]=max(mind);
% end
% 
% function [roipospix,roipos]=getpixcoord(obj,roiposh,roiposold)
% roipospix=roiposold;
%  if mean(roiposh(:,1))>0
%     roipos=2; %right
%     roiposh(:,1)=max(roiposh(:,1),0);
%     roipospix(:,1)=roiposh(:,1)-obj.roisize;
%     roipospix(:,3)=roiposh(:,2); %move to left
% %     roipospix(:,2)=0;
% else
%     roipos=1;
%     roiposh(:,1)=min(roiposh(:,1),0);
%     roipospix(:,1:2)=roiposh(:,1:2);
% %     roipospix(:,3)=0;
% end
% end
% 
% function roiposnm=coord3Dpix2nm(roipos,pos,rotationanglez,polarangle,winsize)
%     xp=roipos(:,1)+winsize/2;
%     yp=roipos(:,2);zp=roipos(:,3);
%     [y2,zi]=rotcoorddeg(yp,zp,polarangle);
%     [xi,yi]=rotcoorddeg(xp,y2,rotationanglez);
%     roiposnm(:,1)=xi+pos(1); 
%     roiposnm(:,2)=-yi+pos(2); %%%%220329 -yi
%     roiposnm(:,3)=zi+pos(3);
% end
% function roipos=coord3Dnm2pix(roiposnm,pos,rotationanglez,polarangle,winsize)
%     xi=roiposnm(:,1)-pos(1); 
%     yi=-(roiposnm(:,2)-pos(2)); %%%%220329 -(...)
%     zi=roiposnm(:,3)-pos(3);
%     [xp,y2]=rotcoorddeg(xi,yi,-rotationanglez); %-
%     [yp,zp]=rotcoorddeg(y2,zi,-polarangle);
%     roipos(:,1)=xp-winsize/2;
%     roipos(:,2)=yp;
%     roipos(:,3)=zp;
% end

function pard=guidef(obj)
% pard.addroi.object=struct('Style','pushbutton','String','add','Callback',@obj.addroi_callback);
% pard.addroi.position=[1,1];
% pard.addroi.Width=1;
% 
% pard.roiform.object=struct('Style','popupmenu','String',{{'polyline','line','free','polygon'}});
% pard.roiform.position=[1,2];
% pard.roiform.Width=2;
% 
% 
% pard.deleterois.object=struct('Style','pushbutton','String','delete all rois','Callback',@obj.deleteroi_callback);
% pard.deleterois.position=[2,1];
% pard.deleterois.Width=2;

pard.polylinesourceT.object=struct('String','3D polyline source:','Style','text');
pard.polylinesourceT.position=[3.25,1];
pard.polylinesourceT.Width=2;

availableProcessors={cellfun(@(x) x.name,obj.locData.SE.processors.eval.processors,'UniformOutput',false)};
in=regexp([availableProcessors{1,1:end}],'(LocMoFitGUI\w*|plot3D_annotate\w*)','match');
pard.polylinesource.object=struct('Style','popupmenu','String',{[in{~cellfun('isempty',in)}]});
pard.polylinesource.position=[3,3];
pard.polylinesource.Width=2;
pard.polylinesource.TooltipString='Specified previously run evaluation that will serve as a  source for the polyline. Avaliable LocMoFit and plot3Dannotation(still in progress).';

pard.polylineexpandT.object=struct('String','Expand polyline to mask by (nm):','Style','checkbox','Value',0);
pard.polylineexpandT.position=[4,1];
pard.polylineexpandT.Width=3;
pard.polylineexpandT.TooltipString='In case there is no user specified mask within annotation tab this will be run irrespective of the status here.';

pard.polylineexpand.object=struct('String','350','Style','edit');
pard.polylineexpand.position=[4,4.5];
pard.polylineexpand.Width=0.5;
pard.polylineexpand.TooltipString='Expand around polyline by nm. Useful when ROI contains many localizations that do not belong to the SC stretch that will be straigthened';


pard.saveLocsforsiteT.object=struct('String','Save straightened localizations for site.','Style','checkbox','Value',0);
pard.saveLocsforsiteT.position=[5,1];
pard.saveLocsforsiteT.Width=5;
pard.saveLocsforsiteT.TooltipString='Straighetend localizations will be saved as /path/to/current/file/straightenedSCs/filename_siteID.csv overwriting the previous version.';

pard.axesangle.object=struct('String','Calculate angle.','Style','checkbox','Value',0);
pard.axesangle.position=[6,1];
pard.axesangle.Width=2;
pard.axesangle.TooltipString='If checked, angle between axes will be calculated along the length of the cspline/midline by fitting a line. Values will be saved within the evaluation output.';

pard.anglemethod.object=struct('Style','popupmenu','String',{{'line','polaranglebinning','covariance'}});%
pard.anglemethod.position=[6,3];
pard.anglemethod.Width=2;
pard.anglemethod.TooltipString='Options are not available. Only line is implemented';

pard.angleplot.object=struct('String','Plot angle.','Style','checkbox','Value',0);
pard.angleplot.position=[8,1];
pard.angleplot.Width=2;
pard.angleplot.TooltipString='If checked when Display is allowed, the change of angle between the axes will be ploted.';


pard.axesLayerT.object=struct('String','Axes Layer:','Style','text');
pard.axesLayerT.position=[7.25,1];
pard.axesLayerT.Width=2;

pard.axesLayer.object=struct('Style','popupmenu','String',{find(obj.locData.getPar('sr_layerson'))});
pard.axesLayer.position=[7,3];
pard.axesLayer.Width=1;
pard.axesLayer.TooltipString='Choose layer that contains chromosome axes localizations.';

pard.createsubsegT.object=struct('String','Straighten subsegments of (nm):','Style','checkbox','Value',1);
pard.createsubsegT.position=[9,1];
pard.createsubsegT.Width=3;
pard.createsubsegT.TooltipString='If checked each segment of the polyline will be devided into subsegments of specified length which will be separately straightened.';


pard.createsubseg.object=struct('String','5','Style','edit');
pard.createsubseg.position=[9,4.5];
pard.createsubseg.Width=0.5;
pard.createsubseg.TooltipString='Specify the length of subsegments.';

pard.COLayerT.object=struct('String','CO marker Layer:','Style','text');
pard.COLayerT.position=[1,1];
pard.COLayerT.Width=2;

pard.COLayer.object=struct('Style','popupmenu','String',{find(obj.locData.getPar('sr_layerson'))});
pard.COLayer.position=[1,2.5];
pard.COLayer.Width=1;
pard.COLayer.TooltipString='Choose layer that contains crossover marker localizations.';

pard.centerCO.object=struct('String','Do not center.','Style','checkbox','Value',0);
pard.centerCO.position=[2,1];
pard.centerCO.Width=2;
pard.centerCO.TooltipString='If checked there will be NO centering around the largest cluster of CO marker localizations.';


% pard.equaldistance.object=struct('Style','checkbox','String','equal distances');
% pard.equaldistance.position=[2,1];
% pard.equaldistance.Width=2;

% pard.fitpeaks.object=struct('Style','pushbutton','String','fit peaks','Callback',@obj.fit_callback);
% pard.fitpeaks.position=[1,1];
% pard.fitpeaks.Width=2;



% pard.dxt.Width=3;
% pard.inputParameters={'numberOfLayers','sr_layerson','se_cellfov','se_sitefov','se_siteroi','se_sitepixelsize'};
pard.plugininfo.type='ROI_Evaluate';
pard.plugininfo.description='';
end

function plot3Dstraigtened(obj,p,locsS)
if nargin<2
    p=[];
    in=true;
elseif nargin<3
    in=true;
end
site=obj.site;
angle=pos2angle(site.annotation.rotationpos.pos);
p.rotationanglez=angle;
if isfield(site.annotation,'polarangle')
    p.polarangle=site.annotation.polarangle;
else
    p.polarangle=0;
end
% p1.sr_size(3)=obj.locData.getPar('se_dz')/2;


% imz=make3Dimages(obj,p);
% end         
                
% function imagez=make3Dimages(obj,p)                
% imagez=[];
% fileind=obj.locData.SE.indexFromID(obj.locData.SE.files,obj.site.info.filenumber);
numlayers=obj.locData.getPar('numberOfLayers');
allfields=[renderSMAP drawerSMAP displayerSMAP];
players=obj.getLayerParameters(1:numlayers,allfields);
if ~iscell(players)
   players={players};
end

if obj.locData.getPar('se_imaxcheck_site')
    p.imaxtoggle=false;
    imax=obj.locData.getPar('se_imax_site');

    for k=1:length(imax)
        players{k}.imax_min=imax(k);
    end
    for k=length(imax)+1:numlayers
        players{k}.imax_min=imax(1);
    end
end
for k=1:numlayers
    prz=copyfields(players{k},p);
    
    prz.sr_pos=site.pos;  %later move all of this out of the loop
    prz.sr_size=ones(2,1)*obj.locData.getPar('se_sitefov')/2;
    prz.sr_pixrec=obj.locData.getPar('se_sitepixelsize');
    prz.sr_sizeRecPix=round((prz.sr_size*2)/prz.sr_pixrec);
    prz.sr_axes=-1;
    prz.normalizeFoV=prz.sr_sizeRecPix(1)/obj.locData.getPar('se_sitefov')*obj.locData.getPar('se_siteroi')/2;
                 
%     posh=[prz.sr_pos(1) prz.sr_pos(2) prz.sr_size(1)*2 prz.sr_size(2)*2];
    locz=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','PSFxnm','intensity_render','phot','numberInGroup',prz.renderfield.selection},'layer',k,'position',obj.site);
    if length(prz.sr_pos)<3
        prz.sr_pos(3)=0;
    end
    xi=locz.xnm-site.pos(1);
    yi=locz.ynm-site.pos(2);
    zi=locz.znm-site.pos(3);
    [x2,y2]=rotcoorddeg(xi,yi,prz.rotationanglez);
    [y3,z3]=rotcoorddeg(y2,zi,prz.polarangle);
    locz.sx=locz.locprecnm;locz.sy=locz.locprecznm;
    prz.normalizeFoV=[];
    prz.sr_pos=[0,0,0];
    prz.sr_size(3)=obj.locData.getPar('se_dz')/2;
    prz.sr_size(2)=prz.sr_size(3);
    locz.x=x2;locz.y=z3;
    rawimagez=renderSMAP(locz,prz,k);
    locz.x=x2;locz.y=y3;
    rawimagezxy=renderSMAP(locz,prz,k);
    layersz(k).images.finalImages=drawerSMAP(rawimagez,prz);
    layerszxy(k).images.finalImages=drawerSMAP(rawimagezxy,prz);
    imax(k)=max(layersz(k).images.finalImages.imax,layerszxy(k).images.finalImages.imax);                   
end
imz{1}=displayerSMAP(layersz,prz);
imz{2}=displayerSMAP(layerszxy,prz);  

imzc.image=horzcat(imz{2}.image(end:-1:1,:,:), imz{1}.image); 
imzc.rangey=imz{1}.rangey;
imzc.rangex=[ imz{2}.rangex(1)-imz{2}.rangex(2) imz{1}.rangex(2)-imz{1}.rangex(1)];

displayimage(imzc,obj.axis);
obj.roisize=imzc.rangex(2)*1000;
obj.roisizey=imzc.rangey(2)*1000;
fl='%2.0f';
title(obj.axis,['\theta=' num2str(p.polarangle,fl) '\circ, \rho=' num2str(p.rotationanglez,fl) '\circ, z= ' num2str(site.pos(3),fl) ' nm'])
obj.images=imzc; 
end



function info(a,b)
text='To rotate in x-y plane: click in the left top view image in the direction you want to point left. \nTo change the z-position, click on the right image, left part, this position will be centered. \nTo change the polar angle, click in the right image, right part. The closer you are to the center, the smaller the change.';
msgbox(sprintf(text));
end

function locI=getFields(loc,indices, FieldList)
    if nargin == 2
        FieldList = fieldnames(loc);
    end 
    for iField = 1:numel(FieldList)
        Field=FieldList{iField};
        locI.(Field) = loc.(Field)(indices);
    end
end
function out=myDoubleGaussian_simetric_old(x,y,coef0)
    %Coef order: C, sigma, centre, shift of each peak
    if nargin<3
      coef0 = [10;30;0;60];
    end
    modelfun=@(coef,x)(coef(1).*(1/(coef(2).*sqrt(2*pi)))*(exp(-0.5.*(((x-(coef(3)-coef(4))).^2)./(coef(2)^2)))+exp(-0.5.*(((x-(coef(3)+coef(4))).^2)./(coef(2).^2)))));
    opts = statset('nlinfit');
    opts.TolFun=1e-10;%,'DerivStep',1
%     opts.RobustWgtFun = 'andrews';
    mdl = fitnlm(x,y,modelfun,coef0,'Options',opts);
    out=mdl.Coefficients.Estimate;
    %[beta,R,J,CovB,MSE,ErrorModelInfo]
    
end

function out=myDoubleGaussian_simetric(x,y,coef0)
    %Coef order: amplitude, sigma, centre, peak shift __III__III__
    if nargin<3
      coef0 = [10;30;0;60];
    end
    modelfun=@(coef,x)(coef(1).*(1/(coef(2).*sqrt(2*pi)))*(exp(-0.5.*(((x-(coef(3)-coef(4))).^2)./(coef(2)^2)))+exp(-0.5.*(((x-(coef(3)+coef(4))).^2)./(coef(2).^2)))));
    obj_fun = @(coefs) sum((modelfun(coefs, x)-y).^2)/size(x,2);%norm(modelfun(coefs, x)-y) currently using MSE
    opts = optimset('TolFun',1e-10,'TolX',1e-10);%,'Display','iter','PlotFcns',@optimplotfval
    [out.estimate,out.fval,out.exitflag,out.output] = fminsearchbnd(obj_fun, coef0,[10; 10; -50; 0],[Inf; Inf; 50; 200],opts);
    %[beta,R,J,CovB,MSE,ErrorModelInfo]
    
end

function out=myDoubleGaussian_asimetric(x,y,coef0)
    %Coef order: amplitude1, amplitude2, sigma, centre, peak shift __iii__III__
    if nargin<3
      coef0 = [10;10;30;0;60];
    end
    modelfun=@(coef,x)(coef(1).*(1/(coef(3).*sqrt(2*pi)))*(exp(-0.5.*(((x-(coef(4)-coef(5))).^2)./(coef(3)^2))))+(coef(2).*(1/(coef(3).*sqrt(2*pi)))*exp(-0.5.*(((x-(coef(4)+coef(5))).^2)./(coef(3).^2)))));
    obj_fun = @(coefs) sum((modelfun(coefs, x)-y).^2)/size(x,2);%norm(modelfun(coefs, x)-y) currently using MSE
    opts = optimset('TolFun',1e-10,'TolX',1e-10);
    [out.estimate,out.fval,out.exitflag,out.output] = fminsearchbnd(obj_fun, coef0,[10; 10; 10; -50; 60],[Inf; Inf; Inf; 50; 80],opts);
    %[beta,R,J,CovB,MSE,ErrorModelInfo]
    
end
function out=restricttoquadrant1(radangles)
    out=radangles;
    out(out>pi/2)=out(out>pi/2)-pi;
    out(out<0)=abs(out(out<0));
end
function out=axisDoubleGaussian_asimetric(x,y,coef0)
    %Coef order: amplitude1, amplitude2, sigma, centre, peak shift __iii__III__
    if nargin<3
      coef0 = [10;10;30;0;60];
    end
    modelfun=@(coef,x)(coef(1).*(1/(coef(3).*sqrt(2*pi)))*(exp(-0.5.*(((x-(coef(4)-coef(5))).^2)./(coef(3)^2))))+(coef(2).*(1/(coef(3).*sqrt(2*pi)))*exp(-0.5.*(((x-(coef(4)+coef(5))).^2)./(coef(3).^2)))));
    obj_fun = @(coefs) sum((modelfun(coefs, x)-y).^2)/size(x,2);%norm(modelfun(coefs, x)-y) currently using MSE
    opts = optimset('TolFun',1e-10,'TolX',1e-10);
    [out.estimate,out.fval,out.exitflag,out.output] = fminsearchbnd(obj_fun, coef0,[10; 10; 10; -40; 40],[Inf; Inf; Inf; 40; 100],opts); %50 and 80 originally
    %[beta,R,J,CovB,MSE,ErrorModelInfo]
    
end

function out=axisSingleGaussian(x,y,coef0)
    %Coef order: amplitude, sigma, centre
    if nargin<3
      coef0 = [10;10;10];
    end
    modelfun=@(coef,x)(coef(1)*(1/(coef(2)*sqrt(2*pi)))*(exp(-0.5*((x-coef(3)).^2)./(coef(2)^2))));
    obj_fun = @(coefs) sum((modelfun(coefs, x)-y).^2)/size(x,2);%norm(modelfun(coefs, x)-y) currently using MSE
    opts = optimset('TolFun',1e-10,'TolX',1e-10);
    [out.estimate,out.fval,out.exitflag,out.output] = fminsearchbnd(obj_fun, coef0,[10; 30; -60],[Inf; Inf; 60],opts);
    
end

function out=mySingleGaussian(x,y,coef0)
    %Coef order: amplitude, sigma, centre
    if nargin<3
      coef0 = [10;10;10];
    end
    modelfun=@(coef,x)(coef(1)*(1/(coef(2)*sqrt(2*pi)))*(exp(-0.5*((x-coef(3)).^2)./(coef(2)^2))));
    obj_fun = @(coefs) sum((modelfun(coefs, x)-y).^2)/size(x,2);%norm(modelfun(coefs, x)-y) currently using MSE
    opts = optimset('TolFun',1e-10,'TolX',1e-10);
    [out.estimate,out.fval,out.exitflag,out.output] = fminsearchbnd(obj_fun, coef0,[10; 10; -50],[Inf; Inf; 50],opts);
    
end

% function out=axAngle_by_polaranglebinning(coord_2D)
% end
% function out=axAngle_by_GMCline(coord_2D)
% end
% function out=axAngle_by_covariance(coord_2D)
% end
function out=straigthen_v2(locs,polyline, anglelayer, calangle, selR)
% Straighten the SC along the 3D polyline 
    %links to the used functions
    %(https://uk.mathworks.com/matlabcentral/answers/473087-how-do-you-calculate-a-trajectory-through-a-series-of-3d-points-using-cubic-splines#answer_384467)
    %interp1 function and spline method: https://uk.mathworks.com/help/matlab/ref/interp1.html#btwp6lt-1-x
    %change with griddedinterpolation later on

%   locs:   localizatiions to be straightened along the polyline or a curve
%   polyline:   array containing 3D coordinates of the midline curve along
                %   which data will be straightened.

    %   pos:   site centre position
    %   selR:   Radius to preselect poinst within the segment %%%This should be
    %   temporary, I need a better preselection function based on the collect
    %   functions used in the point cloud thinner
    
    
    method='pchip';
    str_mth='linear';
    rho_cut=180;
    rho_cutL=5;
    binwidth=5;

    %from interparc - getting spl coefficients and arclength
    t=1:size(polyline(:,3),1);
    nt=size(polyline(:,3),1);
    ndim=size(polyline,2);
    chordlen = sqrt(sum(diff(polyline,[],1).^2,2));
    chordlen = chordlen/sum(chordlen);
    cumarc = [0;cumsum(chordlen)];
    [arclengthPL,seglenPL]=arclength(polyline(:,1),polyline(:,2),polyline(:,3), method);
    %alpha=18
    if nargin<5
      selR = max(mean(seglenPL)*3,200);
    end
    tst=[0;cumsum(seglenPL)];
    %Getting piecewise polynomials
    spl = cell(1,ndim);
    spld = spl;
    diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
    ddiffarray = [2 0;0 1;0 0];
    dddiffarray = [1;0];
    for i = 1:ndim
      switch method
        case 'pchip'
          spl{i} = pchip(cumarc,polyline(:,i));
        case 'spline'
          spl{i} = spline(cumarc,polyline(:,i));
          nc = numel(spl{i}.coefs);
          if nc < 4
            % just pretend it has cubic segments
            spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
            spl{i}.order = 4;
          end
      end
      xp = spl{i};
      xp.coefs = xp.coefs*diffarray;
      xp.order = 3;
      spld{i} = xp;
      xp.coefs=spld{i}.coefs*ddiffarray;
      xp.order = 2;
      spldd{i} = xp;
      xp.coefs=spldd{i}.coefs*dddiffarray;
      xp.order = 1;
      splddd{i} = xp;
    end
    
    % create storage variables
    out.straight=zeros(size(locs.xnmA,1),1);
    useforanglecal=locs.layer==anglelayer;
    out.p=zeros(size(locs.xnmA,1),1);
    out.side=zeros(size(locs.xnmA,1),1);
    out.DGauss=zeros(nt,5);
    out.SGauss=zeros(nt,3);
    points=[locs.xnmA locs.ynmA locs.znmA];
    newpoints=[locs.xnmA locs.ynmA locs.znmA];
    projections=[locs.xnmA locs.ynmA locs.znmA];
    out.correctedPt=polyline; 
        
    finin=zeros(size(points,1),1);
    colour=zeros(size(points,1),1);
    dist2O=ones(size(points,1),1)*selR;
    AXangle=zeros(nt,1);
    m_locprecnm=zeros(nt,1);
    m_locprecznm=zeros(nt,1);
    [Tan,N,B,kappa,tau] = TNB(polyline(:,1),polyline(:,2),polyline(:,3));
    dudt = zeros(nt,ndim);
    d2ud2t = zeros(nt,ndim);
    d3ud3t = zeros(nt,ndim);
    for L = 1:ndim
        dudt(:,L) = ppval(spld{L},cumarc);
        d2ud2t(:,L) = ppval(spldd{L},cumarc);
        d3ud3t(:,L) = ppval(splddd{L},cumarc);
    end
    curvature_2=sqrt(sum(cross(dudt,d2ud2t,2).^2,2))./(sqrt(sum(dudt.^2,2)).^3);
    torsion_2=dot(dudt,cross(d2ud2t,d3ud3t,2),2)./sum(cross(dudt,d2ud2t,2).^2,2); %%doublecheckformula
    for p=2:1:(size(t,2))
        axis=[1, 0, 0];%axis for now it is written only for X so do not change this
        if p==1
            input=[-Inf cumarc(p) cumarc(p+1)];
        elseif p==size(t,2)
            input=[cumarc(p-1) cumarc(p) Inf];
        else
            input=[cumarc(p-1) cumarc(p) cumarc(p+1)];
        end

        intervalP=[ppval(spl{1}, input).',ppval(spl{2}, input).',ppval(spl{3}, input).'];
        O=intervalP(2,:).';
        start=intervalP(1,:).';
        stop=intervalP(3,:).';
        indexH=insphere(points,O,selR); 
        %%%####################ADD unless
        
        if str_mth=='linear'
            T=stop-O;%[ppval(spld{1}, input(2)).',ppval(spld{2}, input(2)).',ppval(spld{3}, input(2)).']'; %tangent of local functionstop-O;%
            n=d2ud2t(p,1:2)*[p; 1]-O;%N(p,:).';%
            if p==nt
                T=[ppval(spld{1}, input(2)).',ppval(spld{2}, input(2)).',ppval(spld{3}, input(2)).']';
            end
        else
            T=[ppval(spld{1}, input(2)).',ppval(spld{2}, input(2)).',ppval(spld{3}, input(2)).']'; %tangent of local function
            n=N(p,:).';%
        end
        r=T/norm(T); %normalised tangent
        n=n/norm(n);
        %both tangent and normal calculated in this way  are given as if O
        %is [0,0,0] to calculate distance of local points i need to
        %translate points by origin. To save computation preselect poins by
        %incircle function with at least 450 nm. Still very intensive,
        %a better solution is needed. Point X to plane distance is given by
        %dot(X-n,r)

        %calculate distance of start and stop to plane
        overhang=round(mean(seglenPL));
        startD=min(dot(start-O-n,r),dot(stop-O-n,r))-overhang;%will be negative, 30 is margin to avoid loosing points
        stopD=max(dot(start-O-n,r),dot(stop-O-n,r))+overhang;%will be positive
        if p==2
            startD=-150;
        elseif p==size(t,1) || p==size(t,1)-1
            stopD=150;
        end
        index=zeros(size(points,1),1);
        indexA=zeros(size(points,1),1);


        for h=1:size(indexH,1)
            if indexH(h)
                D=dot((points(h,:).')-O-n,r); %%%%%%%%%%
                if D>=startD && D<=stopD

                    if finin(h)==1
                        if dist2O(h)>abs(D)
                            index(h)=1;
                            indexA(h)=1;%%%%%%%%%%%
                        else
                            index(h)=0;
                            indexA(h)=0;%%%%%%%%%% keep already straigh points when calculating angle
                        end
                    else
                        index(h)=1;
                        indexA(h)=1;
                    end
                end
                dist2O(h)=abs(D);

            end

        end

        %rotation to specified axis
        angle= -acosd(dot(r,axis));%- means rotate clockwise!! %atan2d(norm(cross(r,axis)),dot(r,axis)) %LINK missing to matlab and wiki page explaining this
        u=cross(axis,r);
        if norm(u)==0
            u=u;
            rotMat=eye(3);%040423
%         elseif dot(r,axis)==0
%             angle=angle+1;
            if dot(axis,r)<0
                rotMat(1,1)=-1;
            end
        else
            
%             ssc = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];%040423 https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
%             rotMat = eye(3) + ssc + (ssc^2)*(1+dot(axis,r))/(norm(u)^2) / 1;%040423
            u=u/norm(u);
            rotMat=eye(3)*cosd(angle(1))+sind(angle(1))*[0,-u(3), u(2);u(3),0,-u(1);-u(2),u(1) 0]+(1-cosd(angle(1)))*  [u(1)^2, u(1)*u(2), u(1)*u(3); u(1)*u(2), u(2)^2, u(2)*u(3); u(1)*u(3), u(2)*u(3), u(3)^2 ];
        
        end
        
        %Add link to the wiki site for the rotation
        %rotMat=eye(3)*cosd(angle(1))+sind(angle(1))*[0,-u(3), u(2);u(3),0,-u(1);-u(2),u(1) 0]+(1-cosd(angle(1)))*  [u(1)^2, u(1)*u(2), u(1)*u(3); u(1)*u(2), u(2)^2, u(2)*u(3); u(1)*u(3), u(2)*u(3), u(3)^2 ];
            %040423 removed rotmat line
        dist=tst(p);%cumarc(p)*arclengthPL;      
        
        
        if mod(p,2)==0
            newpoints(index==1,:)=(rotMat*((points(index==1,:)-(O.')).')).';
        end
        projections(indexA==1,:)=(rotMat*((points(indexA==1,:)-(O.')).')).';%getAxisRotPreview(rotMat*(r),newpoints(index==1,:),[0 0 0]);
        
        %anglesNX(p,j)=dot(r,stop-O);
        %axisangle(p,j)
        if calangle==1
            %restricty=projections(:,2)<200 & projections(:,2)>-200;%& restricty
            [theta, rho]=cart2pol(projections([indexA==1 & useforanglecal],2),projections([indexA==1 & useforanglecal],3));
            
%             theta=mod(theta,2*pi); %%%%%%%%%%%2*
%             theta(theta<0)=theta(theta<0)+1*pi;
            
%             rho_filter=rho<400;
%             %recenter data for better angle calculation
%             [y,z] = pol2cart(theta(rho_filter),rho(rho_filter));
%             [theta, rho]=cart2pol(y-mean(y),z-mean(z));

            rho_filter=rho<rho_cut & rho>rho_cutL;
            weightsA=1./(locs.locprecnm([indexA==1 & useforanglecal]).^2);
            avAXangle=cyclicaverage(theta(rho_filter),pi,weightsA(rho_filter));
            %%%Add locprecision
%             avAXangle=mod(avAXangle,2*pi);
            AXangle(p)=avAXangle;

            if p==size(t,2)-1
                AXangle(p+1)=avAXangle;
            end
%             avAXangle(avAXangle>pi/2)=avAXangle(avAXangle>pi/2)-pi;
%             avAXangle(avAXangle<0)=abs(avAXangle(avAXangle<0));
            curr_ang=90-rad2deg(avAXangle);%
            
            X=projections(indexA==1 & useforanglecal,[2 3]);
            gm = fitgmdist(X,2,'CovarianceType','diagonal','SharedCovariance',true);
            %P = posterior(gm,X);
            idx = cluster(gm,X);
            GMcC=[[1; 1]  gm.mu(:,1)]\gm.mu(:,2);
            correction=[mean(gm.mu(:,1)); GMcC(2)*mean(gm.mu(:,1))+GMcC(1)];
                       
%             curr_ang=-curr_ang;
            rotMat_ax=[1 0 0; 0 cosd(curr_ang) -sind(curr_ang); 0 sind(curr_ang) cosd(curr_ang)];

            if mod(p,2)==0
                newpoints(index==1,:)=(rotMat_ax*newpoints(index==1,:).').';
            end
            rotationforeval=(rotMat_ax*projections(indexA==1 & useforanglecal,:).').';
            fordistance=rotationforeval(:,3);
            bins=linspace(min(fordistance),max(fordistance),round((max(fordistance)-min(fordistance))/binwidth));
            histD=histogram(fordistance,bins);
            binsC=arrayfun(@(i)(bins(i)+bins(i+1))/2,1:(size(bins,2)-1));
            out.DGauss(p,:)=axisDoubleGaussian_asimetric(binsC,histD.Values,[10,10,20,double(mean(fordistance)),50]).estimate;
            forshift=rotationforeval(:,2);
            bins=linspace(min(forshift),max(forshift),round((max(forshift)-min(forshift))/binwidth));
            histD=histogram(forshift,bins);
            binsC=arrayfun(@(i)(bins(i)+bins(i+1))/2,1:(size(bins,2)-1));
            out.SGauss(p,:)=axisSingleGaussian(binsC,histD.Values,[10,20,50]).estimate;
        end
        if mod(p,2)==0
            newpoints(index==1,1)=newpoints(index==1,1)+dist;
%             if abs(out.DGauss(p,4))>5 || abs(out.SGauss(p,3))>5
%                 newpoints(index==1,3)=newpoints(index==1,3)-out.DGauss(p,4);
%                 newpoints(index==1,2)=newpoints(index==1,2)-out.SGauss(p,3);
%                 out.correctedPt(p,:)=(((rotMat.')*(rotMat_ax.')*([0 out.SGauss(p,3) out.DGauss(p,4)].')).'+(O.')).';
%             end
      
            finin(index==1)=1;
        end
%         oldO=O;
        
        colour(index==1)=p;%j
        out.p(index==1)=p;      
        m_locprecznm(p)=mean(locs.locprecznm(index==1));
        m_locprecnm(p)=mean(locs.locprecnm(index==1));
        
        if p==2%size(t,2)-1
%             m_locprecznm(p+1)=mean(locs.locprecznm(index==1));
%             m_locprecnm(p+1)=mean(locs.locprecnm(index==1));
            m_locprecznm(p-1)=mean(locs.locprecznm(index==1));
            m_locprecnm(p-1)=mean(locs.locprecnm(index==1));
        end
    %     % Additonal filtering
    %     AS=alphaShape(double(locs.xnmA),double(locs.ynmA),double(locs.znmA),alpha,'HoleThreshold',200000);
    %     tri = alphaTriangulation(AS,1);
    %     ASfil=AS;
    %     ASfil.Points=AS.Points(unique(tri),:);
    %     %AT=triangulation(tri,AS.Points);
    %     out.ASfilter=inShape(ASfil,double(locs.xnmA),double(locs.ynmA),double(locs.znmA));
%         out.segments(p)=ns;
        
    end
    locs.xnmS=newpoints(:,1);
    locs.ynmS=newpoints(:,2);
    locs.znmS=newpoints(:,3);

    locs.xnmP=projections(:,1); %points without added distance that aligns them in X (along the length of chromatin) - remove this to save memory!!!, can be recreated from straightene
    locs.ynmP=projections(:,2);
    locs.znmP=projections(:,3);
    
    
    out.straight=logical(finin);
    out.segments=colour;
    out.curvature=kappa;
    out.torsion=tau;
    out.curvature_2=curvature_2;
    out.torsion_2=torsion_2;
    out.arclength=arclengthPL;
    out.seglength=seglenPL;
    out.m_locprecnm=m_locprecnm;
    out.m_locprecznm=m_locprecznm;

    if calangle==1
        out.axesangles=AXangle;
    end
    out.locs=locs;
end

function out=straigthen(locs,polyline,  makesubseg, subseglength, anglelayer, calangle, selR)
% Straighten the SC along the 3D polyline 
    %links to the used functions
    %(https://uk.mathworks.com/matlabcentral/answers/473087-how-do-you-calculate-a-trajectory-through-a-series-of-3d-points-using-cubic-splines#answer_384467)
    %interp1 function and spline method: https://uk.mathworks.com/help/matlab/ref/interp1.html#btwp6lt-1-x
    %change with griddedinterpolation later on

%   locs:   localizatiions to be straightened along the polyline or a curve
%   polyline:   array containing 3D coordinates of the midline curve along
                %   which data will be straightened.

    %   pos:   site centre position
    %   selR:   Radius to preselect poinst within the segment %%%This should be
    %   temporary, I need a better preselection function based on the collect
    %   functions used in the point cloud thinner
    
    %alpha=18
    if nargin<7
      selR = 180;
    end

    t=1:size(polyline(:,3),1);
    x=polyline(:,1);
    y=polyline(:,2);
    z=polyline(:,3);
    
    %Getting piecewise polynomials
    polynomials=struct();
    polynomials.x = interp1(t,x,'spline','pp');
    polynomials.y = interp1(t,y,'spline','pp');
    polynomials.z = interp1(t,z,'spline','pp');
    
    % get new coord
    % Tangent, normal, binormal vectors
    locs.xnmS=locs.xnmA;
    locs.ynmS=locs.ynmA;
    locs.znmS=locs.znmA;
    
    locs.xnmP=locs.xnmA;
    locs.ynmP=locs.ynmA;
    locs.znmP=locs.znmA;
    out.straight=zeros(size(locs.xnmA,1),1);

    useforanglecal=locs.layer==anglelayer;
    
    out.p=zeros(size(locs.xnmA,1),1);
    out.j=zeros(size(locs.xnmA,1),1);
    
    points=[locs.xnmA locs.ynmA locs.znmA];
    newpoints=[locs.xnmA locs.ynmA locs.znmA];
    projections=[locs.xnmA locs.ynmA locs.znmA];
    
    
    finin=zeros(size(points,1),1);
    line=[x y z];
    colour=zeros(size(points,1),1);
    dist2O=ones(size(points,1),1)*selR;
    AXangle=[];
    
    for p=1:(size(t,2)-1)
        axis=[1, 0, 0];%axis for now it is written only for X so do not change this
        if makesubseg==1
            ns=round(abs(norm(line(p,:)-line(p+1,:)))/subseglength); %%%%%%%%%%%%%%%%%Need to remove subsegment since I will be using locmofit's cspline
            if ns==0
                ns=1;
            end
        else
            ns=1;
        end
        
        subseg=linspace(t(p),t(p+1),ns+1);%
        
        for j=1:(size(subseg,2)-1)
    
            center=(subseg(j)+subseg(j+1))/2;
            interval=linspace(subseg(j),subseg(j+1),20);
            pX = polyfit(interval,ppval(polynomials.x, interval),3);
            pY = polyfit(interval,ppval(polynomials.y, interval),3);
            pZ = polyfit(interval,ppval(polynomials.z, interval),3);
            d1=[[3 2 1 0].*pX; [3 2 1 0].*pY; [3 2 1 0].*pZ];
            d2=[[6 2 0 0].*pX; [6 2 0 0].*pY; [6 2 0 0].*pZ];
            O=[polyval(pX, center);polyval(pY, center);polyval(pZ, center)];
            start=[polyval(pX, subseg(j));polyval(pY, subseg(j));polyval(pZ, subseg(j))];
            stop=[polyval(pX, subseg(j+1));polyval(pY, subseg(j+1));polyval(pZ, subseg(j+1))];
            indexH=insphere(points,O,selR); %% FIX SELECTION FUNCTION, this is just for testing it is not good!!!!!
    
    
            %%%####################ADD unless
    
            %point selection
            n=d2(:,1:2)*[center; 1]-O; %normal of the local funtion, also a point on the plane perpendicular to the tangent
            n=n/norm(n);
            T=stop-O;%T=d1(:,1:3)*[center*center; center ;1]-O; %tangent of local function
            r=T/norm(T); %normalised tangent
            %both tangent and normal calculated in this way  are given as if O
            %is [0,0,0] to calculate distance of local points i need to
            %translate points by origin. To save computation preselect poins by
            %incircle function with at least 900 nm. Still very intensive,
            %a better solution is needed. Point X to plane distance is given by
            %dot(X-n,r)
    
            %calculate distance of start and stop to plane
            startD=min(dot(start-O-n,r),dot(stop-O-n,r));%will be negative
            stopD=max(dot(start-O-n,r),dot(stop-O-n,r));%will be positive
            index=zeros(size(points,1),1);
    
    
            for h=1:size(indexH,1)
                if indexH(h)
                    D=dot((points(h,:).')-O-n,r);
                    if D>=startD && D<=stopD
    
                        if finin(h)==1
                            if dist2O(h)>abs(D)
                                index(h)=1;
                            else
                                index(h)=0;
                            end
                        else
                            index(h)=1;
                        end
                    end
                    dist2O(h)=abs(D);
    
                end
    
            end
    
            %rotation to specified axis
            angle= -acosd(dot(r,axis));%- means rotate clockwise!! %atan2d(norm(cross(r,axis)),dot(r,axis)) %LINK missing to matlab and wiki page explaining this
            u=cross(axis,r);
            u=u/norm(u);
            
            %Add link to the wiki site for the rotation
            rotMat=eye(3)*cosd(angle(1))+sind(angle(1))*[0,-u(3), u(2);u(3),0,-u(1);-u(2),u(1) 0]+(1-cosd(angle(1)))*  [u(1)^2, u(1)*u(2), u(1)*u(3); u(1)*u(2), u(2)^2, u(2)*u(3); u(1)*u(3), u(2)*u(3), u(3)^2 ];
    
            % I am missing the rotation of the polyline
            if (j==1 && p==1)
                dist=0;
%                 newpoints(index==1,:)=(rotMat*((points(index==1,:)-(O.')).')).';
%                 projections(index==1,:)=newpoints(index==1,:);%getAxisRotPreview(rotMat*(r),newpoints(index==1,:),[0 0 0]);
%                 newpoints(index==1,1)=newpoints(index==1,1)+dist;
    
            else
                dist=dist+abs(norm(O-oldO));
                step=abs(norm(O-oldO));
            end
    
            
            rho_cut=150;
            newpoints(index==1,:)=(rotMat*((points(index==1,:)-(O.')).')).';
            projections(index==1,:)=newpoints(index==1,:);%getAxisRotPreview(rotMat*(r),newpoints(index==1,:),[0 0 0]);
            
            %anglesNX(p,j)=dot(r,stop-O);
            %axisangle(p,j)
            if calangle==1
                %restricty=projections(:,2)<200 & projections(:,2)>-200;%& restricty
                [theta rho]=cart2pol(projections([index==1 & useforanglecal],2),projections([index==1 & useforanglecal],3));
                
                
                %theta(theta<0)=theta(theta<0)+1*pi; %%%%%%%%%%%2
                theta=mod(theta,2*pi);
                theta(theta>pi/2 & theta<3*pi/2)=theta(theta>pi/2 & theta<3*pi/2)+pi;
                theta(theta>=3*pi/2 & theta<2*pi)=theta(theta>=3*pi/2 & theta<=2*pi)-pi;
                theta(theta>pi/2)=theta(theta>pi/2)-pi;
                %theta(theta>=3*pi/2 & theta<2*pi)=2*pi-theta(theta>=3*pi/2 & theta<2*pi);
                %theta(theta>2*pi)=theta(theta>2*pi)-2*pi;
%             theta=mod(theta,2*pi);
%             theta(theta>pi/2 & theta<3*pi/2)=theta(theta>pi/2 & theta<3*pi/2)+pi;
%             theta(theta>2*pi)=theta(theta>2*pi)-2*pi;
%             theta(theta>=3*pi/2 & theta<2*pi)=theta(theta>=3*pi/2 & theta<=2*pi)-pi;
%             theta(theta>pi/2)=theta(theta>pi/2)-pi;
                %theta(theta>=3*pi/2 & theta<2*pi)=2*pi-theta(theta>=3*pi/2 & theta<2*pi);
%             theta(theta>pi/2)=theta(theta>pi/2)-pi;
            %theta(theta>=3*pi/2 & theta<2*pi)=2*pi-theta(theta>=3*pi/2 & theta<2*pi);
            %theta(theta>2*pi)=theta(theta>2*pi)-2*pi;
                weightsA=1./(locs.locprecnm([index==1 & useforanglecal]).^2);
                avAXangle=cyclicaverage(abs(theta(rho<rho_cut)),pi,weightsA(rho<rho_cut));
                AXangle(end+1)=avAXangle;%%%Add locprecision
                curr_ang=90-rad2deg(avAXangle);
                rotMat_ax=[1 0 0; 0 cosd(curr_ang) -sind(curr_ang); 0 sind(curr_ang) cosd(curr_ang)];
                newpoints(index==1,:)=(rotMat_ax*newpoints(index==1,:).').';
            end
            newpoints(index==1,1)=newpoints(index==1,1)+dist;
            oldO=O;
            finin(index==1)=1;
            colour(index==1)=p;%j
            out.p(index==1)=p;
            out.j(index==1)=j;
    
    
        end
%         locs.xnmS=newpoints(:,1);
%         locs.ynmS=newpoints(:,2);
%         locs.znmS=newpoints(:,3);
%     
%         locs.xnmP=projections(:,1); %points without added distance that aligns them in X (along the length of chromatin) - remove this to save memory!!!, can be recreated from straightene
%         locs.ynmP=projections(:,2);
%         locs.znmP=projections(:,3);
%         
%         
%         out.straight=logical(finin);
%         out.segments=colour;
        if calangle==1
            out.axesangles=AXangle;
        end

    %     % Additonal filtering
    %     AS=alphaShape(double(locs.xnmA),double(locs.ynmA),double(locs.znmA),alpha,'HoleThreshold',200000);
    %     tri = alphaTriangulation(AS,1);
    %     ASfil=AS;
    %     ASfil.Points=AS.Points(unique(tri),:);
    %     %AT=triangulation(tri,AS.Points);
    %     out.ASfilter=inShape(ASfil,double(locs.xnmA),double(locs.ynmA),double(locs.znmA));
        out.segments(p)=ns;
        
    end
    locs.xnmS=newpoints(:,1);
    locs.ynmS=newpoints(:,2);
    locs.znmS=newpoints(:,3);

    locs.xnmP=projections(:,1); %points without added distance that aligns them in X (along the length of chromatin) - remove this to save memory!!!, can be recreated from straightene
    locs.ynmP=projections(:,2);
    locs.znmP=projections(:,3);
    
    
    out.straight=logical(finin);
    out.segments=colour;
    out.locs=locs;
end



function index = insphere(data,position, dist)
    index=((data(:,1)-position(1)).^2+(data(:,2)-position(2)).^2+(data(:,3)-position(3)).^2<=dist^2);
end

function out=getAngle(data,method)
%
end

function out=getAxisRotPreview(normal,points,origin)
% This finction is used to project points to a plane that is given by it's
% normal. Normmal corresponds to the tangent of the SC midline, and
% function returns points projected to the plane that is perpendicular to
% the straightened axes
%points is a nx3 matrix
%normal is 3x1 column vector
%origin is 1x3 row vector specifying point on the plane and origin of the
%normal
%out is nX3 matrix

%normalize normal
n=normal/norm(normal);
proj=points - ((points-origin)*n)*n.';
out.xnm=proj(:,1);
out.ynm=proj(:,2);
out.znm=proj(:,3);
end

function out=RotMatrix3DvectortoX(v)
axis=[1,0,0];
r=v/norm(v);
angle= -acosd(dot(r,axis));%- means rotate clockwise!! %atan2d(norm(cross(r,axis)),dot(r,axis)) %LINK missing to matlab and wiki page explaining this
u=cross(axis,r);
u=u/norm(u);
out=eye(3)*cosd(angle(1))+sind(angle(1))*[0,-u(3), u(2);u(3),0,-u(1);-u(2),u(1) 0]+(1-cosd(angle(1)))*  [u(1)^2, u(1)*u(2), u(1)*u(3); u(1)*u(2), u(2)^2, u(2)*u(3); u(1)*u(3), u(2)*u(3), u(3)^2 ];
end

function [ctrlX,ctrlY,ctrlZ] = getCtrlPointsXYZloc(allpars)
% control points:
lx = startsWith(allpars.name,'cx');
ctrlSetID = str2double(regexprep(allpars.name(lx),'\D',''));
lastCtrlPoint_ori = max(ctrlSetID);
par=[];
for f=1:size(allpars.name,1)
    par.(allpars.name{f})=allpars.value(f);
end
      
    for k = 1:lastCtrlPoint_ori
        ctrlX(k,:) = par.(['cx' num2str(k)]); % ->
        ctrlY(k,:) = par.(['cy' num2str(k)]); % ->
        ctrlZ(k,:) = par.(['cz' num2str(k)]); % ->
    end
    
end