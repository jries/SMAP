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

            createsubseg=obj.getSingleGuiParameter('createsubsegT');
            subseglength=obj.getSingleGuiParameter('createsubseg');

            if ~isfield(obj.site.evaluation,polylinesource)
                error('Please choose a polyline source that is available. %s has not been evaulated',polylinesource)
           
            end
            layerson = find(obj.locData.getPar('sr_layerson'));             % check the layers used
            visitFlag = false;                                      % a flag for the conditional loop for the first round
            fieldQ = {'locprecnm','locprecznm','PSFxnm','xnm','znm','ynm','frame','xnmrot','ynmrot','phot1','phot2','filenumber','channel'};    % fields will be used %%%%%%%%%%%%%%%%%%% need to add more for displayer
           
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
            descriptors=[];
            [ctrlX,ctrlY,ctrlZ]=getCtrlPointsXYZloc(obj.site.evaluation.(polylinesource).allParsArg);
            [pt,dudt]= interparc(100, ctrlX, ctrlY, ctrlZ);
            descriptors.pt=pt;
            [L,R,K] = curvature(pt);
            descriptors.curvature = 1./R;
            descriptors.ctrlpoints.ctrlX = ctrlX;
            descriptors.ctrlpoints.ctrlY = ctrlY;
            descriptors.ctrlpoints.ctrlZ = ctrlZ;
            
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
            straightner=straigthen(locs,[descriptors.ctrlpoints.ctrlX descriptors.ctrlpoints.ctrlY descriptors.ctrlpoints.ctrlZ],createsubseg,subseglength, axeslayer, calangle);
            warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale'); 

            out.straightner=straightner;
            out.cspline=descriptors.pt;
            out.curvature=descriptors.curvature;
            out.ctrlpoints=descriptors.ctrlpoints;
%            out.csplineDudt=descriptors.dudt;

            if display
                figure('Name','3dpreview - localizations and polyline');
                hold on
                scatter3(descriptors.ctrlpoints.ctrlX,descriptors.ctrlpoints.ctrlY,descriptors.ctrlpoints.ctrlZ,'green','filled')
                scatter3(locs.xnmA, locs.ynmA, locs.znmA,[2],locs.layer,'filled'); % 
                colormap("flag")
                daspect([1 1 1])
                hold off

                figure('Name','3dpreview - straightened localizations');
                hold on
                scatter3(straightner.locs.xnmS(straightner.straight), straightner.locs.ynmS(straightner.straight),straightner.locs.znmS(straightner.straight),[2],straightner.locs.layer(straightner.straight),'filled'); % 
                colormap("flag")
                daspect([1 1 1])
                hold off
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

               CSV=[newpath '\\' path{end}(1:end-4) obj.site.name(3:end-1) '_' dateC{1} '_straightenedSC.csv'];
               final=struct2table(locsS);
               writetable(final,CSV,'WriteRowNames',true);
           end

           if angleplot==1 && calangle==1
               figure('Name','Angle between axes');
               hold on
               plot(straightner.axesangles); % 
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
      selR = 250;
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
                newpoints(index==1,:)=(rotMat*((points(index==1,:)-(O.')).')).';
                projections(index==1,:)=newpoints(index==1,:);%getAxisRotPreview(rotMat*(r),newpoints(index==1,:),[0 0 0]);
                newpoints(index==1,1)=newpoints(index==1,1)+dist;
    
            else
                dist=dist+abs(norm(O-oldO));
                step=abs(norm(O-oldO));
                newpoints(index==1,:)=(rotMat*((points(index==1,:)-(O.')).')).';
                projections(index==1,:)=newpoints(index==1,:);%getAxisRotPreview(rotMat*(r),newpoints(index==1,:),[0 0 0]);
                newpoints(index==1,1)=newpoints(index==1,1)+dist;
            end
    
    
            oldO=O;
            finin(index==1)=1;
            colour(index==1)=p;%j
            %anglesNX(p,j)=dot(r,stop-O);
            %axisangle(p,j)
            if calangle==1
                restricty=projections(:,2)<200 & projections(:,2)>-200;
                [theta rho]=cart2pol(projections([index==1 & useforanglecal & restricty],2),projections([index==1 & useforanglecal & restricty],3));
                
                theta(theta<0)=theta(theta<0)+1*pi; %%%%%%%%%%%2
                AXangle(end+1)=rad2deg(cyclicaverage(theta(rho<200),pi));%%%Add locprecision
            end

            out.p(index==1)=p;
            out.j(index==1)=j;
    
    
        end
        locs.xnmS=newpoints(:,1);
        locs.ynmS=newpoints(:,2);
        locs.znmS=newpoints(:,3);
    
        locs.xnmP=projections(:,1); %points without added distance that aligns them in X (along the length of chromatin) - remove this to save memory!!!, can be recreated from straightene
        locs.ynmP=projections(:,2);
        locs.znmP=projections(:,3);
        
        
        out.straight=logical(finin);
        out.segments=colour;
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
        out.locs=locs;
    end
end

function index = insphere(data,position, dist)
    index=((data(:,1)-position(1)).^2+(data(:,2)-position(2)).^2+(data(:,3)-position(3)).^2<=dist^2);
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