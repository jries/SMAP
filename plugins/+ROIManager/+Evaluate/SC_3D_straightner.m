classdef SC_3D_straightner<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
        roicoordinates
        roihandles={[]};
        roisize
        roisizey
%         currentroi=0;
        axis
        images
%         roihandle
    end
    methods
        function obj=SC_3D_straightner(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            out=[];
            polylinesource=obj.getSingleGuiParameter('polylinesource');
            expandPL=obj.getSingleGuiParameter('polylineexpandT');
            expandby=obj.getSingleGuiParameter('polylineexpand');

            if ~isfield(obj.site.evaluation,polylinesource)
                error('Please choose a polyline source that is available. %s has not been evaulated',polylinesource)
           
            end
            layerson = find(obj.locData.getPar('sr_layerson'));             % check the layers used
            visitFlag = false;                                      % a flag for the conditional loop for the first round
            fieldQ = {'locprecnm','locprecznm','PSFxnm','xnm','znm','ynm','frame','xnmrot','ynmrot'};    % fields will be used %%%%%%%%%%%%%%%%%%% need to add more for displayer
           
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
            %%polyline

            derivedparameters=obj.site.evaluation.(polylinesource).GuiParameters.fitter.getDerivedPars;
            descriptors=derivedparameters{1,1};%derivePars(obj,obj.site.evaluation.(polylinesource).allParsArg,polylinesource);%
            
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

            figure('Name','3dpreview');
            hold on
            scatter3(descriptors.pt(:,1),descriptors.pt(:,2),descriptors.pt(:,3),[])
            scatter3(locs.xnmA, locs.ynmA, locs.znmA,[2],locs.znmA,'filled'); % 
            colorbar
            daspect([1 1 1])
            hold off
               %out.
            straightner=straigthen(locs,descriptors.pt);


            figure('Name','3dpreview');
            hold on
            scatter3(straightner.locs.xnmS(straightner.straight), straightner.locs.ynmS(straightner.straight),straightner.locs.znmS(straightner.straight),[2],straightner.locs.layer(straightner.straight),'filled'); % 
            colorbar
            colormap("flag")
            daspect([1 1 1])
            hold off


           
            locsOri = locs;
           

            make3Daxis(obj) %make axis,
            plot3Dviews(obj,p);%plot image,
            if ~isfield(obj.site.annotation,'polarangle')
                obj.site.annotation.polarangle=0;
            end
%             obj.currentroi=0;
            if isfield(obj.site.evaluation,obj.name) && isfield(obj.site.evaluation.(obj.name),'roicoordinates')
                obj.roicoordinates=obj.site.evaluation.(obj.name).roicoordinates;
                obj.plotrois;
                out.GuiParameters=obj.site.evaluation.(obj.name).GuiParameters; %ICfix_210515 --- CHECK 
                out.roicoordinates=obj.site.evaluation.(obj.name).roicoordinates; %ICfix_210515 --- CHECK 
                %out=obj.site.evaluation.(obj.name);  - use this above
                %instead of two lines
            else
                obj.roicoordinates=[];
               
            end
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function addroi_callback(obj,a,b)
            roistyle=obj.getSingleGuiParameter('roiform').selection;
            [roifun,roifunmake,vertexbased,par]=roistyle2fun(roistyle);
            
            newroi1=roifun(obj.axis,'DrawingArea',[-obj.roisize -obj.roisizey obj.roisize*2 obj.roisizey*2],'Color','w',par{:});
            [roipospix,roipos]=getpixcoord(obj,newroi1.Position,zeros(size(newroi1.Position,1),3));
          
            roiposnm=coord3Dpix2nm(roipospix,obj.site.pos,obj.site.annotation.rotationpos.angle,obj.site.annotation.polarangle,obj.roisize);
            
%             roiposnm=roipospix; %later: convert with rotation and position
%             obj.currentroi=obj.currentroi+1;  
            obj.roicoordinates(end+1).roistyle=roistyle;
            obj.roicoordinates(end).Position=roiposnm;
            delete(newroi1);
            obj.plotrois;         
        end
        function deleteroi_callback(obj,a,b)
            obj.roicoordinates=[];
%             obj.plotrois; 
            plot3Dviews(obj);
        end
        
        function plotrois(obj)
            if isempty(obj.roicoordinates)
                return
            end
             %loop over roicoords
%              k=obj.currentroi;
            %delete old rois
            for k=1:length(obj.roihandles{1})
                delete(obj.roihandles{1}(k));
                delete(obj.roihandles{2}(k));
            end
            
            colors=lines(255);
            
            for k=1:length(obj.roicoordinates)
                 [roifun,roifunmake,vertexbased,par]=roistyle2fun(obj.roicoordinates(k).roistyle);
                 roiposnm=obj.roicoordinates(k).Position;
                 roipospix=coord3Dnm2pix(roiposnm,obj.site.pos,obj.site.annotation.rotationpos.angle,obj.site.annotation.polarangle,obj.roisize);
                 roipos1=[];roipos2=[];
                 roipos1(:,1:2)=roipospix(:,1:2);
                 roipos2(:,1)=roipospix(:,1)+obj.images.rangex(2)*1000;
                 roipos2(:,2)=roipospix(:,3);

                newroi1(k)=roifunmake(obj.axis,'Position',roipos1,'DrawingArea',[-obj.roisize -obj.roisizey obj.roisize obj.roisizey*2],'Color',colors(k,:),par{:});
                addlistener(newroi1(k),'ROIMoved',@(src,evt) obj.updateposition(src,evt,k,1));
                addlistener(newroi1(k),'DeletingROI',@(src,evt) obj.deleteroi(src,evt,k,1));
                if vertexbased
                    addlistener(newroi1(k),'VertexAdded',@(src,evt) obj.updateposition(src,evt,k,1));
                    addlistener(newroi1(k),'VertexDeleted',@(src,evt) obj.updateposition(src,evt,k,1));
                end

                newroi2(k)=roifunmake(obj.axis,'Position',roipos2,'DrawingArea',[0 -obj.roisizey obj.roisize obj.roisizey*2],'Color',colors(k,:),par{:});
                addlistener(newroi2(k),'ROIMoved',@(src,evt) obj.updateposition(src,evt,k,2));
                addlistener(newroi2(k),'DeletingROI',@(src,evt) obj.deleteroi(src,evt,k,2));
                if vertexbased
                    addlistener(newroi2(k),'VertexAdded',@(src,evt) obj.updateposition(src,evt,k,2));
                    addlistener(newroi2(k),'VertexDeleted',@(src,evt) obj.updateposition(src,evt,k,2));
                end           
            end
            obj.roihandles={newroi1, newroi2};
            obj.site.evaluation.(obj.name).roicoordinates=obj.roicoordinates;
        end
        
        function updateposition(obj,src,evt,currentroi,roiside)
            
            
            oldposnm=obj.roicoordinates(currentroi).Position; %convert
            oldpospix=coord3Dnm2pix(oldposnm,obj.site.pos,obj.site.annotation.rotationpos.angle,obj.site.annotation.polarangle,obj.roisize);
            
            newpospix=src.Position;
            %vertex added or deleted: oldpos has different number
            if size(oldpospix,1) > size(newpospix,1) %vertex deleted
                ind=findextra(oldpospix,newpospix);
                oldpospix(ind,:)=[];
            elseif size(oldpospix,1) < size(newpospix,1) %vertex added
                ind=findextra(newpospix,oldpospix);
                oldpospix=[oldpospix(1:ind-1,:); zeros(1,3); oldpospix(ind:end,:)];
            end
            
            [roipospix,roipos]=getpixcoord(obj,newpospix,oldpospix);
            
            newroiposnm=coord3Dpix2nm(roipospix,obj.site.pos,obj.site.annotation.rotationpos.angle,obj.site.annotation.polarangle,obj.roisize);
            obj.roicoordinates(currentroi).Position=newroiposnm;
            plotrois(obj)

        end
        
        function deleteroi(obj,src,evt,k,roiside)
            obj.roicoordinates(k)=[];
            obj.site.evaluation.(obj.name).roicoordinates(k)=[];
            plotrois(obj)
        end

    end

end

function ind=findextra(A,B)
xA=A(:,1);
xB=B(:,1);
mind=min((xA-xB').^2,[],2);
[~,ind]=max(mind);
end

function [roipospix,roipos]=getpixcoord(obj,roiposh,roiposold)
roipospix=roiposold;
 if mean(roiposh(:,1))>0
    roipos=2; %right
    roiposh(:,1)=max(roiposh(:,1),0);
    roipospix(:,1)=roiposh(:,1)-obj.roisize;
    roipospix(:,3)=roiposh(:,2); %move to left
%     roipospix(:,2)=0;
else
    roipos=1;
    roiposh(:,1)=min(roiposh(:,1),0);
    roipospix(:,1:2)=roiposh(:,1:2);
%     roipospix(:,3)=0;
end
end

function roiposnm=coord3Dpix2nm(roipos,pos,rotationanglez,polarangle,winsize)
    xp=roipos(:,1)+winsize/2;
    yp=roipos(:,2);zp=roipos(:,3);
    [y2,zi]=rotcoorddeg(yp,zp,polarangle);
    [xi,yi]=rotcoorddeg(xp,y2,rotationanglez);
    roiposnm(:,1)=xi+pos(1); 
    roiposnm(:,2)=-yi+pos(2); %%%%220329 -yi
    roiposnm(:,3)=zi+pos(3);
end
function roipos=coord3Dnm2pix(roiposnm,pos,rotationanglez,polarangle,winsize)
    xi=roiposnm(:,1)-pos(1); 
    yi=-(roiposnm(:,2)-pos(2)); %%%%220329 -(...)
    zi=roiposnm(:,3)-pos(3);
    [xp,y2]=rotcoorddeg(xi,yi,-rotationanglez); %-
    [yp,zp]=rotcoorddeg(y2,zi,-polarangle);
    roipos(:,1)=xp-winsize/2;
    roipos(:,2)=yp;
    roipos(:,3)=zp;
end

function pard=guidef(obj)
pard.addroi.object=struct('Style','pushbutton','String','add','Callback',@obj.addroi_callback);
pard.addroi.position=[1,1];
pard.addroi.Width=1;

pard.roiform.object=struct('Style','popupmenu','String',{{'polyline','line','free','polygon'}});
pard.roiform.position=[1,2];
pard.roiform.Width=2;


pard.deleterois.object=struct('Style','pushbutton','String','delete all rois','Callback',@obj.deleteroi_callback);
pard.deleterois.position=[2,1];
pard.deleterois.Width=2;

pard.polylinesourceT.object=struct('String','3D polyline source:','Style','text');
pard.polylinesourceT.position=[3,1];
pard.polylinesourceT.Width=2;

pard.polylinesource.object=struct('String','LocMoFitGUI_2','Style','edit');
pard.polylinesource.position=[3,2.5];
pard.polylinesource.Width=2;
pard.polylinesource.TooltipString='Specified previously run evaluation that will serve as a  source for the polyline. Avaliable LocMoFit and plot3Dannotation(still in progress).';

pard.polylineexpandT.object=struct('String','Expand polyline to mask by (nm):','Style','checkbox','Value',0);
pard.polylineexpandT.position=[4,1];
pard.polylineexpandT.Width=3;
pard.polylineexpand.TooltipString='In case there is no user specified mask within annotation tab this will be run irrespective of the status here.';

pard.polylineexpand.object=struct('String','350','Style','edit');
pard.polylineexpand.position=[4,4.5];
pard.polylineexpand.Width=0.5;
pard.polylineexpand.TooltipString='Expand around polyline by nm. Useful when ROI contains many localizations that do not belong to the SC stretch that will be straigthened';


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

function plot3Dstraigtened(obj,p,in)
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



function plot3Dviews(obj,p)
if nargin<2
    p=[];
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


function make3Daxis(obj)
ax=obj.setoutput('orthogonal_views');
ax.NextPlot='replacechildren';
ax.ButtonDownFcn={@sideview_click,obj};
f=ax.Parent;
f.Units='normalized';
uicontrol(f,'Style','pushbutton','String','reset','Units','normalized','Position',[0.9,0.05,0.08,0.05],'Callback',{@resetview,obj,obj.site})
uicontrol(f,'Style','pushbutton','String','Info','Units','normalized','Position',[0.05,0.05,0.08,0.05],'Callback',@info)
obj.axis=ax;
end
        

function sideview_click(hax,dat,obj)
site=obj.site;
dx=(obj.images.rangex(2)-obj.images.rangex(1))/2*1000;
dy=(obj.images.rangey(2)-obj.images.rangey(1))*1000;
pos=dat.IntersectionPoint;
if pos(1)>0 && pos(1)<dx*0.75 %only right side: side view%
    site.pos(3)=site.pos(3)+pos(2);
elseif pos(1)>dx*.75  %outside: rotate polarangle
    if ~isfield(site.annotation,'polarangle')
        site.annotation.polarangle=0;
    end
    anglenew=site.annotation.polarangle+pos(2)/dy/pi/2*180;
    anglered=mod(anglenew+180,360)-180;
    site.annotation.polarangle=anglered;
    obj.setPar('se_currentPolarAngle',anglered);
elseif pos(1)<0 %rotate rotation angle
    angle=atan2d(pos(2),(dx/2+pos(1)))+180;
    anglenew=site.annotation.rotationpos.angle-angle;
    anglered=mod(anglenew+180,360)-180;
    site.setlineangle(0,anglered);
end
    
plot3Dviews(obj)
obj.plotrois
end

function displayimage(img,hax)
if isempty(img)
    return
end
 imagesc(img.rangex*1000,img.rangey*1000,img.image,'Parent',hax,'Pickable','none','HitTest','off')

set(hax,'Xlim',double(img.rangex*1000))
set(hax,'Ylim',double(img.rangey*1000))
set(hax,'YDir','reverse')
hax.HitTest='on';
hax.PickableParts='all';

axis(hax,'equal')
set(hax,'YDir','normal')
axis(hax,'tight')
hax.XTick=0;
hax.YTick=0;
hax.XTickLabel={};
hax.YTickLabel={};
hax.TickDir='out';
hax.YAxisLocation='right';
hax.Box='on';

end


function [roifun,roifunmake,vertexbased,par]=roistyle2fun(roistyle)
vertexbased=false;
par={};
switch roistyle
    case 'rectangle'
        roifun=@drawrectangle;
        roifunmake=@images.roi.Rectangle;
       
    case 'ellipse'
        roifun=@drawellipse;
        roifunmake=@images.roi.Ellipse;
    case 'polygon'
        roifun=@drawpolygon;
        roifunmake=@images.roi.Polygon;
        vertexbased=true;
    case 'free'
        roifun=@drawfreehand;
        roifunmake=@images.roi.Freehand;
         par={'Closed',false};
    case 'polyline'
        roifun=@drawpolyline;  
        roifunmake=@images.roi.Polyline;
                vertexbased=true;
    case 'line'
        roifun=@drawline; 
        roifunmake=@images.roi.Line;

end
end
function resetview(a,b,obj,site)
site.setlineangle(0,0);
site.annotation.polarangle=0;
site.pos(3)=0;
obj.setPar('se_currentPolarAngle',0)
plot3Dviews(obj)
obj.plotrois
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

function out=straigthen(locs,polyline, selR)
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
    if nargin<4
      selR = 500;
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
    
    out.p=zeros(size(locs.xnmA,1),1);
    out.j=zeros(size(locs.xnmA,1),1);
    
    points=[locs.xnmA locs.ynmA locs.znmA];
    newpoints=[locs.xnmA locs.ynmA locs.znmA];
    projections=[locs.xnmA locs.ynmA locs.znmA];
    
    
    finin=zeros(size(points,1),1);
    line=[x y z];
    colour=zeros(size(points,1),1);
    dist2O=ones(size(points,1),1)*selR;
    
    
    for p=1:(size(t,2)-1)
        axis=[1, 0, 0];%axis for now it is written only for X so do not change this
        ns=round(abs(norm(line(p,:)-line(p+1,:)))/5); %%%%%%%%%%%%%%%%%Need to remove subsegment since I will be using locmofit's cspline
        if ns==0
            ns=1;
        end
        subseg=linspace(t(p),t(p+1),ns);
        
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
            anglesNX(p,j)=dot(r,stop-O);
            out.p(index==1)=p;
            out.j(index==1)=j;
    
    
        end
        locs.xnmS=newpoints(:,1);
        locs.ynmS=newpoints(:,2);
        locs.znmS=newpoints(:,3);
    
        locs.xnmP=projections(:,1); %ponits without added distance that aligns them in X (along the length of chromatin)
        locs.ynmP=projections(:,2);
        locs.znmP=projections(:,3);
        
        
        out.straight=logical(finin);
        out.segments=colour;
        
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

% function derivedPars = derivePars(obj, allpars,locmofit) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Only a temporary fix so that I do not have to run locmofit again every time
%          
%         % control points:
%         par=[];
%         for f=1:size(allpars.name,1)
%             par.(allpars.name{f})=allpars.value(f);
%         end
%         [ctrlX,ctrlY,ctrlZ] = getCtrlPointsXYZ_loc(par, obj.site.evaluation.(locmofit).fitInfo.modelPar_internal{1, 1}.numOfCtrlPointSet);
%         
%         if isempty(obj.site.evaluation.(locmofit).GuiParameters.fitter.model{1, 1}.locsPrecFactor)
%             locsPrecFactor = 1;
%         else
%             locsPrecFactor = obj.site.evaluation.(locmofit).GuiParameters.fitter.model{1, 1}.locsPrecFactor;
%         end
%         
%         minD = locsPrecFactor*0.75;
%         
%         arcLen = arclength(ctrlX, ctrlY, ctrlZ,'linear');
%         samplingFactor = round(arcLen/minD);
%         [pt,dudt] = interparc(round(samplingFactor/2), ctrlX, ctrlY, ctrlZ);
%         
%         [L,R,K] = curvature(pt);
%         derivedPars.pt=pt;
%         derivedPars.dudt=dudt;
%         derivedPars.curvature = 1./R;
%         
%         derivedPars.curvatureVector = K;
%         derivedPars.curvatureRadius = R;
%         derivedPars.arclength = L;
%         
%         derivedPars.avgCurvature = mean(derivedPars.curvature,'omitnan');
%         derivedPars.p75Curvature = prctile(derivedPars.curvature, 75);
% %             figure; plot3(pt(:,1),pt(:,2),pt(:,3))
% %             % plot3([pt(:,1) pt(:,1)+50000*K(:,1)]',[pt(:,2) pt(:,2)+50000*K(:,2)]',[pt(:,3) pt(:,3)+50000*K(:,3)]')
% %             axis equal
% %             hold on; arrow3(pt,pt+K.*50000,'-k',0.2,0.4)
% 
% 
% end
% function [ctrlX,ctrlY,ctrlZ] = getCtrlPointsXYZ_loc(par, numOfCtrlPointSet)
% % control points:
%       
%     for k = 1:numOfCtrlPointSet
%         ctrlX(k,:) = par.(['cx' num2str(k)]); % ->
%         ctrlY(k,:) = par.(['cy' num2str(k)]); % ->
%         ctrlZ(k,:) = par.(['cz' num2str(k)]); % ->
%     end
%     
% end