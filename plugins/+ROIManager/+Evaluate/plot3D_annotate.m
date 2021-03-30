classdef plot3D_annotate<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
        roicoordinates
        roihandles
%         currentroi=0;
        axis
        images
%         roihandle
    end
    methods
        function obj=plot3D_annotate(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            out=[];
            make3Daxis(obj) %make axis,
            plot3Dviews(obj,p);%plot image,
            
%             obj.currentroi=0;
            if isfield(obj.site.evaluation,obj.name) && isfield(obj.site.evaluation.(obj.name),'roicoordinates')
                obj.roicoordinates=obj.site.evaluation.(obj.name).roicoordinates;
                obj.plotrois;
            else
                obj.roicoordinates=[];
            end
        end
     
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function addroi_callback(obj,a,b)
            roistyle=obj.getSingleGuiParameter('roiform').selection;
            [roifun,roifunmake,vertexbased]=roistyle2fun(roistyle);
            
            newroi1=roifun(obj.axis);
            [roipospix,roipos]=getpixcoord(obj,newroi1.Position,zeros(size(newroi1.Position,1),3));
            roiposnm=roipospix; %later: convert with rotation and position
%             obj.currentroi=obj.currentroi+1;  
            obj.roicoordinates(end+1).roistyle=roistyle;
            obj.roicoordinates(end).Position=roiposnm;
            delete(newroi1);
            obj.plotrois;         
        end
        
        function plotrois(obj)
             %loop over roicoords
%              k=obj.currentroi;
            %delete old rois
            for k=1:length(obj.roihandles{1})
                delete(obj.roihandles{1}(k));
                delete(obj.roihandles{2}(k));
            end
            
            for k=1:length(obj.roicoordinates)
                 [roifun,roifunmake,vertexbased]=roistyle2fun(obj.roicoordinates(k).roistyle);
                 roiposnm=obj.roicoordinates(k).Position;
                 roipospix=roiposnm; %later convert
                 roipos1=[];roipos2=[];
                 roipos1(:,1:2)=roipospix(:,1:2);
                 roipos2(:,1)=roipospix(:,1)+obj.images.rangex(2);
                 roipos2(:,2)=roipospix(:,3);

                newroi1(k)=roifunmake(obj.axis,'Position',roipos1);
                addlistener(newroi1(k),'ROIMoved',@(src,evt) obj.updateposition(src,evt,k,1));
                addlistener(newroi1(k),'DeletingROI',@(src,evt) obj.deleteroi(src,evt,k,1));
                if vertexbased
                    addlistener(newroi1(k),'VertexAdded',@(src,evt) obj.updateposition(src,evt,k,1));
                    addlistener(newroi1(k),'VertexDeleted',@(src,evt) obj.updateposition(src,evt,k,1));
                end

                newroi2(k)=roifunmake(obj.axis,'Position',roipos2);
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
            oldpospix=obj.roicoordinates(currentroi).Position; %convert
            [roipospix,roipos]=getpixcoord(obj,src.Position,oldpospix);
            
            newroiposnm=roipospix;
            obj.roicoordinates(currentroi).Position=newroiposnm;
            plotrois(obj)

        end
%         function plotdistances(obj)
%             modality=obj.getSingleGuiParameter('modality').selection;
%             pos=obj.site.evaluation.(obj.name).(modality).Position;
%             
%             period2=mean(diff(pos(:,1)));
%             form='%2.2f';
%             period=(pos(end,1)-pos(1,1))/(size(pos,1)-1);
%             
%             title(obj.axis,['Period: ' num2str(period,form)])
%             
%             obj.site.evaluation.(obj.name).(modality).Period=period;
%             
%         end
%         function fit_callback(obj,a,b)
%             pos=obj.site.evaluation.(obj.name).Position;
%             period=mean(diff(pos(:,1)));
%             
%             
%         end
    end

end

function [roipospix,roipos]=getpixcoord(obj,roiposh,roiposold)
roipospix=roiposold;
 if mean(roiposh(:,1))>0
    roipos=2; %right
    roiposh(:,1)=max(roiposh(:,1),0);
    roipospix(:,1)=roiposh(:,1)-obj.images.rangex(2);
    roipospix(:,3)=roiposh(:,2); %move to left
%     roipospix(:,2)=0;
else
    roipos=1;
    roiposh(:,1)=min(roiposh(:,1),0);
    roipospix(:,1:2)=roiposh(:,1:2);
%     roipospix(:,3)=0;
end
end



function pard=guidef(obj)
pard.addroi.object=struct('Style','pushbutton','String','add','Callback',@obj.addroi_callback);
pard.addroi.position=[1,1];
pard.addroi.Width=1;

pard.roiform.object=struct('Style','popupmenu','String',{{'line','polyline','rectangle','ellipse','polygon','free'}});
pard.roiform.position=[1,2];
pard.roiform.Width=2;

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
                 
    posh=[prz.sr_pos(1) prz.sr_pos(2) prz.sr_size(1)*2 prz.sr_size(2)*2];
    locz=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','PSFxnm','intensity_render','phot','numberInGroup',prz.renderfield.selection},'layer',k,'position',obj.site);
    if length(prz.sr_pos)<3
        prz.sr_pos(3)=0;
    end
    xi=locz.xnm-posh(1);
    yi=locz.ynm-posh(2);
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
dx=(obj.images.rangex(2)-obj.images.rangex(1))/2;
dy=obj.images.rangey(2)-obj.images.rangey(1);
pos=dat.IntersectionPoint;
if pos(1)>0 && pos(1)<dx*0.75 %only right side: side view%
    site.pos(3)=site.pos(3)+pos(2)*1000;
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
end

function displayimage(img,hax)
if isempty(img)
    return
end
 imagesc(img.rangex,img.rangey,img.image,'Parent',hax,'Pickable','none','HitTest','off')

set(hax,'Xlim',double(img.rangex))
set(hax,'Ylim',double(img.rangey))
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


function [roifun,roifunmake,vertexbased]=roistyle2fun(roistyle)
vertexbased=true;
switch roistyle
    case 'rectangle'
        roifun=@drawrectangle;
        roifunmake=@Rectangle;
    case 'ellipse'
        roifun=@mydrawellipse;
    case 'polygon'
        roifun=@drawpolygon;
    case 'free'
        roifun=@drawfreehand;
    case 'polyline'
        roifun=@drawpolyline;  
        roifunmake=@images.roi.Polyline;
    case 'line'
        roifun=@drawline; 
        vertexbased=false;
end
end
function resetview(a,b,obj,site)
site.setlineangle(0,0);
site.annotation.polarangle=0;
site.pos(3)=0;
obj.setPar('se_currentPolarAngle',0)
plot3Dviews(obj)
end
function info(a,b)
text='To rotate in x-y plane: click in the left top view image in the direction you want to point left. \nTo change the z-position, click on the right image, left part, this position will be centered. \nTo change the polar angle, click in the right image, right part. The closer you are to the center, the smaller the change.';
msgbox(sprintf(text));
end