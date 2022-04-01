classdef plot3D_annotate<interfaces.SEEvaluationProcessor
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
        function obj=plot3D_annotate(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            out=[];
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