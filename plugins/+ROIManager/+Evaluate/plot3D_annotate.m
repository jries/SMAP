classdef plot3D_annotate<interfaces.SEEvaluationProcessor
%     Calculates the line-profile along user-defined direction and fits it
%     with a Gaussian or double-Gaussian model
    properties
        roicoordinates
        roihandles
        currentroi
        axis
    end
    methods
        function obj=plot3D_annotate(varargin)        
                obj@interfaces.SEEvaluationProcessor(varargin{:});
        end
        function out=run(obj, p)
            make3Daxis(obj) %make axis,
            
            plot3Dviews(obj,p)
            %plot image, add callbacks
            %plot ROIs: in same function
            %function:
               %plotimage
               %plotrois
               %newroi_callback
               

%         if isfield(obj.site.evaluation,obj.name) && isfield(obj.site.evaluation.(obj.name),modality) && ~isempty(obj.site.evaluation.(obj.name).(modality).Position)
%             obj.roihandle=images.roi.Polyline(obj.axis,'Position',obj.site.evaluation.(obj.name).(modality).Position);
%             out=obj.site.evaluation.(obj.name);
%             addlistener(obj.roihandle,'ROIMoved',@(src,evt) obj.updateposition(src,evt));
%             addlistener(obj.roihandle,'VertexAdded',@(src,evt) obj.updateposition(src,evt));
%             addlistener(obj.roihandle,'VertexDeleted',@(src,evt) obj.updateposition(src,evt));
%             obj.plotdistances;
%         end

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
            
            title(obj.axis,['Period: ' num2str(period,form)])
            
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



function plot3Dviews(obj,p)
site=obj.site;
angle=pos2angle(site.annotation.rotationpos.pos);
p1.rotationanglez=angle;
if isfield(site.annotation,'polarangle')
    p1.polarangle=site.annotation.polarangle;
else
    p1.polarangle=0;
end
% p1.sr_size(3)=obj.locData.getPar('se_dz')/2;


imz=make3Dimages(obj,p1);
imzc.image=horzcat(imz{2}.image(end:-1:1,:,:), imz{1}.image); 
imzc.rangey=imz{1}.rangey;
imzc.rangex=[ imz{2}.rangex(1)-imz{2}.rangex(2) imz{1}.rangex(2)-imz{1}.rangex(1)];

displayimage(imzc,obj.axis);
fl='%2.0f';
title(obj.axis,['\theta=' num2str(p1.polarangle,fl) '\circ, \rho=' num2str(p1.rotationanglez,fl) '\circ, z= ' num2str(p1.sr_pos(3),fl) ' nm'])
end         
                
function imagez=make3Dimages(obj,p)                
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
     posh=[prz.sr_pos(1) prz.sr_pos(2) prz.sr_size(1)*2 prz.sr_size(2)*2];
%     locz=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','PSFxnm','intensity_render','phot','numberInGroup',pr.renderfield.selection},'layer',k,'position',posh);
    locz=obj.locData.getloc({'xnm','ynm','znm','locprecnm','locprecznm','PSFxnm','intensity_render','phot','numberInGroup',prz.renderfield.selection},'layer',k,'position',obj.site);
    %     if strcmpi('tiff', pr.rendermode.selection)||strcmpi('raw', pr.rendermode.selection)
%         rawimage=renderSMAP(obj.locData,pr,k);
%     else
%     rawimage=renderSMAP(locz,pr,k);
%     end
%     
%     prz=pr;
    if length(prz.sr_pos)<3
        prz.sr_pos(3)=0;
    end
    xi=locz.xnm-posh(1);
    yi=locz.ynm-posh(2);
    zi=locz.znm-prz.sr_pos(3);
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
imagez{1}=displayerSMAP(layersz,prz);
imagez{2}=displayerSMAP(layerszxy,prz);  
    
%           plotz=false;
%           imax=zeros(numlayers,1);
           
%                 pr=obj.getLayerParameters(k, obj.processors.renderer.inputParameters);   
%                 pd=obj.getLayerParameters(k, obj.processors.drawer.inputParameters); 
%                 pr=copyfields(pr,p);pd=copyfields(pd,p);

                
%                 if nargin>3&&~isempty(pl)&&length(pl)>=k
%                     pr=copyfields(pr,pl{k});
%                 end
%                 if pr.layercheck
%                     pr.ch_filelist.Value=fileind;
%                     pr.ch_filelist.selection=pr.ch_filelist.String{fileind};

%                     obj.processors.renderer.setParameters(pr)
%                     obj.processors.drawer.setParameters(pr);
                
                %filter filenumber
%                     groupc=pr.groupcheck;
%                     filterold=obj.locData.getFilter(k);
% %                     filternew=filterold;
% %                     locs=obj.locData.getloc({'filenumber','xnm','ynm'},'grouping',groupc);
%                     obj.locData.filter('filenumber',k,'inlist',filenumber)
%                     obj.locData.filter('xnm',k,'minmax',[p.sr_pos(1)-p.sr_size(1),p.sr_pos(1)+p.sr_size(1)])
%                     obj.locData.filter('ynm',k,'minmax',[p.sr_pos(2)-p.sr_size(2),p.sr_pos(2)+p.sr_size(2)])
                    
                
%                     filternew.filenumber=(locs.filenumber==filenumber);
%                     filternew.xnm=rec.LocalizationFilter.minMaxFilter(locs.xnm,p.sr_pos(1)-p.sr_size(1),p.sr_pos(1)+p.sr_size(1));
%                     filternew.ynm=rec.LocalizationFilter.minMaxFilter(locs.ynm,p.sr_pos(2)-p.sr_size(2),p.sr_pos(2)+p.sr_size(2));
%                 filternew=myrmfield(filternew,'xnm');
%                 filternew=myrmfield(filternew,'ynm');
%                 obj.locData.setFilter(filternew,k)
                
                    

                    

%                 end
%            end
%            pd=obj.getLayerParameters(k, obj.processors.displayer.inputParameters); 
%            pd=copyfields(pd,p);

%            image=displayerSMAP(layers,pr);
%            image.parameters=pr;
%            image.parameters.layerparameters=players;
%            image.layers=layers;
%            image.imax=imax;
           
%            if plotz

%                imagez=vertcat(imagezxy,imagez);
%            end
%            axis equal
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
site=obj.SE.currentsite;
dx=site.image.rangex(2)-site.image.rangex(1);
dy=site.image.rangey(2)-site.image.rangey(1);
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
    
site.image=[];
plotsite(obj,site)


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

function resetview(a,b,obj,site)
site.setlineangle(0,0);
site.annotation.polarangle=0;
site.pos(3)=0;
obj.setPar('se_currentPolarAngle',0)
redrawsite_callback(a,b,obj)
end
function info(a,b)
text='To rotate in x-y plane: click in the left top view image in the direction you want to point left. \nTo change the z-position, click on the right image, left part, this position will be centered. \nTo change the polar angle, click in the right image, right part. The closer you are to the center, the smaller the change.';
msgbox(sprintf(text));


end