classdef Viewer3DV01<interfaces.DialogProcessor
    % Viewer3DV01 Interactive 3D viewerfor localization data. Same render
    % options as the standard 2D renderer, includes stereoscopic rendering
    % and export of animations.
    properties
        axis
        timer
        currentimage
        commandfig
        recpar={};
        stereopar
        locDataL
        posL=[0,0;0,0];
    end
    methods
        function obj=Viewer3DV01(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.inputParameters=[renderSMAP drawerSMAP displayerSMAP];
            obj.inputParameters{end+1}='sr_roihandle';
            obj.inputParameters{end+1}='linewidth_roi';
            obj.inputParameters{end+1}='layers';
            obj.inputParameters{end+1}='numberOfLayers';
            obj.inputParameters=unique(obj.inputParameters);
             obj.showresults=false;
             obj.guiselector.show=true;
        end
%         function makeGui(obj)
%             makeGui@interfaces.DialogProcessor(obj);
%             h=obj.guihandles;           
%         end
        
        function showpanel_callback(obj,a,b)
            if isempty(obj.commandfig)||~isvalid(obj.commandfig)               
                fp=getParentFigure(obj.handle);
                posfig=fp.Position;
                posfig(1)=posfig(1)+posfig(3);
                posfig(3:4)=[200,200];
                obj.commandfig=figure('MenuBar','none','ToolBar','none','Position',posfig);               
            end
            figure(obj.commandfig);
             h.ptranslation=makepanel([0 0.5 0.5 0.5],'translation','');
            h.protation=makepanel([0.5 0.5 0.5 0.5],'rot (command/strg)','command');
            h.pzoom=makepanel([0.5 0.0 0.5 0.5],'zoom (alt)','alt');

            function h=makepanel(pos,title,modifier)
                h=uipanel('Parent',obj.commandfig,'Units','normalized','Position',pos,'Title',title);
                h.Units='normalized';
                uicontrol('Parent',h,'Units','normalized','Position',[1 1 1 1]/3,'String','0','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','0','Character','0')})
                uicontrol('Parent',h,'Units','normalized','Position',[0 1 1 1]/3,'String','<-','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','leftarrow','Character','')})
                uicontrol('Parent',h,'Units','normalized','Position',[2 1 1 1]/3,'String','->','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','rightarrow','Character','')})
                uicontrol('Parent',h,'Units','normalized','Position',[1 2 1 1]/3,'String','^','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','uparrow','Character','')})
                uicontrol('Parent',h,'Units','normalized','Position',[1 0 1 1]/3,'String','v','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','downarrow','Character','')})   
                uicontrol('Parent',h,'Units','normalized','Position',[0 2 1 1]/3,'String','x (,)','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','comma','Character',',')})
                uicontrol('Parent',h,'Units','normalized','Position',[2 2 1 1]/3,'String','(.)','Callback',{@obj.keypress,struct('Modifier',modifier,'Key','period','Character','.')})                  
            end
        end
        function initGui(obj)

        end
        function out=run(obj,p)
            out=[];
            obj.locDataL=[];
            obj.makelocDatacopy;
            obj.addSynchronization('sr_roiposition',[],[],{@obj.redraw}); 
            if isempty(obj.axis)||isstruct(obj.axis)||~isvalid(obj.axis)
                figure;
                obj.axis=gca;
            end
            
            
            
             set(obj.axis,'NextPlot','replacechildren','PickableParts','all','Units','pixels')
            fig=getParentFigure(obj.axis);

            obj.axis.Units='normalized';
             if strcmp(p.stereo.selection,'goggles')
                 fig.ToolBar='none';
                 fig.MenuBar='none';
                 obj.axis.Position=[0 0 1 1];
                 fig.Color=[0 0 0];
             else
                 obj.axis.Position=[0.05 0.05 .94 .9];
                 fig.ToolBar='figure';
                 fig.MenuBar='figure';
                 fig.Color=[0.94 0.94 0.94];
             end
             figure(fig);
             axis(obj.axis,'tight');
             if p.fillimage
                 axis(obj.axis,'normal');
             else
                axis(obj.axis,'equal');
             end
             axis(obj.axis,'ij');
             set(fig,'WindowKeyPressFcn',{@obj.keypress,[]})
             set(fig,'WindowButtonDownFcn',{@obj.mousebutton,1})
             fig.BusyAction='cancel';
             fig.Interruptible='off';
%              set(obj.axis,'PickableParts','all');
             obj.timer=uint64(0);
             obj.redraw
        end
        
        function makelocDatacopy(obj)

            hroi=obj.getPar('sr_roihandle');
            if ~isvalid(hroi)
                return
            end
            posrx=hroi.getPosition;
            len=sqrt(sum((posrx(2,:)-posrx(1,:)).^2));
            meanpos=mean(posrx,1);
           
            dpos=obj.posL-posrx;
            dpos(1,:)=-dpos(1,:);
            if isempty(obj.locDataL)||any(dpos(:)<0)
                lps=obj.getLayerParameters;
                for k=length(lps):-1:1
                    rfields{k}=lps{k}.renderfield.selection;
                end
                zfield=obj.getSingleGuiParameter('zfield').selection;
                zfielderr=obj.getSingleGuiParameter('zfielderr').selection;
                rfields(end+1:end+2)={zfield, zfielderr};
                rfields=horzcat(unique(rfields),{'xnm','ynm','locprecnm','phot','numberInGroup'});
                lenL=len*2;
                obj.posL=[meanpos(1)-lenL, meanpos(2)-lenL;meanpos(1)+lenL, meanpos(2)+lenL];
                posLnm=obj.posL*1000;
                inx=obj.locData.loc.xnm<posLnm(2,1)&obj.locData.loc.xnm>posLnm(1,1);
                iny=obj.locData.loc.ynm<posLnm(2,2)&obj.locData.loc.ynm>posLnm(1,2);
%                 sum(inx&iny)
                
                 obj.locDataL=obj.locData.copy(rfields,inx&iny);
                 obj.locDataL.loc.znm=single(obj.locDataL.loc.(zfield));
                 obj.locDataL.grouploc.znm=single(obj.locDataL.grouploc.(zfield));
                 obj.locDataL.loc.locprecznm=single(obj.locDataL.loc.(zfielderr));
                 obj.locDataL.grouploc.locprecznm=single(obj.locDataL.grouploc.(zfielderr));
                 obj.locDataL.filter({zfield,zfielderr});
            end
        end

        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function mousebutton(obj,src,callbackdata,calltype)
            
            src.WindowButtonMotionFcn = @motion;
            src.WindowButtonUpFcn = @up;
            disteye=100;
            roih=obj.getPar('sr_roihandle');
            pos=roih.getPosition;
            theta=obj.getSingleGuiParameter('theta');
%             roivec=pos(2,:)-pos(1,:);
%             phi=atan2(roivec(2),roivec(1));
%             vold=sph2cart(phi,theta,1);
             oldpos=obj.axis.CurrentPoint;
     
            oldanglex=atan2(oldpos(1,1),disteye);
            oldangley=atan2(oldpos(1,2),disteye);
            axposx=obj.axis.XLim;axposy=obj.axis.YLim;
            function motion(src,callbackdata)
                 posh=obj.axis.CurrentPoint;
                 if (posh(1,1)<axposx(1)||posh(1,1)>axposx(2) || posh(1,2)<axposy(1)||posh(1,2)>axposy(2))
                     up(src,callbackdata)
                 end
%                 vnew=vold+(posh-oldpos)/disteye;
                
                newanglex=atan2(posh(1,1),disteye);
                newangley=atan2(posh(1,2),disteye);
                dax=newanglex-oldanglex;
                day=newangley-oldangley;
%                 theta=obj.getSingleGuiParameter('theta');
                thetan=theta+day;
                obj.setGuiParameters(struct('theta',thetan));
                posn=rotpos(pos,dax);
                if isvalid(roih)
                roih.setPosition(posn);
                end
            end
            function up(src,callbackdata)
                src.WindowButtonMotionFcn='';
%                 src.WindowButtonUpFcn='';
            end
            
        end
        
        
        
        function keypress(obj,a,d2,data)
            if isempty(data)
                data=d2;
            end
            
           %1.up, 2.down, 3.left, 4.right, 5.back, 6.front, 0.reset
            switch data.Character
                case {'w','W','8',30}
                    dir=1;
                case {'s','S','2',31}
                    dir=2;
                case {'a','A','4',28}
                    dir=3;
                case {'d','D','6',29}
                    dir=4;
                case {',','q','Q','7'}
                    dir=5;
                case {'.','e','E','9'}
                    dir=6;
                case {'0'}
                    dir=0;   
                case 32 %space bar: rotate
                     obj.guihandles.rotateb.Value=~ obj.guihandles.rotateb.Value;
                    if obj.guihandles.rotateb.Value
                        obj.rotate_callback;
                    end
                    return
                otherwise 
                    switch data.Key
                        case 'comma'
                            dir=5;
                        case 'period'
                            dir=6;
                        otherwise 
                            return;
                    end
                    
            end
            
            p=obj.getGuiParameters;
            roih=obj.getPar('sr_roihandle');
            pos=roih.getPosition;
            roivec=pos(2,:)-pos(1,:);
            roivecp(2)=roivec(1);
            roivecp(1)=-roivec(2);
            step=0.1;
            stepl=0.3;
            dphi=pi/16;
            dtheta=pi/16;
            if any(strcmp(data.Modifier,'shift'))
                stepfac=0.2;
            else
                stepfac=1;
            end
%             data.Modifier
            if any(strcmp(data.Modifier,'command'))||any(strcmp(data.Modifier,'control'))
                %rotate
                phi=0;
                theta=p.theta;
                switch dir
                    
                    case 1
                       %tilt up down
                       theta=theta+dtheta*stepfac;
                       if theta>pi
                           theta=theta-2*pi;
                       end
                       if theta<-pi
                           theta=theta+2*pi;
                       end
                       phi=0;
                    case 2
                        phi=0;
                       theta=theta-dtheta*stepfac;
                       if theta>pi
                           theta=theta-2*pi;
                       end
                       
                       if theta<-pi
                           theta=theta+2*pi;
                       end
                    case 3
                       phi=dphi*stepfac;
                    case 4
                        phi=-dphi*stepfac;
                    case 0
                        if strcmp(data.Character,'0')
                        theta=0;
                        end
                end
                 mpos=mean(pos,1);
                [dx,dy]=rotcoord(roivec(1)/2,roivec(2)/2,phi);
                pos(1,1)=mpos(1)-dx;
                pos(2,1)=mpos(1)+dx;
                pos(1,2)=mpos(2)-dy;
                pos(2,2)=mpos(2)+dy;
                obj.setGuiParameters(struct('theta',theta));
                
            elseif any(strcmp(data.Modifier,'alt'))
                %change size
                
                switch dir
                    case 6
                        lw=obj.getPar('linewidth_roi');
                        lw2=lw*(1+step*stepfac);
                        obj.setPar('linewidth_roi',lw2);
                       %tilt up down
                    case 5
                        lw=obj.getPar('linewidth_roi');
                        lw2=lw*(1-step*stepfac);
                        obj.setPar('linewidth_roi',lw2);
                      
                    case 3
                        
                       pos(1,:)=pos(1,:)+roivec/2*step*stepfac;
                       pos(2,:)=pos(2,:)-roivec/2*step*stepfac;
                    case 4
                        pos(1,:)=pos(1,:)-roivec/2*step*stepfac;
                       pos(2,:)=pos(2,:)+roivec/2*step*stepfac;
                    case 1
                        po.zmin=p.zmin-stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax+stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 2
                        po.zmin=p.zmin+stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax-stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                end
            else
                switch dir
                    case 6
                        lw=obj.getPar('linewidth_roi')/1000;
                        pos(1,:)=pos(1,:)+stepl*roivecp./norm(roivecp)*stepfac*lw;
                        pos(2,:)=pos(2,:)+stepl*roivecp./norm(roivecp)*stepfac*lw;
                    case 5
                        lw=obj.getPar('linewidth_roi')/1000;
                        pos(1,:)=pos(1,:)-stepl*roivecp./norm(roivecp)*stepfac*lw;
                        pos(2,:)=pos(2,:)-stepl*roivecp./norm(roivecp)*stepfac*lw;
                    case 3
                        pos(1,:)=pos(1,:)+step*roivec*stepfac;
                        pos(2,:)=pos(2,:)+step*roivec*stepfac;
                    case 4
                        pos(1,:)=pos(1,:)-step*roivec*stepfac;
                        pos(2,:)=pos(2,:)-step*roivec*stepfac;
                    case 1
                        po.zmin=p.zmin+stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax+stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 2
                        po.zmin=p.zmin-stepfac*step*(p.zmax-p.zmin);
                        po.zmax=p.zmax-stepfac*step*(p.zmax-p.zmin);
                        obj.setGuiParameters(po);
                    case 0
                        if strcmp(data.Character,'0')
                        po.zmin=p.zmin-(p.zmax+p.zmin)/2;
                        po.zmax=p.zmax-(p.zmax+p.zmin)/2;
                        obj.setGuiParameters(po);
                        end
                end
            end
            roih.setPosition(pos);
        end
        
        
        function redraw(obj)
                if isempty(obj.axis)||isstruct(obj.axis)||~isvalid(obj.axis)
                    return
                end
            renderfield={};
            roih=obj.getPar('sr_roihandle');
            if ~isa(roih,'imline')
                return
            end
            obj.makelocDatacopy;
            p=obj.getAllParameters;
            stereo=p.stereo.Value>2;
            locCopy=obj.locDataL; %maybe not needed
            lo=logical(obj.getPar('sr_layerson'));
            layerson=find(lo);
            indg=0;indu=0;
            if sum(lo)==0
                return
            end
            nl=obj.getPar('numberOfLayers');
            g=locCopy.isgrouped(1:nl);
            
            gt=g(lo(1:nl));
            group=zeros(2,1);
            if any(gt==1)
                group(2)=1;
            end
            if any(gt==0)
                group(1)=1;    
            end
            
            for k=nl:-1:1
                renderfield{k}=p.(['layer' num2str(k) '_']).renderfield.selection;
            end
            
            zmean=(p.zmin+p.zmax)/2;
            if p.setpixelsize
                ph.sr_pixrec=p.pixrecset;
            else
                ph.sr_pixrec=p.sr_pixrec;
            end           

            ph.rangey=[p.zmin p.zmax];
            
            rpos=roih.getPosition;

            lr=sqrt(sum((rpos(2,:)-rpos(1,:)).^2));
            rx=[-lr/2 lr/2]*1000;
            ax=obj.axis;
            obj.timer=tic;
           
            ph.sr_roihandle=obj.getPar('sr_roihandle');       
            ph.rangex=rx;
%             
%             if group(1)
%                 [loc,indu,sortind]=getlocrot('ungrouped','inlayeru');  
%                
%             end
%             
%             if group(2)
%                 [locg,indg,sortindg]=getlocrot('grouped','inlayerg');
% %                 locg.ballradius=0*locg.xnmline+p.transparencypar(2);
%             end
%            
%             if sum(indg)==0&&sum(indu)==0
%                 return
%             end
            %transparency
            transparency.parameter=p.transparencypar;
            transparency.mode=p.transparencymode.Value;
            switch p.transparencymode.Value
                case 1 %MIP
                case 2 %transparent
                    transparency.parameter=p.transparencypar(1)/ph.sr_pixrec(1)^2*10;
                case 3 %balls
                    transparency.parameter=[p.transparencypar(1)/ph.sr_pixrec(1)^2*10 p.transparencypar(end)];
                    
            end
            
            ph.sr_axes=[];
            for k=p.numberOfLayers:-1:1
                pl=p.(['layer' num2str(k) '_']);
                if pl.layercheck
                    if length(obj.recpar)>=k
                        rp=obj.recpar{k};
                    else
                        rp=[];
                    end
                     pr=copyfields(copyfields(copyfields(p,pl),ph),rp);
                     loc=getlocrot(k,pl); 
                     if stereo
                         pr=getstereosettings(pr,1);
                         layer1(k).images=renderplotlayer(pr,1);
                         %same intensity scaling
                         pr.imax_min=layer1(k).images.finalImages.imax;
                         obj.currentimage.imax(k)=pr.imax_min;
                         pr.imaxtoggle=0;
                         pr=getstereosettings(pr,2);
                         layer2(k).images=renderplotlayer(pr,2);
                     else
                        layer(k).images=renderplotlayer(pr,0);
                        if ~isempty(layer(k).images.finalImages.imax)
                        obj.currentimage.imax(k)=layer(k).images.finalImages.imax;
                        end
                     end
                end
            end
            if stereo
                srim1=displayerSMAP(layer1,pr);
                srim2=displayerSMAP(layer2,pr);
                switch p.stereo.Value
                    case 3
                        srim=srim1;
                        srim.image=srim1.image+srim2.image;
                    case {4,5}
                        srim=assembleSideviews(srim1,srim2,p);
                    case {6}
                        srim=assembleSideviews(srim1,srim2,p);
                end
            else
                srim=displayerSMAP(layer,pr);
            end
            if isempty(srim)
                return
            end
            obj.currentimage=copyfields(obj.currentimage,srim);
            
            him=imagesc(srim.rangex*1000,srim.rangey*1000,srim.image,'Parent',ax);
             ax.HitTest='on';
            ax.PickableParts='all';
            ax.YDir='normal';
            him.PickableParts='none';
            if  p.stereo.Value==5 %goggles
                axis(ax,'off')
                ax.Parent.ToolBar='none';
                ax.Parent.MenuBar='none';
                ax.Units='normalized';
                ax.Position=[0.0 0.0 1 1];
            end
           drawnow limitrate 
           
           
           
            function srim=assembleSideviews(srim1,srim2,p)
                srim=srim1;
                mpar=1/p.stereomagnification;
                widths=obj.axis.Parent.Position(3);
                
                ar=obj.axis.Parent.Position(3)/obj.axis.Parent.Position(4);
                sx=ceil(round(widths*mpar/2))*2;
                sy=ceil(round(widths*mpar/ar/2))*2;
                image=zeros(sy,sx,3);
                s=size(srim1.composite);
                mpx=ceil(sx/2);
                mpy=ceil(sy/2);
                shalfimy=floor(s(1)/2);
                shalfimx=floor(s(2)/2);
                eyepx=round(obj.stereopar.eyemm/p.monitorpmm*mpar);
                mpxh=mpx-eyepx;
                shalfy=min(shalfimy,mpy);
                shalfx=min(min(shalfimx,mpxh),eyepx);
%                 shalfy=100;shalfx=100;
%                 im1=srim1.composite;
                rim1=shalfimy-shalfy+1:shalfimy+shalfy;
                rim2=shalfimx+1-shalfx:shalfimx+shalfx;
                imcut1=srim1.composite(rim1,rim2,:);
                imcut2=srim2.composite(rim1,rim2,:);
                
                image(mpy-shalfy+1:mpy+shalfy,mpxh-shalfx+1:mpxh+shalfx,:)=imcut1;
                mpxh=mpx+eyepx;
                image(mpy-shalfy+1:mpy+shalfy,mpxh-shalfx+1:mpxh+shalfx,:)=imcut2;
%                 image(mpy-shalf(1)+1:mpy-shalf(1)+s(1),mpxh-shalf(2)+1:mpxh-shalf(2)+s(2),:)=imcut2;                
                
                fx=sx/s(2);fy=sy/s(1);
                srim.image=image;
                srim.rangex=srim1.rangex*fx;
                srim.rangey=srim1.rangey*fy;
            end
            function images=renderplotlayer(pr,stereochannel)
                if stereochannel>0
%                     if pr.groupcheck
%                         locg.x=locg.(['x' num2str(stereochannel)]);
%                     else
                        loc.x=loc.(['x' num2str(stereochannel)]);
%                     end
                end
                pr.shiftxy_min=0;
                pr.shiftxy_max=0;
                indin=true(length(loc.x),1);
%                  if pr.groupcheck
%                         ind=find(layerson==k);
%                         indroi=locg.inlayerg{ind};
%                         indh=(indroi(indg));
%                         images.srimage=renderSMAP(locg,pr,k,indh(sortindg),transparency);
%                  else 
%                      ind=find(layerson==k);
%                      indroi=loc.inlayeru{ind};
%                      indh=(indroi(indu));
                     images.srimage=renderSMAP(loc,pr,k,indin,transparency);
%                  end
                images.finalImages=drawerSMAP(images.srimage,pr);        
                
            end
            function [loc,indu,sortind]=getlocrot(layer,pl)

                [loc,indu]=locCopy.getloc({'xnmline','ynmline','znm','locprecnm','locprecznm',renderfield{:},'numberInGroup','phot'},...
                    'position','roi','layer',layer,'shiftxy',[pl.shiftxy_min,pl.shiftxy_max,pl.shiftxy_z]);   
                loc.znm=loc.znm+pl.shiftxy_z;
                if strcmp(p.animatemode.selection,'Translate')&&strcmp(p.raxis.selection,'vertical')
                    thetaoffset=pi/2;
%                     induf=find(indu);
                    indz=(loc.znm>p.zpos-p.zdist/2 & loc.znm<=p.zpos+p.zdist/2);
                    indu(indu)=indz;
                    loc=copystructReduce(loc,indz);
       
                else
                    thetaoffset=0;
                end
                [yrot,depth]=rotcoord(loc.znm-zmean,loc.ynmline,p.theta+thetaoffset);
                [sortdepth,sortind]=sort(-depth);
                sortdepth=depth(sortind);
%                 sortdepth=sortdepth-max(sortdepth); %reference point on plane
                
                    eyemm=35;
                    dplanemm=p.dplanemm;
                    pixmm=p.monitorpmm;
                    obj.axis.Units='pixel';
                    widthpix=obj.axis.Position(3);
%                     pixnm=ph.sr_pixrec;
                    widthnm=ph.rangex(2)-ph.rangex(1);
                    heightpix=obj.axis.Position(4);
                    widthmm=widthpix*pixmm;
                    eyenm=eyemm*widthnm/widthmm;
                    dplanenm=dplanemm*widthnm/widthmm;
                    xe1=eyenm;xe2=-eyenm;
                    x=loc.xnmline(sortind);
                    y=yrot(sortind)+zmean;
                    if length(p.transparencypar)<2
                        p.transparencypar(2)=5;
                    end
               if stereo   
%                     dplanenm=max(abs(sortdepth))*1.5;
                    fac=(sortdepth+dplanenm)/dplanenm;
                    loc.x2=(x-xe1)./(fac)+xe1;
                    loc.x1=(x-xe2)./(fac)+xe2;
                    loc.y=(y)./(fac);

%                     loc.x1=(x-xe1)./(1+sortdepth/dplanenm)+xe1;
%                     loc.x2=(x-xe2)./(1+sortdepth/dplanenm)+xe2;
%                     loc.y=(y)./(1+sortdepth/dplanenm);
                    obj.stereopar=struct('eyemm',eyemm,'eyenm',eyenm,'widthnm',widthnm,'widthmm',widthmm,'widthpix',widthpix,'heightpix',heightpix);
                elseif strcmp(p.stereo.selection,'perspective')
                    loc.x=(x)./(1+sortdepth/dplanemm);
                    loc.y=(y)./(1+sortdepth/dplanemm);
                    maxr=4;
                    minr=0.1;
                    radius=1./(1+sortdepth/dplanemm);
                    radius(radius<0)=0;
                    radius=min(radius,maxr);
                    loc.perspective=radius;
                    loc.ballradius=0*loc.x+p.transparencypar(2);
                    loc.intensity_render=radius.^2;
%                     loc.locprecnm=loc.locprecnm.*radius;
%                     loc.locprecznm=loc.locprecznm.*radius;
                else
                    loc.x=x;
                    loc.y=y;
                    loc.ballradius=0*loc.x+p.transparencypar(2);
                end
    
                
                %change later:
                sx=loc.locprecnm(sortind);
                if ~isempty(loc.locprecznm)
                sy=loc.locprecznm(sortind);
                else
                    sy=sx;
                end
                
                th=p.theta+thetaoffset;
                sxr=sx*sin(th);
                syr=sy*cos(th);
                
                loc.sx=sx;
                loc.sy=sqrt(sxr.^2+syr.^2);
%                 loc.y=yrot(sortind)+zmean;
                loc.znm=loc.znm(sortind);
                loc.numberInGroup=loc.numberInGroup(sortind);
                loc.phot=loc.phot(sortind);
                for kc=1:length(renderfield)
                    if ~isempty(loc.(renderfield{kc}))
                        loc.(renderfield{kc})=loc.(renderfield{kc})(sortind);
                    end
                end
        end

        end
        
        function rotate_callback(obj,button,b,savemovie)
            if nargin<4
                savemovie=[];
            end
            global SMAP_stopnow
            bh=obj.guihandles.rotateb;
            if bh.Value
                bh.FontWeight='bold';
            else
                bh.FontWeight='normal';
            end
            p=obj.getGuiParameters;
            roih=obj.getPar('sr_roihandle');
            
            %initialize movie
            if ~isempty(savemovie)
                indframe=savemovie.frames;
                            s=size(obj.currentimage.image);
                if length(s)==2
                    s(3)=1;
                end
                outim=zeros(s(1),s(2),s(3),savemovie.frames);
                 obj.redraw;
            end
                
           obj.axis.Parent.CurrentCharacter='x';
            while bh.Value && ~SMAP_stopnow && strcmp(p.raxis.selection,obj.getSingleGuiParameter('raxis').selection) && (isempty(savemovie)||indframe>0)
                if obj.axis.Parent.CurrentCharacter==32
                    bh.Value=0;
                end
                if ~isempty(savemovie)
                    outim(:,:,:,savemovie.frames-indframe+1)=obj.currentimage.image;
                    indframe=indframe-1;
                end
                
                switch p.animatemode.selection
                    case 'Rotate'
                        switch p.raxis.selection
                            case 'vertical'
                                    pos=roih.getPosition;
                                    posr=rotpos(pos,obj.getSingleGuiParameter('dangle')*pi/180);
                                    roih.setPosition(posr);           
                            case 'horizontal'
                                    theta=obj.getSingleGuiParameter('theta');
                                    theta=theta-obj.getSingleGuiParameter('dangle')*pi/180;
                                    theta=mod(theta,2*pi);                       
                                    obj.setGuiParameters(struct('theta',theta));
                                    obj.redraw;   
                        end
                    case 'Translate'
                       switch p.raxis.selection
                           
                            case 'vertical'
                                zpos=obj.getSingleGuiParameter('zpos');
                                dz=obj.getSingleGuiParameter('dangle');
                                zrange=obj.getSingleGuiParameter('anglerange');
                                znew=zpos+dz;
                                if znew>zrange(2)
                                    znew=zrange(1);
                                end
                                if znew<zrange(1);
                                    znew=zrange(2);
                                end
                                obj.setGuiParameters(struct('zpos',znew));
                                obj.redraw;
%                                     pos=roih.getPosition;
%                                     roivec=pos(2,:)-pos(1,:);
                                    
%                                     posr=rotpos(pos,obj.getSingleGuiParameter('dangle')*pi/180);
%                                     roih.setPosition(posr);           
                            case 'horizontal'
                                
                                    xpos=obj.getSingleGuiParameter('zpos');
                                    step=obj.getSingleGuiParameter('dangle');
                                    xrange=obj.getSingleGuiParameter('anglerange');
                                    
                                  
                                    if xpos+step>xrange(2)
%                                         xnew=xrange(1);
                                        step=-(xrange(2)-xrange(1));
                                        xnew=xpos+step;
                                    
                                    elseif xpos+step<xrange(1);
%                                         xnew=xrange(2);
                                        step=(xrange(2)-xrange(1));
                                        xnew=xpos+step;
                                    else
                                        xnew=xpos+step;
                                    end
                                    obj.setGuiParameters(struct('zpos',xnew));
                                    
                                    pos=roih.getPosition;
                                    roivec=pos(2,:)-pos(1,:);
                                    roivecp(2)=roivec(1);
                                    roivecp(1)=-roivec(2);
                                    roivec=roivecp/norm(roivecp);
                                    
                                    pos(1,:)=pos(1,:)+step*roivec/1000;
                                    pos(2,:)=pos(2,:)+step*roivec/1000;
%                                     posz=pos+roivec/norm(roivec)*obj.getSingleGuiParameter('dangle');
                                    roih.setPosition(pos);
%                                     theta=obj.getSingleGuiParameter('theta');
%                                     theta=theta-obj.getSingleGuiParameter('dangle')*pi/180;
%                                     theta=mod(theta,2*pi);                       
%                                     obj.setGuiParameters(struct('theta',theta));
%                                     obj.redraw;   
                        end
                end
                
%                 pause(0.01)
            end
            if ~isempty(savemovie)
                options.color=true;
                options.message=true;
                options.comp='lzw';

                imout=uint8(outim*(2^8-1));
                saveastiff(imout,savemovie.file,options)
            end
            
            obj.recpar={};
            if  ~ strcmp(p.raxis.selection,obj.getSingleGuiParameter('raxis').selection)
                rotate_callback(obj,button,b)
            end 
            
        end
        function savemovie_callback(obj,a,b)
            [path,fo]=fileparts(obj.locData.files.file(1).name);
            [file,path]=uiputfile([path filesep fo '.tif']);
            if ~file
                return
            end
            savemovie.file=[path file];
            
            
            p=obj.getGuiParameters(false,true);
            savemovie.frames=ceil((p.anglerange(2)-p.anglerange(1))/p.dangle);
            
            obj.redraw;
            for k=1:length(obj.currentimage.imax)
                obj.recpar{k}.imaxtoggle=false;
                obj.recpar{k}.imax_min=obj.currentimage.imax(k);
            end
            
            %initialize
            switch p.animatemode.selection
                    case 'Rotate'
                        switch p.raxis.selection
                            case 'vertical'          
                            case 'horizontal'
                                theta=p.anglerange(1)*pi/180;                      
                                obj.setGuiParameters(struct('theta',theta));  
                        end
                    case 'Translate'
                       switch p.raxis.selection                           
                            case 'vertical'
                                znew=p.anglerange(1);
                                obj.setGuiParameters(struct('zpos',znew));          
                            case 'horizontal'
                                xnew=p.anglerange(1);
                                obj.setGuiParameters(struct('zpos',xnew));
                                xpos=p.zpos;
                                step=xnew-xpos;
                                roih=obj.getPar('sr_roihandle');
                                pos=roih.getPosition;
                                roivec=pos(2,:)-pos(1,:);
                                roivecp(2)=roivec(1);
                                roivecp(1)=-roivec(2);
                                roivec=roivecp/norm(roivecp);
                                pos(1,:)=pos(1,:)+step*roivec/1000;
                                pos(2,:)=pos(2,:)+step*roivec/1000;
                                roih.setPosition(pos); 
                        end
            end           
            button=obj.guihandles.rotateb;
            button.Value=1;
            obj.rotate_callback(0,0,savemovie)
            obj.recpar={};
            button.Value=0;
            
            
%             global SMAP_stopnow
%             [path,fo]=fileparts(obj.locData.files.file(1).name);
%             [file,path]=uiputfile([path filesep fo '.tif']);
%             if ~file
%                 return
%             end
%             p=obj.getGuiParameters(false,true);
%             if length(p.anglerange)==1
%                 p.anglerange(2)=p.anglerange(1);
%                 p.anglerange(1)=0;
%             end
%             if p.dangle<0
%                 angles=(p.anglerange(2):p.dangle:p.anglerange(1))*pi/180;
%             else
%                 angles=(p.anglerange(1):p.dangle:p.anglerange(2))*pi/180;
%             end
%             s=size(obj.currentimage.image);
%             if length(s)==2
%                 s(3)=1;
%             end
%             outim=zeros(s(1),s(2),s(3),length(angles));
%             obj.redraw;
%             for k=1:length(obj.currentimage.imax)
%                 obj.recpar{k}.imaxtoggle=false;
%                 obj.recpar{k}.imax=obj.currentimage.imax(k);
%             end
% 
%             switch p.raxis.selection
%                 case 'vertical'
%                     roih=obj.getPar('sr_roihandle');
%                     pos=roih.getPosition;
%                     for k=1:length(angles) 
%                         posr=rotpos(pos,angles(k));
%                         roih.setPosition(posr);
%                         outim(:,:,:,k)=obj.currentimage.image;
%                         drawnow
%                         if SMAP_stopnow
%                             break
%                         end
%                     end
% 
%                 case 'horizontal'  
%                     for k=1:length(angles) 
%                         obj.setGuiParameters(struct('theta',angles(k)));
%                         obj.redraw;
%                         outim(:,:,:,k)=obj.currentimage.image;
%                         drawnow
%                         if SMAP_stopnow
%                             break
%                         end
%                     end    
%             end
%             options.color=true;
%             options.message=true;
%             options.comp='lzw';
% 
%             imout=uint8(outim*(2^8-1));
%             saveastiff(imout,[path,file],options)
%             obj.recpar={};
        end
        
        function resetazimuth(obj, a,b)
            obj.setGuiParameters(struct('theta',0));
            obj.redraw;
        end
        
        function savesideview_callback(obj,a,b)
            f=obj.getPar('lastSMLFile');
            fn=strrep(f,'sml.mat','3dxz.tif');
            [file,path]=uiputfile(fn);
            if file
                imwrite(obj.currentimage.image,[path file],'tif')
            end
            
        end
                   
    end
end

function pos=rotpos(pos,angle)
roivec=pos(2,:)-pos(1,:);
mpos=mean(pos,1);
[dx,dy]=rotcoord(roivec(1)/2,roivec(2)/2,angle);
pos(1,1)=mpos(1)-dx;
pos(2,1)=mpos(1)+dx;
pos(1,2)=mpos(2)-dy;
pos(2,2)=mpos(2)+dy;
end

function p=getstereosettings(p,channel)
switch p.stereo.Value
    case 3
        p.render_colormode.selection='normal';
        p.render_colormode.Value=1;
        p.colorfield_max=1;
        p.colorfield_min=0;
        if channel==1
            p.lut.selection='cyan';
            
        else
            p.lut.selection='red';
        end
        
    case 3
    case 4
end
end

function zfield_callback(obj, a,b)
str=obj.guihandles.zfield.String;
ind=find(strcmp(str,'znm'));
if ~isempty(ind)
    obj.guihandles.zfield.Value=ind(1);
end

ind=find(strcmp(str,'locprecznm'));
if ~isempty(ind)
    obj.guihandles.zfielderr.Value=ind(1);
else
    ind=find(strcmp(str,'locprecnm'));
    if ~isempty(ind)
        obj.guihandles.zfielderr.Value=ind(1);
    end
end

end

function pard=guidef(obj)
pard.zfieldt.object=struct('String','z, zerr','Style','text');
pard.zfieldt.position=[1,1];
pard.zfieldt.Width=0.5;
pard.zfieldt.Optional=true;

pard.zfield.object=struct('String',' ','Style','popupmenu');
pard.zfield.position=[1,1.3];
pard.zfield.Width=.7;
pard.zfield.Optional=true;

pard.zfielderr.object=struct('String',' ','Style','popupmenu');
pard.zfielderr.position=[1,1.9];
pard.zfielderr.Width=.75;
pard.zfielderr.Optional=true;

pard.text2.object=struct('String','zmin/zmax','Style','text');
pard.text2.position=[2,1];
pard.text2.Width=0.6;
pard.text2.Optional=false;
% pard.text3.object=struct('String','zmax','Style','text');
% pard.text3.position=[2,2];
% pard.text3.Width=0.3;

pard.setpixelsize.object=struct('String','pixelsize x,z','Style','checkbox','Value',1);
pard.setpixelsize.position=[4,1];
pard.setpixelsize.Width=1;
pard.setpixelsize.Optional=true;
pard.pixrecset.object=struct('Style','edit','String','2 2'); 
pard.pixrecset.position=[4,1.8];
pard.pixrecset.Width=0.4;
pard.pixrecset.Optional=true;

pard.fillimage.object=struct('Style','checkbox','String','fill','Value',0); 
pard.fillimage.position=[4,2.2];
pard.fillimage.Width=0.4;
pard.fillimage.Optional=true;

pard.transparencymode.object=struct('String',{{'projection', 'transparency','balls'}} ,'Style','popupmenu');
pard.transparencymode.position=[6,1];
pard.transparencymode.Width=1.5;
pard.transparencymode.TooltipString=sprintf('maximum intensity, \n partial transparency (parameter is related to transparency), \n render as ball (parameter is ball diamter in reconstructed pixels');
pard.transparencymode.Optional=true;

pard.thetat.object=struct('String','Polar angle theta ','Style','pushbutton','Callback',@obj.resetazimuth);
pard.thetat.position=[3,1];
pard.thetat.Width=1.1;
pard.thetat.TooltipString='Push to set polar angle to zero';
pard.thetat.Optional=false;

pard.zmin.object=struct('Style','edit','String','-400'); 
pard.zmin.position=[2,1.6];
pard.zmin.Width=0.5;
pard.zmin.Optional=false;

pard.zmax.object=struct('Style','edit','String','400'); 
pard.zmax.position=[2,2.1];
pard.zmax.Width=0.5;
pard.zmax.Optional=false;

pard.theta.object=struct('Style','edit','String','0'); 
pard.theta.position=[3,2.1];
pard.theta.Width=0.5;
pard.theta.Optional=false;

pard.transparencypar.object=struct('Style','edit','String','1'); 
pard.transparencypar.position=[6,2.5];
pard.transparencypar.Width=0.5;
pard.transparencypar.TooltipString=pard.transparencymode.TooltipString;
pard.transparencypar.Optional=true;



pard.showcontrols.object=struct('String','Show Controls','Style','pushbutton','Callback',@obj.showpanel_callback);
pard.showcontrols.position=[8,3];
pard.showcontrols.Width=1;
pard.showcontrols.TooltipString='opens control panel to move, rotate, zoom';
pard.showcontrols.Optional=false;

pard.rotateb.object=struct('String','Animate','Style','togglebutton','Callback',@obj.rotate_callback);
pard.rotateb.position=[1,3];
pard.rotateb.TooltipString='Start continuous animation';
pard.rotateb.Optional=false;

pard.animatemode.object=struct('String',{{'Rotate','Translate'}},'Style','popupmenu');
pard.animatemode.position=[1,4];
pard.animatemode.TooltipString='Select rotation or translation';
pard.animatemode.Optional=false;

pard.raxis.object=struct('String',{{'horizontal','vertical'}},'Style','popupmenu');%,'Callback',@obj.axischange_callback);
pard.raxis.position=[2,4];
pard.raxis.TooltipString='axis for continuous rotation / direction of translation';
pard.raxis.Optional=false;

pard.danglet.object=struct('String','step','Style','text');
pard.danglet.position=[3,3];
pard.danglet.Width=0.5;
pard.danglet.TooltipString='step in angle between frames of continuous rotation';
pard.danglet.Optional=true;

pard.dangle.object=struct('String','3','Style','edit');
pard.dangle.position=[3,3.5];
pard.dangle.Width=0.5;
pard.dangle.TooltipString=pard.danglet.TooltipString;
pard.dangle.Optional=true;

pard.zpost.object=struct('String','position','Style','text');
pard.zpost.position=[3,4];
pard.zpost.Width=0.5;
pard.zpost.TooltipString='slice thickness (nm) for z-movie';
pard.zpost.Optional=true;

pard.zpos.object=struct('String','0','Style','edit');
pard.zpos.position=[3,4.5];
pard.zpos.Width=0.5;
pard.zpos.TooltipString=pard.danglet.TooltipString;
pard.zpos.Optional=true;

pard.zdistt.object=struct('String','dz','Style','text');
pard.zdistt.position=[5,4];
pard.zdistt.Width=0.5;
pard.zdistt.TooltipString='slice thickness (nm) for z-movie';
pard.zdistt.Optional=true;

pard.zdist.object=struct('String','10','Style','edit');
pard.zdist.position=[5,4.5];
pard.zdist.Width=0.5;
pard.zdist.TooltipString=pard.danglet.TooltipString;
pard.zdist.Optional=true;

pard.savemovie.object=struct('String','save movie','Style','pushbutton','Callback',@obj.savemovie_callback);
pard.savemovie.position=[2,3];
pard.savemovie.TooltipString='Save rotating movie. Uses min - max angle';
pard.savemovie.Optional=false;

pard.tx.object=struct('String','min max','Style','text');
pard.tx.position=[4,3];
pard.tx.TooltipString='start and stop angle/position. Uses step from above';
pard.tx.Optional=true;

pard.anglerange.object=struct('String','0 360','Style','edit');
pard.anglerange.position=[4,4];
pard.anglerange.TooltipString=pard.tx.TooltipString;
pard.anglerange.Optional=true;

pard.stereo.object=struct('Style','popupmenu','String',{{'no stereo','perspective','anaglyphs (color)','free-view','cross-eyed','goggles'}}); 
pard.stereo.position=[7,1];
pard.stereo.Width=1;
pard.stereo.TooltipString='mode for stereo reconstruction. anaglyphs need cyan/red glasses.';
pard.stereo.Optional=true;

pard.tx2.object=struct('String','pixel (mm)','Style','text');
pard.tx2.position=[7,2];
pard.tx2.Width=.6;
pard.tx2.Optional=true;

pard.monitorpmm.object=struct('Style','edit','String','0.12'); 
pard.monitorpmm.position=[7,2.6];
pard.monitorpmm.Width=.4;
pard.monitorpmm.TooltipString='Size of a monitor pixel in nm. Used to calculate eye distance. Adjust freely to optimize 3D reconstruction. iPhone: 0.12; iMac: 0.3.';
pard.monitorpmm.Optional=true;
pard.tx2.TooltipString=pard.monitorpmm.TooltipString;

pard.tx3.object=struct('String','d plane','Style','text');
pard.tx3.position=[8,1];
pard.tx3.Width=.6;
pard.tx3.Optional=true;
pard.dplanemm.object=struct('Style','edit','String','1000');  
pard.dplanemm.position=[8,1.6];
pard.dplanemm.Width=.4;
pard.dplanemm.TooltipString='Distance to plane of reconstruction. Smaller values result in stronger 3D effect';
pard.dplanemm.Optional=true;
pard.tx3.TooltipString=pard.dplanemm.TooltipString;

pard.tx4.object=struct('String','Mag','Style','text');
pard.tx4.position=[8,2];
pard.tx4.Width=.6;
pard.tx4.Optional=true;
pard.stereomagnification.object=struct('Style','edit','String','2');  
pard.stereomagnification.position=[8,2.6];
pard.stereomagnification.Width=.4;
pard.stereomagnification.TooltipString='Used only for side-by-side reconstructions of left and right image. You can adjust the magnification with this paramter';
pard.stereomagnification.Optional=true;
pard.tx4.TooltipString=pard.stereomagnification.TooltipString;


pard.savesideview.object=struct('String','save tif','Style','pushbutton','Callback',@obj.savesideview_callback);
pard.savesideview.position=[8,4];
pard.savesideview.Width=1;
pard.savesideview.Optional=false;


pard.syncParameters={{'locFields','zfield',{'String'},{@zfield_callback,obj}},{'locFields','zfielderr',{'String'}}};
 
pard.plugininfo.name='Viewer 3D';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description=sprintf(['localization based 3D viewer for superresolution data.\n'...
    'Linked to linear Roi. Updates automatically. \n Keyboard shortcuts: '...
    '\n left (<-, a,4), right (->,d,6), up (uparrow,w,8),down (downarrow,s,2),front (comma,q,7),'...
    'back (period,e,9), reset (0)\n '...
    'translate: no modifier, rotate: command/strg, zoom: alt']);
end