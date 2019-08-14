classdef GuiFormat<interfaces.GuiModuleInterface & interfaces.LocDataInterface
    properties 
        roihandle
        roiposition
        roicallbackid
        ovboxhandle
        sraxis
        timer
    end
    methods
        function obj=GuiFormat(varargin)
            obj@interfaces.GuiModuleInterface(varargin{:})     
            obj.outputParameters={'sr_layersseparate','sr_plotlayernames'};
        end
        function makeGui(obj)
            widtht=105;
            width=102;
            fontsize=obj.guiPar.fontsize;
            fieldheight=16;       
            h.hroi=uipanel('Parent',obj.handle,'Title','ROI','Units','pixel','Position',[1 2.5*fieldheight widtht 5.7*fieldheight]);
%               h.lwtxt=uicontrol('Parent',h.hroi,'Style','text','String','|<->|','Position',[0, 0,width/3,fieldheight*1.2],'FontSize',fontsize*.85);
            
            swidth=widtht/4-1;
             h.linewidth_roi=uicontrol('Style','edit','Parent',h.hroi,'String','150','Position',[2*swidth,3.3*fieldheight,swidth*2,fieldheight*1.5],'FontSize',fontsize,'Callback',{@lw_callback,obj});
%          
            h.roishow=uicontrol('Parent',h.hroi,'Style','checkbox','String','show','FontSize',fontsize,'Position',[0*swidth,0.3*fieldheight+2,width-fieldheight,fieldheight*1.2],'Callback',{@roi_callback,obj,0},'Value',1);
            h.roidelete=uicontrol('Parent',h.hroi,'Style','pushbutton','String','X','FontSize',fontsize,'Position',[width-fieldheight,0.3*fieldheight+2,fieldheight,fieldheight*1.2],'Callback',{@roi_callback,obj,-1});
           
            
%             icon=rand(round(fieldheight*1.5),round(fieldheight*1.5),3);
            h.roi5=uicontrol('Parent',h.hroi,'Style','pushbutton','String','','FontSize',fontsize,'Position',[0,3.3*fieldheight,swidth,fieldheight*1.5],'Callback',{@roi_callback,obj,5});
            h.roi5.CData=imread('+gui/icons/point.tif');            
            h.roi4=uicontrol('Parent',h.hroi,'Style','pushbutton','String','','FontSize',fontsize,'Position',[swidth,3.3*fieldheight,swidth,fieldheight*1.5],'Callback',{@roi_callback,obj,4});
            h.roi4.CData=imread('+gui/icons/line.tif');
            
            h.roi1=uicontrol('Parent',h.hroi,'Style','pushbutton','String','','FontSize',fontsize,'Position',[0*swidth,1.8*fieldheight,swidth,fieldheight*1.5],'Callback',{@roi_callback,obj,1});
            h.roi1.CData=imread('+gui/icons/rec.tif');
            h.roi2=uicontrol('Parent',h.hroi,'Style','pushbutton','String','','FontSize',fontsize,'Position',[swidth,1.8*fieldheight,swidth,fieldheight*1.5],'Callback',{@roi_callback,obj,2});
            h.roi2.CData=imread('+gui/icons/ellipse.tif');
            h.roi3=uicontrol('Parent',h.hroi,'Style','pushbutton','String','','FontSize',fontsize,'Position',[2*swidth,1.8*fieldheight,swidth,fieldheight*1.5],'Callback',{@roi_callback,obj,3});
            h.roi3.CData=imread('+gui/icons/free.tif');


            h.roi6=uicontrol('Parent',h.hroi,'Style','pushbutton','String','','FontSize',fontsize,'Position',[3*swidth,1.8*fieldheight,swidth,fieldheight*1.5],'Callback',{@roi_callback,obj,6});
            h.roi6.CData=imread('+gui/icons/polygon.tif');
            %format
            h.hformat=uipanel('Parent',obj.handle,'Title','format','Units','pixel','Position',[1 8.5*fieldheight+5 widtht 6*fieldheight]);
             h.picrt=uicontrol('Parent',h.hformat,'Style','text','String','Pixrec (nm)','Position',[0 4.1*fieldheight,width,fieldheight*1.2]);
            
            h.hplus=uicontrol('Parent',h.hformat,'Style','togglebutton','String','+','FontSize',fontsize*1.5,'Position',[0,2.8*fieldheight,width/3,fieldheight*1.4],'Callback',{@format_callback,obj,1});
            h.hminus=uicontrol('Parent',h.hformat,'Style','togglebutton','String','-','FontSize',fontsize*1.5,'Position',[width/3*2,2.8*fieldheight,width/3,fieldheight*1.4],'Callback',{@format_callback,obj,2});
            h.pixrec=uicontrol('Style','edit','Parent',h.hformat,'String','20','Position',[width/3,3*fieldheight,width/3,fieldheight*1.2],'FontSize',fontsize,'Callback',{@format_callback,obj,3});         
            h.pref1=uicontrol('Style','pushbutton','Parent',h.hformat,'String','2','Position',[0,1.6*fieldheight,width/2,fieldheight*1.2],'FontSize',fontsize,'Callback',{@predef_callback,obj},'KeyPressFcn',{@predefkey_callback,obj});
            h.pref2=uicontrol('Style','pushbutton','Parent',h.hformat,'String','10','Position',[width/2,1.6*fieldheight,width/2,fieldheight*1.2],'FontSize',fontsize,'Callback',{@predef_callback,obj},'KeyPressFcn',{@predefkey_callback,obj});
            h.resetview=uicontrol('Style','pushbutton','Parent',h.hformat,'String','Reset','Position',[width/2,0.1*fieldheight,width/2,fieldheight*1.4],'FontSize',fontsize,'Callback',{@resetview_callback,obj});
            h.parformat=uicontrol('Style','pushbutton','Parent',h.hformat,'String','Par','Position',[0,0.1*fieldheight,width/2,fieldheight*1.4],'FontSize',fontsize,'Callback',{@formatpardialog,obj});
            
            h.hlayers=uipanel('Parent',obj.handle,'Title','layers','Units','pixel','Position',[1 15*fieldheight widtht 3*fieldheight+15]);
            h.layeron3=uicontrol('Style','checkbox','Parent',h.hlayers,'String','3','Position',[width/3,fieldheight,width/3,fieldheight],'FontSize',fontsize);
            h.layeron6=uicontrol('Style','checkbox','Parent',h.hlayers,'String','6','Position',[width/3*2,0,width/3,fieldheight],'FontSize',fontsize);
            h.layeron2=uicontrol('Style','checkbox','Parent',h.hlayers,'String','2','Position',[0,0,width/3,fieldheight],'FontSize',fontsize);
            h.layeron5=uicontrol('Style','checkbox','Parent',h.hlayers,'String','5','Position',[width/3*2,fieldheight,width/3,fieldheight],'FontSize',fontsize);
            h.layeron1=uicontrol('Style','checkbox','Parent',h.hlayers,'String','1','Value',1,'Position',[0,fieldheight,width/3,fieldheight],'FontSize',fontsize);
            h.layeron4=uicontrol('Style','checkbox','Parent',h.hlayers,'String','4','Position',[width/3,0,width/3,fieldheight],'FontSize',fontsize);
            h.sr_layersseparate=uicontrol('Style','checkbox','Parent',h.hlayers,'String','split','Value',0,'Position',[0,2*fieldheight,width*.6,fieldheight],'FontSize',fontsize*.8);
            h.sr_plotlayernames=uicontrol('Style','checkbox','Parent',h.hlayers,'String','label','Value',0,'Position',[width/2,2*fieldheight,width*.6,fieldheight],'FontSize',fontsize*0.8);
            
            
             h.layeron1.Callback={@obj.layer_update,1};
            h.layeron2.Callback={@obj.layer_update,2};
            h.layeron3.Callback={@obj.layer_update,3};
            h.layeron4.Callback={@obj.layer_update,4};
            h.layeron5.Callback={@obj.layer_update,5};
            h.layeron6.Callback={@obj.layer_update,6};
            
            callobj=obj;
            obj.addSynchronization('layer1_layercheck',h.layeron1,'Value',{@callobj.layer_update,1})
            obj.addSynchronization('layer2_layercheck',h.layeron2,'Value',{@callobj.layer_update,2})
            obj.addSynchronization('layer3_layercheck',h.layeron3,'Value',{@callobj.layer_update,3})
            obj.addSynchronization('layer4_layercheck',h.layeron4,'Value',{@callobj.layer_update,4})
            obj.addSynchronization('layer5_layercheck',h.layeron5,'Value',{@callobj.layer_update,5})
            obj.addSynchronization('layer6_layercheck',h.layeron6,'Value',{@callobj.layer_update,6})
            obj.addSynchronization('sr_pixrec',h.pixrec,'String',{@callobj.pixrec_callback,3})
            obj.addSynchronization('linewidth_roi',h.linewidth_roi,'String')
           
            
            %overview axes
            hrecg=obj.getPar('mainGuihandle');
            h.overviewimage=uipanel('Parent',hrecg,'Units','pixel','Position',[0 400 435 350],'Tag','OV');
            
            h.ov_axes=axes('Parent',h.overviewimage,'Units','normalized','Position',[0.05 0.07 .95 .93]);
            h.ov_axes.FontSize=obj.guiPar.fontsize-3;
            set(h.ov_axes,'NextPlot','replacechildren','PickableParts','all')
            set(h.ov_axes,'ButtonDownFcn',{@clickOnSrImage,obj})
            obj.setPar('ov_axes',h.ov_axes);
            h.redrawov=uicontrol('Style','pushbutton','Parent',h.overviewimage,'Position',[ 2 1 45 18],'String','update','Callback',{@redrawov_callback,obj});
            h.detachov=uicontrol('Style','pushbutton','Parent',h.overviewimage,'Position',[ 404 1 25 18],'String','~>','Callback',{@detach_callback,obj,h.overviewimage});
            h.redrawov.Units='normalized';
            h.detachov.Units='normalized';
            h.detachov.TooltipString='detach overview image';
            %viewtogglebutton
            h.overview_select=uicontrol('Style','togglebutton','Parent',obj.handle,'String','OV -> filter','Position',[1,1,widtht,fieldheight*1.8],...
                'FontSize',fontsize*1.2,'Callback',{@viewselect_callback,obj});            
            obj.guihandles=h;

%             %menu to detach Format Gui
            f=getParentFigure(obj.handle);
            c2=uicontextmenu(f); 
            obj.handle.UIContextMenu=c2;
             m5 = uimenu(c2,'Label','detach','Callback',{@detach_callback,obj,obj.handle});
           
            h.hplus.TooltipString='zoom out, increase pixelsize';
            h.hminus.TooltipString='zoom in, decrease pixelsize';
            h.pixrec.TooltipString='pixel size for reconstruction (nm)';
            h.pref1.TooltipString=sprintf('preset pixel size. To change this value:  \n 1. click on button and leave mouse over it \n 2. type pixelsize in nm, finishe with Enter');
            h.pref2.TooltipString=h.pref1.TooltipString;
            h.resetview.TooltipString='adjust pixelsize to fit all width or height';
            h.parformat.TooltipString='additional global parameters for rendering';
            h.redrawov.TooltipString='redraw overview image';
            h.overview_select.TooltipString='toggels between overview image and histogram view';
            h.linewidth_roi.TooltipString='width of the roi when using the line or point';
            h.roi4.TooltipString='line roi, width set to the right';
            h.roi5.TooltipString='point roi, rectangular roi around with width and height set to the right';
            h.roi1.TooltipString='rectangular roi';
            h.roi2.TooltipString='circular roi';
            h.roi6.TooltipString='polynome roi';
            h.roi3.TooltipString='free roi';
            h.roishow.TooltipString='redraw roi after updating image. Untick if Roi is in the way.';
            h.roidelete.TooltipString='Delete current ROI';
            makeGui@interfaces.GuiModuleInterface(obj);
            obj.initGui;
        end
        
        function initGui(obj)
            callobj=obj;
            obj.addSynchronization('currentfileinfo',[],[],{@callobj.loaded_notify})
            obj.addSynchronization('sr_figurenumber',[],[],{@callobj.makesrfigure})
            obj.setPar('sr_pos',[0 0]);
            obj.setPar('sr_imsizecheck',false);
            obj.setPar('sr_pixfactor',1);
            obj.setPar('sr_imagesize',1000);
            obj.setPar('sr_layerson',1);
            obj.setPar('sr_colorbarthickness',4);
            addlistener(obj.P,'sr_render',@obj.plotovbox);
            lw_callback(0,0,obj)
            obj.updateFormatParameters;
            delete(obj.getPar('sr_figurehandle'))
            obj.timer=timer('StartDelay',.2);
            stop(obj.timer);
            obj.timer.TimerFcn=@obj.rerender;
     
        end

        function loaded_notify(obj,~,~)
            if isempty(obj.locData.loc),return;end
            updateFormatParameters(obj) 
            redrawov_callback(0,0,obj) 
        end
        function updateFormatParameters(obj)
            ax=obj.getPar('sr_axes');
            if (isempty(ax)||~isvalid(ax))
                hfg=figure;
                obj.setPar('sr_figurenumber',hfg.Number);
                obj.makesrfigure;
            end
            pixrec=obj.getPar('sr_pixrec');
            if obj.getPar('sr_imsizecheck')
                ims=obj.getPar('sr_imagesize');
                if length(ims)==1
                    ims=[ims ims];
                end             
            else
                ax=obj.getPar('sr_axes');
                pos=ax.Position;
                ims=pos(3:4)/obj.getPar('sr_pixfactor'); 
                ims=round(ims-0.001);
            end
            obj.setPar('sr_sizeRecPix',ims);
            obj.setPar('sr_size',(ims/2*pixrec));
            p=obj.getGuiParameters;
            lon=[p.layeron1,p.layeron2,p.layeron3,p.layeron4,p.layeron5,p.layeron6];
            loo=obj.getPar('sr_layerson');
            loo(1:6)=lon;
            obj.setPar('sr_layerson',loo);   
        end

        function sr_axes=makesrfigure(obj,fignumber)
            if nargin<2           
                fignumber=obj.getPar('sr_figurenumber');
                if isempty(fignumber)
                    ff=figure;
                    fignumber=ff.Number;
                    obj.setPar('sr_figurenumber',fignumber)
                end
                    
            end
            clf(fignumber);
            hf=figure(fignumber);

            pos=obj.getPar('mainGuihandle').Position;
            scrsz = get(groot,'ScreenSize');  
            posim=abs([pos(3)+pos(1)+10 pos(2) max(1,min(pos(4),scrsz(3)-pos(3)-pos(1)-30)), max(1,pos(4)-30)]);
            if posim(3)>200 && posim(4)>200 %too small
                set(hf,'Units','pixels','Position', posim);
            end
            hg.hsr=hf;

            hg.sr_axes=axes('Parent',hg.hsr);%,'Units','normalized','Position',[0.10 0.09,.7,.83]);
            set(hg.sr_axes,'NextPlot','replacechildren','PickableParts','all','Units','pixels')
            set(hg.hsr,'WindowButtonDownFcn',{@clickOnSrImageW,obj})
%             set(hg.sr_axes,'ButtonDownFcn',{@clickOnSrImage,obj})
            obj.setPar('sr_figurehandle',hg.hsr);
            set(hg.hsr,'SizeChangedFcn',{@srSizeChange,obj})
            sr_axes=hg.sr_axes;
            hg.hsr.WindowScrollWheelFcn={@scroll_wheel,obj};
            hg.hsr.BusyAction='cancel';
            hg.hsr.ButtonDownFcn=@obj.bringGuiToFront;   
            obj.setPar('sr_axes',sr_axes); 
            obj.sraxis=sr_axes;
            

        end
        
        function bringGuiToFront(obj,a,b)
            figure(obj.getPar('mainGuihandle'));
        end
        
        function plotovbox(obj,a,b)
            pos=obj.getPar('sr_pos');
            size=obj.getPar('sr_size');
            x(1)=(pos(1)-size(1))/1000;
            x(2)=(pos(1)+size(1))/1000;
            y(1)=(pos(2)-size(2))/1000;
            y(2)=(pos(2)+size(2))/1000;

            oax=obj.getPar('ov_axes');        
            oax.NextPlot='add';
            if ishandle(obj.ovboxhandle)
                delete(obj.ovboxhandle)
            end
            obj.ovboxhandle=plot([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'Parent',oax,'Color',[0 1 1],'Pickable','none');
            oax.NextPlot='replacechildren';
        end
    
        function layer_update(obj,layer,a,layer2)   
            if nargin>2
                layer=layer2;
            end
            state=obj.getSingleGuiParameter(['layeron' num2str(layer)]);
            layerson=obj.getPar('sr_layerson');
            layerson(layer)=state;
            obj.setPar('sr_layerson',layerson);
            notify(obj.P,'sr_display')
        end
        
        function pixrec_callback(obj,par)
            obj.updateFormatParameters;
            notify(obj.P,'sr_render')
        end
        
        function pout=roiset(obj,p)
            if nargin<2 %get current roi
                pout.isvalid=~isempty(obj.roihandle)&&isvalid(obj.roihandle);
                if pout.isvalid
                    pout.roimode=class(obj.roihandle);
                    pout.position=obj.roihandle.getPosition;
                end
            else
                pout=[];
                if p.isvalid
                    roi_callback(0,0,obj,p.roimode,p.position)
                end
            end
        end
        
        function linecallback(obj,pos)
            global  roimodecallback 
            persistent lineh
            ax=obj.getPar('sr_axes');
            f=ax.Parent;
%             oldbdf=f.WindowButtonDownFcn;
%             f.WindowButtonDownFcn='';
            if isempty(lineh) || ~isvalid(lineh.text)
                lineh.handle1=line([0 0],[0 0],'Parent',ax);
                lineh.handle2=line([0 0],[0 0],'Parent',ax);
                lineh.ax=ax;
                lineh.text=text(0,0,' ','Color','w','FontSize',16,'VerticalAlignment','bottom','HitTest','off','Parent',ax);
            end  
            switch roimodecallback
                case {4,'imline'}
                    lw=obj.getPar('linewidth_roi');
                    len=sqrt((pos(2,1)-pos(1,1))^2+(pos(2,2)-pos(1,2))^2);                   
                    mpos=mean(pos,1);
                    set(lineh.text,'String',[ num2str(len*1000,'%4.0f') ' nm'],'Position',mpos);
                    roivec=pos(2,:)-pos(1,:);
                    roivecp(2)=roivec(1);
                    roivecp(1)=-roivec(2);
                    roivecp=roivecp/norm(roivecp);
                    posn=0*pos;
                    posn(1,:)=pos(1,:)+roivecp*lw/2000;
                    posn(2,:)=pos(2,:)+roivecp*lw/2000;
                     posn2(1,:)=pos(1,:)-roivecp*lw/2000;
                    posn2(2,:)=pos(2,:)-roivecp*lw/2000;
                    set(lineh.handle1,'XData',posn(:,1),'YData',posn(:,2))
                    set(lineh.handle2,'XData',posn2(:,1),'YData',posn2(:,2))
                case {1,2,'imrect','imellipse'}
                w=pos(3);
                l=pos(4);
%                 linetexth=text(pos(1),pos(2)+pos(4),[ num2str(w*1000,'%4.0f') ' x ' num2str(l*1000,'%4.0f') ' nm'],'Color','w','FontSize',16,'VerticalAlignment','bottom','HitTest','off','Parent',ax);
                set(lineh.text,'Position',[pos(1),pos(2)+pos(4)],'String',[ num2str(w*1000,'%4.0f') ' x ' num2str(l*1000,'%4.0f') ' nm'],'Color','w','FontSize',16,'VerticalAlignment','bottom','HitTest','off','Parent',ax);
                case {3,6,'imfreehand','impoly'} %free
                    w=max(pos(:,1))-min(pos(:,1));
                    l=max(pos(:,2))-min(pos(:,2));
                    set(lineh.text,'Position',[min(pos(:,1)),min(pos(:,2))],'String',[ num2str(w*1000,'%4.0f') ' x ' num2str(l*1000,'%4.0f') ' nm'],'Color','w','FontSize',16,'VerticalAlignment','bottom','HitTest','off','Parent',ax);
%                     linetexth=text(min(pos(:,1)),min(pos(:,2)),[ num2str(w*1000,'%4.0f') ' x ' num2str(l*1000,'%4.0f') ' nm'],'Color','w','FontSize',16,'VerticalAlignment','bottom','HitTest','off','Parent',ax);   
            end
            obj.setPar('sr_roiposition',pos);
            obj.roiposition=pos;
%             f.WindowButtonDownFcn=oldbdf;
        end
        function rerender(obj,varargin)
            notify(obj.P,'sr_render')
        end
    end
end

function scroll_wheel(a,eventdata,obj)
stop(obj.timer)
persistent timercount 
% if isempty(totalscroll)
%     totalscroll=0;
% end
mint=0.01;
% totalscroll=1+totalscroll;
if isempty(timercount)||toc(timercount)>mint
    vs=eventdata.VerticalScrollCount;
    if vs>0
        eventcase=1;
    else
        eventcase=2;
    end
    obj.setPar('fastrender',true)
    format_callback(0,0,obj,eventcase)
    obj.setPar('fastrender',false)
%     totalscroll=0;
end
timercount=tic;
start(obj.timer);
end

function format_callback(handle,action,obj,eventcase)
% if nargin<5
%     totalscroll=1;
% end
h=obj.guihandles;
pixrec=str2num(get(h.pixrec,'String'));
switch eventcase
    case 1
        set(h.hminus,'Value',0)
        set(h.hplus,'Value',1)
        pixrec=pixrec*2;
        
    case 2
        set(h.hminus,'Value',1)
        set(h.hplus,'Value',0)
        pixrec=pixrec/2;
    case 3
        
        disp('pixelsize changed')       
end
pixrec=round(pixrec,2,'significant');
obj.setPar('sr_pixrec',num2str(pixrec));
obj.pixrec_callback(obj)
% set(h.pixrec,'String',num2str(pixrec));
% obj.updateFormatParameters;


% obj.updateParameters;
% disp('pixelsize')
notify(obj.P,'sr_render')
% obj.guiPar2Par;
end

function clickOnSrImage(handle, action,obj)
% obj.checkForSRFigure(0,0);
if isempty(obj.locData.loc) %no localizations loaded
    return
end
% if action.Button==3&&(handle==obj.guihandles.ov_axes)
%     return
% end
if action.Button==3
     resetview_callback(0,0,obj)
else
    
pos=action.IntersectionPoint*1000;
obj.setPar('sr_pos',pos);
obj.updateFormatParameters;
hfig=obj.getPar('sr_figurehandle');
figure(hfig);
notify(obj.P,'sr_render')
end
end

function clickOnSrImageW(src, action,obj)
% obj.checkForSRFigure(0,0);
if isempty(obj.locData.loc) %no localizations loaded
    return
end
% src.SelectionType
% src.CurrentObject==obj.sraxis
if ~(src.CurrentObject==obj.sraxis)
    return
end

switch src.SelectionType
    case 'normal'
        pos=obj.sraxis.CurrentPoint*1000;
        oldpos=pos(1,1:2);
        src.WindowButtonMotionFcn = @motion;
        src.WindowButtonUpFcn = @up;
        if ~isempty(obj.roicallbackid)&&~isempty(obj.roihandle)&&isvalid(obj.roihandle)
            obj.roihandle.removeNewPositionCallback(obj.roicallbackid);
            obj.roicallbackid=[];
        end
         obj.setPar('fastrender',true)
    case 'open'
        resetview_callback(0,0,obj)
    case 'alt'
        pos=obj.sraxis.CurrentPoint(1,1:2)*1000;
        updatepos(pos)
end
    function motion(src,callbackdata)
         posh=obj.sraxis.CurrentPoint;
         newpos=posh(1,1:2)*1000;
         dpos=newpos-oldpos;
         srpos=obj.getPar('sr_pos');
         if sum((dpos).^2)>1
             posax=srpos(1:2)-dpos;
             updatepos(posax)
         end
    end
    function up(src,callbackdata)
        obj.setPar('fastrender',false)
%         obj.getPar('fastrender')
%         obj.roicallbackid=obj.roihandle.addNewPositionCallback(@obj.linecallback);
        src.WindowButtonMotionFcn='';
        src.WindowButtonUpFcn='';
        if ~isempty(obj.roihandle)&&isvalid(obj.roihandle)
         obj.roicallbackid=obj.roihandle.addNewPositionCallback(@obj.linecallback);
        end
        updatepos()
    end
    function updatepos(posh)
        if nargin>0
            
        obj.setPar('sr_pos',posh);
        obj.updateFormatParameters;
        end
%         hfig=obj.getPar('sr_figurehandle');
%         figure(hfig);
        notify(obj.P,'sr_render')
    end
end


function srSizeChange(handle, action,obj)
hf=handle.Position;
rim=75;
set(handle.CurrentAxes,'Position',[rim rim hf(3)-2*rim hf(4)-2*rim])
obj.updateFormatParameters;
end

function predefkey_callback(callobj,event,obj)
global editon
    if event.Character>='0'&&event.Character<='9'
        if editon==0
            s=event.Character;
        editon=1;
        else
            s=[callobj.String event.Character];
        end
        set(callobj,'String',s)
    elseif event.Character>1
        editon=0;
    end
end

function predef_callback(callobj,event,obj)
global edition
edition=0;
s=get(callobj,'String');
set(obj.guihandles.pixrec,'String',s)
format_callback('','',obj,3)
end

function viewselect_callback(callobj,event,obj)
hfilter=obj.getPar('filterpanel');

if strcmp(obj.guihandles.overviewimage.Tag,'detached')||...
        strcmp(hfilter.Tag,'detached')
    hfilter.Visible='on';
    obj.guihandles.overviewimage.Visible='on';
    
else

    if obj.guihandles.overview_select.Value 
        hfilter.Visible='on';
        obj.guihandles.overviewimage.Visible='off';
        obj.guihandles.overview_select.String='filter->OV';
    else
        obj.guihandles.overviewimage.Visible='on';
        obj.guihandles.overview_select.String='OV->filter';
        hfilter.Visible='off';
    end
end
end


function redrawov_callback(callobj,event,obj)
% render overview image, click in overview image centers superresolution
% image on that spot
if isempty(obj.locData.loc)
    return
end
hax=obj.getPar('ov_axes');
hax.Units='pixels';
pos=hax.Position;
hax.Units='normalized';

fi=obj.getPar('currentfileinfo');
pixrec=fi.cam_pixelsize_um*1000; 
roi=fi.roi;
xext=[roi(1) roi(1)+roi(3)]* pixrec(1);
yext=[roi(2) roi(2)+roi(4)]* pixrec(end);  

% sron=obj.getPar('sr_layerson');
files=obj.locData.files.file;
for k=1:length(files)
    roi=files(k).info.roi;
xexth=[roi(1) roi(1)+roi(3)]* pixrec(1);
yexth=[roi(2) roi(2)+roi(4)]* pixrec(end); 
xext(1)=min(xext(1),xexth(1));xext(2)=max(xext(2),xexth(2));
yext(1)=min(yext(1),yexth(1));yext(2)=max(yext(2),yexth(2));
end
    
  axis(hax,'equal')



% p.sr_size=[roi(3) roi(4) ]/2* pixrec;
p.sr_size=[xext(2)-xext(1) yext(2)-yext(1)]/2;

px=(xext(2)-xext(1))/pos(3);
py=(yext(2)-yext(1))/pos(4);
p.sr_pixrec=max(px,py);

p.sr_sizeRecPix=round(p.sr_size/p.sr_pixrec*2);
p.sr_pos=[mean(xext) mean(yext)];
p.sr_axes=hax;

p.sr_plotlayernames=false;
pall=obj.getLayerParameters;
for k=1:length(pall)
    pall{k}=copyfields(pall{k},p);
end
try
TotalRender(obj.locData,pall,{'xnm','ynm'});
catch err
    err
    obj.status('Error in reconstruction. Maybe some render settings are incompatible with current data');
end
end




function paro=formatpardialog(callobj,event,obj)
persistent settings
if isempty(settings)
    settings.customcheck=false;
    settings.Pixelsize='1x1';
    settings.imsize='2000 2000';
%     settings.layerssep=false;
    settings.colorbarthickness=4;
    settings.plotscalebar=true;
    settings.customcheck=false;
%     settings.plotlayernames=false;
end
[settings, button] = settingsdlg(...
    'Description', 'SR format parameters',... 
    'title'      , 'Par',... 
    'Pixelsize',{'1x1','2x2','3x3','4x4'},...
    {'thickness of colorbar (pix)','colorbarthickness'},settings.colorbarthickness,...
    {'Plot scale bar';'plotscalebar'},[true ],...
    {'Custom image size';'customcheck'},[false true],...
    {'Imagesize (pixel) or magnification (if <20)','imsize'},num2str(settings.imsize) );
     

if strcmpi(button,'ok')
    obj.setPar('sr_imsizecheck',settings.customcheck);
    obj.setPar('sr_pixfactor',str2num(settings.Pixelsize(1)));
    if ischar(settings.imsize)
        settings.imsize=str2num(settings.imsize);
    end
    if settings.customcheck
        ims=settings.imsize;
        if ims(1)<20 %magnification
            if length(ims)==1
                ims(2)=ims(1);
            end
            imold=obj.getPar('sr_sizeRecPix');
            imnew=imold.*ims;
            mag=ims(1);
        else
            imnew=ims;
            imold=obj.getPar('sr_sizeRecPix');
            mag=ims(1)/imold(1);
        end
        obj.setPar('sr_imagesize',imnew);
        settings.imsize=imnew;
        pixelsize=obj.getPar('sr_pixrec');
        obj.setPar('sr_pixrec',pixelsize/mag);
    end
%     obj.setPar('sr_layersseparate',settings.layerssep);
    obj.setPar('sr_colorbarthickness',settings.colorbarthickness);
    obj.setPar('sr_plotscalebar',settings.plotscalebar);
%     obj.setPar('sr_plotlayernames',settings.plotlayernames);
%     if settings.newfig
%         obj.setPar('sr_axes',obj.makesrfigure((settings.fignumber)));
%     end
    obj.updateFormatParameters;
else
    paro=[];
end
end

function roi_callback(callobj,data,obj,roimode,roiposition)
% Define ROI in reconstructed superresolution image (point, linear,
% elliptical, rectengular, polynomial or free hand). Used in numerous
% plugins.
global roimodecallback
roimodecallback=roimode;
p=obj.getGuiParameters;
sr_axes=obj.sraxis;
f=sr_axes.Parent;
% oldbdf=f.WindowButtonDownFcn;

if nargin<5
    roiposition=[];
    f.WindowButtonDownFcn='';
elseif ~p.roishow %if restore and not show: dont do anything
    return
end
delete(obj.roihandle)


xlim=sr_axes.XLim;
ylim=sr_axes.YLim;
switch roimode
    case {1,'imrect'}
        h=imrect(sr_axes,roiposition);
    case {2,'imellipse'}
        h=imellipse(sr_axes,roiposition);
    case {3,'imfreehand'}
        h=imfreehand(sr_axes,roiposition);
    case {4,'imline'} %line
        h=imline(sr_axes,roiposition);        
    case {5,'impoint'}
        h=impoint(sr_axes,roiposition);
    case {6,'impoly'}
        h=impoly(sr_axes,roiposition);
    case 0
        if p.roishow
            roi_callback(callobj,data,obj,class(obj.roihandle),obj.roiposition);
            return
        else
            h=[];
        end
    case -1
        h=[];
        obj.roihandle=h;
        obj.setPar('sr_roihandle',obj.roihandle);
        mg=obj.getPar('mainGui');
        mg.children.guiRender.temproi=h;
end
f.WindowButtonDownFcn={@clickOnSrImageW,obj};
if ~isempty(h)
    obj.roicallbackid=addNewPositionCallback(h,@obj.linecallback);
    obj.roihandle=h;
    obj.setPar('sr_roihandle',obj.roihandle);
    obj.linecallback(obj.roihandle.getPosition)
end
sr_axes.XLim=xlim;
sr_axes.YLim=ylim;

end

function resetview_callback(oject,data,obj)
  si=obj.getPar('sr_sizeRecPix');
  if ~isempty(obj.locData.loc)&&~isempty(obj.locData.loc.xnm)
    mx=myquantilefast(obj.locData.loc.xnm,[0.9995,0.0005],100000);
    maxx=mx(1);minx=mx(2);
    my=myquantilefast(obj.locData.loc.ynm,[0.9995,0.0005],100000);
    maxy=my(1);miny=my(2);
  else
      disp('cannot find size of image, no reset')
      return
  end
  obj.setPar('sr_pos',[(maxx+minx)/2 (maxy+miny)/2]);
  pixrec=round(max((maxx-minx)/si(1),(maxy-miny)/si(2)));
  obj.setPar('sr_pixrec',pixrec);
  obj.pixrec_callback(obj)
end

function lw_callback(oject,data,obj)
lw=obj.getSingleGuiParameter('linewidth_roi');
obj.setPar('linewidth_roi',lw);
if ~isempty(obj.roihandle)
obj.linecallback(obj.roihandle.getPosition);
end
end

function detach_callback(a,b,obj,handle)
if strcmp(handle.Tag,'detached')
    maingui=obj.getPar('mainGuihandle');
    fold=handle.Parent;
    handle.Parent=maingui;
    handle.Units='pixel';
    handle.Position=[0 400 430 350];
    handle.Tag='attached';
    obj.guihandles.ov_axes.Units='normalized';
    obj.guihandles.ov_axes.Position=[0 0 1 1];
    close(fold);
else
f=figure('MenuBar','none','Toolbar','none');
handle.Parent=f;
handle.Position(1)=0;
handle.Position(2)=0;
f.Position(3:4)=handle.Position(3:4);
handle.Tag='detached';
handle.Units='normalized';
obj.guihandles.ov_axes.Units='normalized';
end
end
