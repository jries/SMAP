classdef GuiHistSlider< interfaces.LayerInterface
%     Show histogram of values of localization field and use sliders to
%     adjust filter range
    properties
        field
        data
        histogram={};
        backup
        

    end
    methods
        function obj=GuiHistSlider(varargin)
            obj@interfaces.LayerInterface(varargin{:});
        end
        function makeGui(obj)
            fontsize=obj.guiPar.fontsize;
            height=0.1;
            h.histax=axes('Parent',obj.handle,'Position',[0.1,0.34,0.68,0.58],'FontSize',fontsize*0.75);
            h.smin=uicontrol('Parent',obj.handle,'Units','normalized','Style','slider','Position',[0.05,0.13,0.8,height*0.7]);
            h.smax=uicontrol('Parent',obj.handle,'Units','normalized','Style','slider','Position',[0.05,0.03,0.8,height*0.7]);
            h.field=uicontrol('Parent',obj.handle,'Units','normalized','Style','text','Position',[0.8,0.85,0.2,height],'FontSize',fontsize);
            h.filteron=uicontrol('Parent',obj.handle,'Units','normalized','Style','checkbox','Position',[0.8,0.7,0.2,height],'String','filter','FontSize',fontsize*0.7,'Callback',{@showchanged_callback,obj});
            h.filteron.TooltipString='if checked, this field is filtered';
            
            h.lockrange=uicontrol('Parent',obj.handle,'Units','normalized','Style','checkbox','Position',[0.8,0.5,0.2,height],'String','range fix','FontSize',fontsize*0.7);
            h.lockrange.TooltipString='if checked, the range will be restricted to the width given here';
            
            h.range=uicontrol('Parent',obj.handle,'Units','normalized','Style','edit','Position',[0.8,0.4,0.2,height],'String','0','FontSize',fontsize);
            h.range.TooltipString='width of range of values';
            
            h.vmin=uicontrol('Parent',obj.handle,'Units','normalized','Style','edit','Position',[0.85,0.16,0.15,height],'String','0','FontSize',fontsize);
            h.vmax=uicontrol('Parent',obj.handle,'Units','normalized','Style','edit','Position',[0.85,0.03,0.15,height],'String','0','FontSize',fontsize);
            set(h.vmin, 'Callback',{@vchanged_callback,obj,'min',1});
            set(h.vmax, 'Callback',{@vchanged_callback,obj,'max',1});
            
            set(h.smin, 'Callback',{@slider_callback,obj,'min'});
            set(h.smax, 'Callback',{@slider_callback,obj,'max'});
            
            h.autoupdate=uicontrol('Parent',obj.handle,'Units','normalized','Style','checkbox','Position',[0.8,0.6,0.2,height],'String','auto update','FontSize',fontsize*0.7);
            
            h.invert=uicontrol('Parent',obj.handle,'Units','normalized','Style','checkbox','Position',[0.8,0.3,0.2,height],'String','invert','FontSize',fontsize*.7,'Callback',{@showchanged_callback,obj});
            h.invert.TooltipString='Invert selection. Good for trying out filters.';
            
            
            h.handle=obj.handle;
            obj.guihandles=h;
            callobj=obj;
            obj.addSynchronization([obj.layerprefix 'selectedField'],[],[],{@callobj.selectedField_callback})   
            callobj=obj;
            obj.addSynchronization('locFields',[],[],{@callobj.updateGui}) 
        end
        function updateGui(obj)
%             obj.field=[];
        end
        function selectedField_callback(obj,direct)
            if nargin<2
                direct=false;
            end
            sfield=obj.getPar('selectedField','layer',obj.layer);
            if length(sfield)<5||~sfield{5}
                return
            end
            filenumber=obj.getPar(['layer' num2str(obj.layer) '_']).ch_filelist.Value;
            fieldh=sfield{1};
            if isempty(obj.locData.loc)
                return
            end
%             if ~isfield(obj.locData.loc,fieldh)
                switch fieldh
                    case 'shiftxy'
                        obj.guihandles.smin.Min=-500;
                        obj.guihandles.smin.Max=500;
                        obj.guihandles.smax.Min=-500;
                        obj.guihandles.smax.Max=500;
                    case 'imax'
                        lp=obj.getPar('','layer',obj.layer);
                        imaxtoggle=lp.imaxtoggle;
                        if imaxtoggle
                        obj.guihandles.smin.Min=-5;
                        obj.guihandles.smin.Max=-1;
                        obj.guihandles.smax.Min=-1;
                        obj.guihandles.smax.Max=-1+0.001;
                        obj.guihandles.smax.Value=-1;
                        obj.guihandles.vmax.String=num2str(-1+.001);
                        obj.guihandles.field.String='Contrast (q)';
                        else
                            img=obj.locData.layer(obj.layer).images.srimage.image;
                            maximg=max(img(:))*1.4;
                        obj.guihandles.smin.Min=0;
                        obj.guihandles.smin.Max=maximg;
                        obj.guihandles.smax.Min=maximg;
                        obj.guihandles.smax.Max=maximg*1.001;
                        obj.guihandles.smax.Value=maximg*1.001;
                        obj.guihandles.vmax.String=num2str(maximg*1.001);
                        obj.guihandles.field.String='Contrast (Imax)';
                        end
                    otherwise
                        obj.guihandles.field.String=fieldh;
%                         if sfield{2}==sfield{3}
%                             sfield{2}=sfield{2}-1;
%                         end
                        obj.guihandles.vmin.String=num2str(sfield{2});
                        obj.guihandles.vmax.String=num2str(sfield{3});
                end

            if length(sfield)>3&&~isempty(sfield{4})
                obj.guihandles.filteron.Value=sfield{4}(1);
                if length(sfield{4})>1
                    obj.guihandles.invert.Value=sfield{4}(2);
                end
            end
            
            if isfield(obj.backup,fieldh)&& myisfield(obj.locData,fieldh)&&~strcmpi(fieldh,'colorfield')...
                    &&obj.backup.(fieldh).len==length(obj.locData.loc.(fieldh))...
                    &&obj.backup.(fieldh).histogram.filenumber == filenumber; 
                restorestruc=obj.backup.(fieldh);
                restore=true;
            else
                restore=false;
            end
             %backup
            if ~isempty(obj.field) && myisfield(obj.locData.loc,fieldh) %first tiem might be empty
                obj.backup.(obj.field).lockrange=obj.guihandles.lockrange.Value;
%                 obj.backup.(obj.field).quantile=obj.data.quantile;
                obj.backup.(obj.field).range=obj.guihandles.range.String;
                obj.backup.(obj.field).histogram=obj.histogram;
                obj.backup.(obj.field).len=length(obj.locData.loc.(fieldh));
            end
%             if obj.field==fieldh
%                 restore=false;
%             end
            obj.field=fieldh;
            
%             if strcmp(fieldh,'shiftxy')
                
            if restore  && ~isempty(restorestruc.histogram)&&~direct
%                 obj.guihandles.lockrange.Value=restorestruc.lockrange;
                obj.guihandles.range.String=restorestruc.range;
%                 obj.data.quantile=restorestruc.quantile;
                q=restorestruc.histogram.x([1 end]);
                obj.guihandles.smin.Min=q(1);
                obj.guihandles.smin.Max=q(2);
                obj.guihandles.smax.Min=q(1);
                obj.guihandles.smax.Max=q(2);
                obj.histogram=restorestruc.histogram;        
            elseif isfield(obj.locData.loc,fieldh) && isfield(obj.locData.loc,'filenumber')
                    v=double(obj.locData.loc.(fieldh)(obj.locData.loc.filenumber==filenumber));
                    v=real(v);
                    v(isinf(v))=[];
                    q=getquantile(v);
                    if q(2)==q(1)
                        q(1)=q(1)-1;
                        q(2)=q(2)+1;
                    end
%                    if imag(q)~=0
%                        q=real(q);
%                    end
%                     obj.data.quantile=q;
%                     obj.guihandles.lockrange.Value=0;
%                     obj.guihandles.autoupdate.Value=0;
                    set(obj.guihandles.smin,'Min',q(1),'Max',q(2),'Value',q(1))
                    set(obj.guihandles.smax,'Min',q(1),'Max',q(2),'Value',q(2))
                    [obj.histogram.x,obj.histogram.hist]=makeHist(v,q);
                    obj.histogram.filenumber=filenumber;

            else
%                 obj.guihandles.lockrange.Value=0;
                dx=(-obj.guihandles.smin.Min+obj.guihandles.smin.Max)/50;
                obj.histogram.x=obj.guihandles.smin.Min:dx:obj.guihandles.smin.Max;
                obj.histogram.hist=ones(size(obj.histogram.x));
%                 obj.histogram.x
            end
            
            if ~isempty(obj.locData.loc)&&isfield(obj.locData.loc,fieldh)
                obj.data.values=double(obj.locData.loc.(fieldh));
            else
                plot(0,0,'Parent',obj.guihandles.histax)
%                 obj.data.quantile=[0 1];
                obj.data.values=0;
            end
            vchanged_callback([],0, obj,'minmax',false)
            
             
        end
    end
end
    

function q=getquantile(v)
mv=max(v);
if min(v)==mv
    q=[0.9 1.1 ]*double(mv);
else
    try
%          q=myquantilefast(v,[0.015,0.985],10000);
         q=myquantilefast(v,[0.05,0.95],10000);
         q=q+abs(q(2)-q(1))*[-1 1]*.05;
    catch
    end
    if isempty(q)||any(isnan(q))
        q=double([min(v) mv]);
        if isempty(q)
            q=[0 1];
        end
    end
end
end

function vchanged_callback(object,data, obj,minmax,direct)
if isempty(obj.locData.loc)
    return
end

p=obj.getGuiParameters;
filenumber=obj.getPar(['layer' num2str(obj.layer) '_']).ch_filelist.Value;
if isnan(p.vmin)
    ph.vmin=-inf;
else
    ph.vmin=p.vmin;
end
if isnan(p.vmax)
    ph.vmax=inf;
else
    ph.vmax=p.vmax;
end
if p.lockrange
    switch minmax
        case {'min'}
            ph.vmax=ph.vmin+p.range;
        case 'max'
            ph.vmin=ph.vmax-p.range;
        case 'minmax'
            ph.vmin=ph.vmax-p.range;         
    end            
else
   ph.range=ph.vmax-ph.vmin; 
end

if ph.vmin==ph.vmax
    ph.vmax=ph.vmin+1e-4;
end
hsmin=obj.guihandles.smin;
ph.smin=min(max(ph.vmin,hsmin.Min),hsmin.Max);
hsmax=obj.guihandles.smax;
ph.smax=min(max(ph.vmax,hsmax.Min),hsmax.Max);

obj.setGuiParameters(ph);

% q=obj.data.quantile;
if isempty(obj.histogram)|| (myisfield(obj.histogram,'filenumber')&&obj.histogram.filenumber ~=filenumber)%check if filenumber has changed   
   if length(obj.data.values)==length(obj.locData.loc.filenumber)
    v=double(obj.data.values(obj.locData.loc.filenumber==filenumber));
   else v=0:0.1:1;
   end
    q=getquantile(v);
    [obj.histogram.x,obj.histogram.hist]=makeHist(v,q);
    obj.histogram.filenumber=filenumber;
else
    q=obj.histogram.x([1 end]);
end
x=obj.histogram.x;his=obj.histogram.hist;

ind1=find(x>p.vmin,1,'first');
ind2=find(x<=p.vmax,1,'last');
h2=his(ind1:ind2);
x2=x(ind1:ind2);
set(obj.guihandles.histax,'NextPlot','replace')

stairs(x,his,'k','Parent',obj.guihandles.histax,'LineWidth',1)

set(obj.guihandles.histax,'NextPlot','add')
bar(x2,h2,1,'Parent',obj.guihandles.histax,'g')
if q(1)~=q(2)&&~isnan(q(1))
    obj.guihandles.histax.XLim=[q(1) q(2)];
end

if strcmp(obj.field,'colorfield')
    vlayer=obj.locData.loc.colorfield;
else
    vlayer=obj.locData.getloc(obj.field,'layer',obj.layer).(obj.field);
end
hlayer=double(hist(vlayer,x));
hlayer(1,1)=0;hlayer(1,end)=0;
plot(x(:),hlayer(:)/max(hlayer)*max(his),'r','Parent',obj.guihandles.histax)

if direct 
    sf={obj.field,ph.vmin,ph.vmax,p.filteron,true};
    obj.setPar('selectedField',sf,'layer',obj.layer)
end

if direct
    if obj.guihandles.autoupdate.Value 
        notify(obj.P,'sr_render')
    end
end
end


function [x,h]=makeHist(v,q)
d=(q(2)-q(1))/50;
if d<2&&d>0.2
    d=1;
end
if d==0
    x=q(1)-1:q(2)+1;
else
x=q(1):d:q(2);
end
h=hist(v,x);
end

function slider_callback(object,data,obj,minmax)
h=obj.guihandles;
if h.smin.Value==h.smin.Min
    if h.smin.Value>0
        h.smin.Value=0;
    else
    h.smin.Value=-inf;
    end
end
if h.smax.Value==h.smax.Max
    h.smax.Value=inf;
end
switch minmax
    case 'min'
        h.vmin.String=num2str(h.smin.Value);
    case 'max'
        h.vmax.String=num2str(h.smax.Value);
end
vchanged_callback(object,data, obj,minmax,1)
 obj.selectedField_callback(1);
end

function showchanged_callback(object,data,obj)
p=obj.getGuiParameters;
sf={obj.field,p.vmin,p.vmax,[p.filteron,p.invert],true};
obj.setPar('selectedField',sf,'layer',obj.layer)
    if obj.guihandles.autoupdate.Value 
        notify(obj.P,'sr_render')
    end
end
% 
% function checkbox_callback(object,data,name,obj)
% %     obj.(name).(obj.field)=object.Value;
% end