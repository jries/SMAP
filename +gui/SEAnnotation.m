classdef SEAnnotation< interfaces.SEProcessor
%     Annotate ROIs using lists of user-defined tags or up to four
%     user-defined linear or arbitrary image-ROIs.
    properties
%         SEpreview
%         list
    end
    methods
        function obj=SEAnnotation(varargin)   
            obj@interfaces.SEProcessor(varargin{:})
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.SEProcessor(obj);
            set(obj.guihandles.loadlist,'Callback',{@loadlist_callback,obj})
            obj.guihandles.list1.KeyPressFcn={@keypress,obj,1};
            obj.guihandles.list2.KeyPressFcn={@keypress,obj,2};
            obj.guihandles.list3.KeyPressFcn={@keypress,obj,3};
            obj.guihandles.list4.KeyPressFcn={@keypress,obj,4};
            
            obj.guihandles.list1.Callback={@list_callback,obj,1};
            obj.guihandles.list2.Callback={@list_callback,obj,2};
            obj.guihandles.list3.Callback={@list_callback,obj,3};
            obj.guihandles.list4.Callback={@list_callback,obj,4};
            obj.guihandles.comments.Callback={@comments_callback,obj};
            
            obj.guihandles.line1.Callback={@line_callback,obj,1};
            obj.guihandles.line2.Callback={@line_callback,obj,2};
            obj.setPar('ROI_lineannotation_handle_1',obj.guihandles.line1);
            obj.setPar('ROI_lineannotation_handle_2',obj.guihandles.line2);
            obj.makeinfobutton('ne');
%             set(obj.guihandles.redrawall,'Callback',{@redrawall_callback,obj})
%             set(obj.guihandles.clearall,'Callback',{@clearall_callback,obj})
%             addlistener(obj.SE.locData,'loaded',@obj.loaded_notify);
        end
        function attachLocData(obj,locData)
            attachLocData@interfaces.SEProcessor(obj,locData);
            maing=obj.getPar('mainGuihandle');
            maing.KeyPressFcn={@keypress,obj,0};
%             addlistener(obj.SE.locData,'loaded',@obj.loaded_notify);
        end
        
        function loaded_notify(obj,lb,eventdata)
            a=obj.SE.sites(1).annotation;
            if ~isempty(a)
            entries{1}=a.list1.string;
            entries{2}=a.list2.string;
            entries{3}=a.list3.string;
            entries{4}=a.list4.string;
            setlist(obj,entries)
            end
        end
        
        function sitechange(obj,site)        
            obj.guihandles.list1.Value=site.annotation.list1.value;
            obj.guihandles.list2.Value=site.annotation.list2.value;
            obj.guihandles.list3.Value=site.annotation.list3.value;
            obj.guihandles.list4.Value=site.annotation.list4.value;
            setlistentries(obj,site);
            obj.guihandles.comments.String=site.annotation.comments;
            
            
            alphaimage=site.image.angle;
            pos=site.annotation.line1.pos;
            angle=pos2angle(pos);%+alphaimage;
            len=sqrt(sum((pos(1,:)-pos(2,:)).^2))*1000;
             obj.guihandles.line1.String=[num2str(angle,'%3.1f') ', ' num2str(len,'%3.0f')];
             
             pos=site.annotation.line2.pos;
            angle=pos2angle(pos);%+alphaimage;
            len=sqrt(sum((pos(1,:)-pos(2,:)).^2))*1000;
             obj.guihandles.line2.String=[num2str(angle,'%3.1f') ', ' num2str(len,'%3.0f')];
             if isfield(site.annotation,'use')
                obj.guihandles.usesite.Value=site.annotation.use;
             end
            
        end
        function usesite_callback(obj,object,b)
            selected=obj.SE.processors.preview.guihandles.sitelist.Value;
            if length(selected)<=1
            
            site=obj.SE.currentsite;
            site.annotation.use=object.Value;
            else
                for k=1:length(selected)
                   site=obj.SE.sites(selected(k));
                    site.annotation.use=object.Value;
                end
            end
            obj.SE.processors.preview.updateSitelist;
        end
        function drawroi_callback(obj,a,b)
            p=obj.getSingleGuiParameter('roiselect');
             obj.SE.processors.preview.lineannotation(3,p.selection);
        end
%         function updateSingleParameter(obj, data,actionData,field)
%             val=obj.getSingleGuiParameter(field);
% %             obj.SE.sePar.(data.Parent.Title).(field)=val;
%         end
    end
end

function list_callback(a,b,obj,list)
    prev=obj.SE.processors.preview;
    siteind=prev.guihandles.sitelist.Value;
    if length(siteind)>1
        for k=1:length(siteind)
            site=obj.SE.sites(siteind(k));
            listname=['list' num2str(list)];
            site.annotation.(listname).string=obj.guihandles.(listname).String;
            site.annotation.(listname).value=obj.guihandles.(listname).Value;
        end
    else
         site=obj.SE.currentsite;
        listname=['list' num2str(list)];
        site.annotation.(listname).string=obj.guihandles.(listname).String;
        site.annotation.(listname).value=obj.guihandles.(listname).Value;  
    end
    prev.updateSitelist;
end

function comments_callback(a,b,obj)
    site=obj.SE.currentsite;
    site.annotation.comments=obj.guihandles.comments.String;
end

function keypress(object,data,obj,focus)
if strcmpi(data.Key ,'rightarrow')
    if strcmpi(data.Modifier,'shift')
        newlist=focus+1; 
        if newlist<1
            newlist=4;
        end
        if newlist>4
            newlist=1;
        end
        uicontrol(obj.guihandles.(['list' num2str(newlist)]))
        
    else
        obj.SE.processors.preview.nextsite(1);
    end
end
if strcmpi(data.Key ,'leftarrow')
    if strcmpi(data.Modifier,'shift')
        newlist=focus-1; 
        if newlist<1
            newlist=4;
        end
        if newlist>4
            newlist=1;
        end
        uicontrol(obj.guihandles.(['list' num2str(newlist)]))
        
    else
        obj.SE.processors.preview.nextsite(-1);
    end
end
end

function loadlist_callback(a,b,obj)
[f p]=uigetfile('settings/*.txt');
if f
    loadlist(obj,[p f])
end
end

function loadlist(obj,f)
% [~,fn]=fileparts(f);
% entries=eval(fn);
fid=fopen(f);
ind=1;
tline = fgets(fid);
while ischar(tline)
    entries(ind)=textscan(tline,'%s','Delimiter',',');
    tline = fgets(fid);
    ind=ind+1;
end

setlist(obj,entries);
end

function setlist(obj,entries)
obj.guihandles.list1.String=entries{1};
obj.guihandles.list2.String=entries{2};
obj.guihandles.list3.String=entries{3};
obj.guihandles.list4.String=entries{4};
obj.guihandles.list1.Value=1;
obj.guihandles.list2.Value=1;
obj.guihandles.list3.Value=1;
obj.guihandles.list4.Value=1;

for k=1:length(obj.SE.sites)
    setlistentries(obj,obj.SE.sites(k))
end


end

function setlistentries(obj,site)
    for l=1:4
    	site.annotation.(['list' num2str(l)]).string=obj.guihandles.(['list' num2str(l)]).String;
    end
end

function line_callback(a,b,obj,linenumber)
 obj.SE.processors.preview.lineannotation(linenumber,obj.guihandles.(['line' num2str(linenumber)]));
end

function reverseline(a,b,obj,linenumber)
    site=obj.SE.currentsite;
    pos=site.annotation.(['line' num2str(linenumber)]).pos;
    site.annotation.(['line' num2str(linenumber)]).pos=pos([2 1],:);
    site.annotation.(['line' num2str(linenumber)]).angle=mod(site.annotation.(['line' num2str(linenumber)]).angle+180+180,360)-180;
    obj.SE.processors.preview.lineannotation(linenumber,obj.guihandles.(['line' num2str(linenumber)]));
end

function pard=guidef(obj)
pard.list1.object=struct('Style','listbox','String','1|2|3|4|5|6|7|8');
pard.list1.position=[6,1.];
pard.list1.Height=5;

pard.list2.object=struct('Style','listbox','String','l1|l2');
pard.list2.position=[12,1.];
pard.list2.Height=5;

pard.list3.object=struct('Style','listbox','String',{'empty'});
pard.list3.position=[6,2.];
pard.list3.Height=5;

pard.list4.object=struct('Style','listbox','String',{'empty'});
pard.list4.position=[12,2.];
pard.list4.Height=5;

pard.list1_title.object=struct('Style','text','String','list1');
pard.list1_title.position=[1,1.];

pard.list2_title.object=struct('Style','text','String','list2');
pard.list2_title.position=[7,1.];

pard.list3_title.object=struct('Style','text','String','list3');
pard.list3_title.position=[1,2];

pard.list4_title.object=struct('Style','text','String','list4');
pard.list4_title.position=[7,2];


pard.loadlist.object=struct('Style','pushbutton','String','<-load');
pard.loadlist.position=[2,3];
pard.loadlist.Width=0.5;

pard.line1.object=struct('Style','pushbutton','String','line 1');
pard.line1.position=[4,3];
pard.line1.Height=1.5;
pard.line1.Width=0.8;

pard.line1switch.object=struct('Style','pushbutton','String','<>','Callback',{{@reverseline,obj,1}});
pard.line1switch.position=[4,3.8];
pard.line1switch.Height=1.5;
pard.line1switch.Width=0.2;


pard.line2.object=struct('Style','pushbutton','String','line 2');
pard.line2.position=[4,4];
pard.line2.Height=1.5;
pard.line2.Width=0.8;

pard.line2switch.object=struct('Style','pushbutton','String','<>','Callback',{{@reverseline,obj,2}});
pard.line2switch.position=[4,4.8];
pard.line2switch.Height=1.5;
pard.line2switch.Width=0.2;

pard.roiselect.object=struct('Style','popupmenu','String',{{'rectangle','ellipse','polygon','polyline','free'}});
pard.roiselect.position=[6,3];
pard.roiselect.Height=1.5;
pard.roi.object=struct('Style','pushbutton','String','ROI','Callback',@obj.drawroi_callback);
pard.roi.position=[6,4];
pard.roi.Height=1.5;

pard.usesite.object=struct('Style','checkbox','String','use site','Callback',@obj.usesite_callback);
pard.usesite.position=[2,4];
pard.usesite.Height=1;


pard.t1.object=struct('Style','text','String','Comments:');
pard.t1.position=[9,3];

pard.comments.object=struct('Style','edit','String','','Max',5);
pard.comments.position=[12,3];
pard.comments.Height=3;
pard.comments.Width=2;
pard.helpfile='SMAP.Gui.ROIannotation.txt';
pard.plugininfo.name='ROIannotation';
end