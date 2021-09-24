classdef DisplayImageTagsTiff<interfaces.DialogProcessor
%     Shows raw camera images saved with the localization data.
    properties
        imageloader
    end
    methods
        function obj=DisplayImageTagsTiff(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
            obj.showresults=true;
        end
        
        function out=run(obj,p)
            reader=obj.imageloader.reader;
            img=reader.getImage(0,0,1,0);
            imgposn=3;
            if isempty(img)
                img=reader.getImage(0,1,0,0);
                imgposn=2;
            end
            if isempty(img)
                img=reader.getImage(0,0,0,1);
                imgposn=4;
            end
            
            readtags=split(obj.getSingleGuiParameter('fields'),',');

            frame=1;
            imgpos={0,0,0,0}; imgpos{imgposn}=frame-1;
            imgmeta=reader.getImageTags(imgpos{:});
            disp('read data, this might take a while')
            while ~isempty(imgmeta)
                for k=1:length(readtags)
                    val=imgmeta.get(readtags{k});
                    if ischar(val)
                        val=str2double(val);
                    end
                    data(k,frame)=val;
                end
                frame=frame+1;
                imgpos{imgposn}=frame-1;
                imgmeta=reader.getImageTags(imgpos{:});
            end
            disp('meta data reading done')
            for k=1:length(readtags)
                    name=readtags{k}(1:min(length(readtags{k}),7));
                    ax=obj.initaxis([name ':' num2str(k)]);
%                     tab=uitab(tg,'Title',imagetags.tags{k});
%                     ax=axes('Parent',tab);
                    dat=data(k,:);
%                     if ischar(dat{1})
%                         datm=str2double(dat);
%                     else
%                         datm=cell2mat(dat);
%                     end
                    frames=(1:length(dat))';
                    frameind=dat~=0;
                    plot(ax,frames(frameind), dat(frameind))
                    xlabel('frame')
                    ylabel(readtags{k})
                    xlim([min(frames(frameind)) max(frames(frameind))])
            end
           out.data=data;
           out.tags=readtags;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function loadtiff_callback(a,b,obj)
lastf=obj.getPar('lastSMLFile');
if ~isempty(lastf)
    path=fileparts(lastf);
else
    path='';
end
[f path]=uigetfile([path filesep '*.tif']);
if f
    obj.setGuiParameters(struct('filename',[path f]));
     obj.imageloader=imageloaderAll(obj.getSingleGuiParameter('filename'));
end
end

function selectfield_callback(a,b,obj)

md=obj.imageloader.getmetadatatags;
tag=gettag(md);

if ~isempty(tag)
    fields=obj.getSingleGuiParameter('fields');
    if ~isempty(fields)
        fields=[fields ','];
    end
    fields=[fields tag{1}];
    obj.setGuiParameters(struct('fields',fields))
end
end



function tag = gettag(ma)
f=figure;
tc=uitable(f);
[~,ind]=sortrows(ma(:,1));


tc.Data=ma(ind,:);
tc.Position(3)=f.Position(3)-30;
tc.Position(2)=100;
tc.CellSelectionCallback=@cellselecth;
w=tc.Position(3);
tc.ColumnWidth={w*.75,.2*w};
uicontrol('Style','pushbutton','String','Ok','Position',[200 10 100 20],'Callback',@buttoncallback)
uicontrol('Style','pushbutton','String','Cancel','Position',[10 10 100 20],'Callback',@buttoncallback)
pos=[];
waitfor(f)

    function buttoncallback(a,b)
        if ~isempty(pos)&&strcmp(a.String,'Ok')
        tag=tc.Data(pos(1),:);
        else
            tag={};
        end
        close(f)
    end
    function cellselecth(t,d)
        pos=d.Indices;
    end
end


function pard=guidef(obj)
pard.filename.object=struct('Style','edit','String','');
pard.filename.position=[1,1];
pard.filename.Width=3;

pard.load_tiff.object=struct('Style','pushbutton','String','Load','Callback',{{@loadtiff_callback,obj}});
pard.load_tiff.position=[1,4];
pard.load_tiff.Width=1;

pard.fields.object=struct('Style','edit','String','');
pard.fields.position=[2,1];
pard.fields.Width=3;

pard.fieldselect.object=struct('Style','pushbutton','String','select Field','Callback',{{@selectfield_callback,obj}});
pard.fieldselect.position=[2,4];
pard.fieldselect.Width=1;


pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Show image tags for MM tiff file.';
end

