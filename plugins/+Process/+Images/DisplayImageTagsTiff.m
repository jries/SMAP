classdef DisplayImageTagsTiff<interfaces.DialogProcessor
%     Shows raw camera images saved with the localization data.
    properties
        imageloader
    end
    methods
        function obj=DisplayImageTagsTiff(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
            obj.showresults=true;
            initMM(obj);
        end
        
        function out=run(obj,p)
            readtags=split(obj.getSingleGuiParameter('fields'),',');
            disp('read data, this might take a while')
            md=getMetadataTifMM(obj.getSingleGuiParameter('filename'),readtags,true(1,length(readtags)));
            disp('meta data reading done')
            for k=1:length(readtags)
                    name=readtags{k}(1:min(length(readtags{k}),7));
                    ax=obj.initaxis([name ':' num2str(k)]);
                    dat=md{k};
                    frames=(1:length(dat))';
                    frameind=dat~=0;
                    plot(ax,frames(frameind), dat(frameind))
                    xlabel('frame')
                    ylabel(readtags{k})
                    xlim([min(frames(frameind)) max(frames(frameind))])
            end
           out.data=md;
           out.tags=readtags;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end


function loadtiff_callback(a,b,obj)
filename=obj.getSingleGuiParameter('filename');
if isempty(filename)
    lastf=obj.getPar('lastSMLFile');
    if ~isempty(lastf)
        path=fileparts(lastf);
    else
        path='';
    end
else
    path=fileparts(filename);
end
[f path]=uigetfile([path filesep '*.tif']);
if f
    obj.setGuiParameters(struct('filename',[path f]));
end
end

function selectfield_callback(a,b,obj)
md=getMetadataTifMM(obj.getSingleGuiParameter('filename'));
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

