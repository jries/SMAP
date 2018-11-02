classdef fileFeeder<interfaces.WorkflowModule;
    properties
        
    end
    methods
       function obj=fileFeeder(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=0; 
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            
        end
        function output=run(obj,data,p)
            %afterwards run list
            filelist=obj.guihandles.filelist.String;
            if iscell(filelist)
                for k=1:length(filelist)
                    data=interfaces.WorkflowData;
                    data.data=filelist{k};
                    
                    obj.output(data);
                end
                data.eof=true;
                data.data=[];
                output=data;
            end
        end
    
    end
end


function add_callback(a,b,obj)
filelist=obj.guihandles.filelist.String;
if isempty(filelist)
    fs=obj.getPar('lastSMLFile');
    filelist={};
else
    if ~iscell(filelist)
        filelist={filelist};
    end
    fs=filelist{1};
end
if isempty(fs)
    fs='*.*';
end
[f,p]=uigetfile(fs,'MultiSelect','on');
if f
    if ~iscell(f)
        f={f};
    end
    for k=1:length(f);
        filelist{end+1}=[p f{k}];
    end
    obj.guihandles.filelist.String=filelist;
    v=min(max(1,obj.guihandles.filelist.Value),length(filelist));
   obj.guihandles.filelist.Value=v;
end

end
function remove_callback(a,b,obj)
v=obj.guihandles.filelist.Value;
filelist=obj.guihandles.filelist.String;
filelist(v)=[];
v=min(max(1,obj.guihandles.filelist.Value),length(filelist));
obj.guihandles.filelist.Value=v;
obj.guihandles.filelist.String=filelist;

end
function clear_callback(a,b,obj)
obj.guihandles.filelist.String={};
end
function pard=guidef(obj)
pard.filelist.object=struct('Style','listbox','String','');
pard.filelist.position=[5,1];
pard.filelist.Height=5;
pard.filelist.Width=4;

pard.add.object=struct('Style','pushbutton','String','add','Callback',{{@add_callback,obj}});
pard.add.position=[6,1];
pard.add.Width=0.5;

pard.remove.object=struct('Style','pushbutton','String','remove','Callback',{{@remove_callback,obj}});
pard.remove.position=[6,2];
pard.remove.Width=0.5;

pard.clear.object=struct('Style','pushbutton','String','clear','Callback',{{@clear_callback,obj}});
pard.clear.position=[6,4];
pard.clear.Width=0.5;

end