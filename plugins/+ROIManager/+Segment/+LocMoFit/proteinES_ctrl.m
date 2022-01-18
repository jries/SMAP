classdef proteinES_ctrl<interfaces.DialogProcessor&interfaces.SEProcessor
    % PLUGIN_TEMPLATE Summary of this plugin goes here
    % put a description of your plugin here.
        %replace Plugin_Template by filename   
    properties
        %define class properties if needed
    end
    methods
        function obj=proteinES_ctrl(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=false; %set true, if results are shown by default
            obj.history=false; %if set true, every time the plugin is called, its parameters are added to the history. Makes sense only if plugin changes the data
        end       
        function initGui(obj)
            %is called after the GUI, defined in guidef, is made. Here you
            %can add additional GUI components, e.g. some that are not
            %defined by simple uicontrols (e.g. uitable), set additional
            %synchronizations and callbacks.
        end
        function out=run(obj,p)
            bu = basicUnit(p.capDepth, p.capShape, p.bottomDepth, p.bottomShape, p.radius);
            pes = proteinES(p.proteinID, p.proteinName, bu, p.capShape, p.bottomShape, p.velocityEPM, p.zmt0);
            timePar = readtable(p.timeParPath);
            pes.addTimePoint(timePar);
            %Here you implement the functionality
            %p contains:
            %global parameters defined in pard.inputParameters
            %uicontrol parameters, already parsed (e.g. edit with numerical
            %value returns a double, not a text. popupmenu etc return three
            %fields: String (all entries), Value (selected item), selection
            %(string of selected item).
                       
            %set a global parameter
            if ~isfield(obj.P.par, 'proteinES')
                obj.setPar('proteinES', [])
            end
             obj.P.par.proteinES.content.(pes.proteinName)=pes;
            %this example updates all lists in other plugins which are
            %linked to the fields of localization data.
            
            out=[]; %no output
            out.clipboard={'results1',3,'text1'}; % out.clipboard is copied to clipboard, separated by tabs.
            out.error='this error occured because you did something wrong'; %if an error occoured, you can output it in the status bar with this command.

        end
        function pard=guidef(obj)
            %pard structure can be =[]; All fields are optional.
            
            %define your GUI: for every GUI control define a structure as
            %follows, replace 'guiobject' by a name by which you want to
            %access the parameters
            % required fields: .object=struct(...). Defines a uicontrol
            % object. Just pass on all arguments you would otherwise pass
            % on to uicontrol. Careful: if you want to pass on a cell
            % array, you need to put it into double brackets: {{'a','b'}}.
            % Use i.e. for 'String' property of popupmenus or for
            % 'Callback' with options.
            %please note, that if you define your own callback function,
            %you need to take care of synchronizaiton with outputParameters
            %yourself, e.g. by calling obj.obj.updateGuiParameter(0,0,guihandle);
            %.position: relative position in a 4x11 grid. You can use
            %non-integer values
            %optional fields: Width, Height. In relative units.
            
            pard.t_Depth.object=struct('String','Depth','Style','text','Value',0);
            pard.t_Depth.position=[1,1.5];
            pard.t_Depth.Width=0.5;
            
            pard.t_Shape.object=struct('String','Shape','Style','text','Value',0);
            pard.t_Shape.position=[1,2];
            pard.t_Shape.Width=0.5;
            
            pard.t_Cap.object=struct('String','Cap','Style','text','Value',0);
            pard.t_Cap.position=[2,1];
            pard.t_Cap.Width=1;
            
            pard.t_Bottom.object=struct('String','Bottom','Style','text','Value',0);
            pard.t_Bottom.position=[3,1];
            pard.t_Bottom.Width=1;
            
            pard.capDepth.object=struct('String','25','Style','edit','Value',0);
            pard.capDepth.position=[2,1.5];
            pard.capDepth.Width=0.5;
            
            pard.bottomDepth.object=struct('String','110','Style','edit','Value',0);
            pard.bottomDepth.position=[3,1.5];
            pard.bottomDepth.Width=0.5;
            
            pard.capShape.object=struct('String','ellipse','Style','edit','Value',0);
            pard.capShape.position=[2,2];
            pard.capShape.Width=0.5;
            
            pard.bottomShape.object=struct('String','bar','Style','edit','Value',0);
            pard.bottomShape.position=[3,2];
            pard.bottomShape.Width=0.5;
            
            pard.t_radius.object=struct('String','Radius','Style','text','Value',0);
            pard.t_radius.position=[4,1];
            pard.t_radius.Width=1;
            
            pard.radius.object=struct('String','25','Style','edit','Value',0);
            pard.radius.position=[4,1.75];
            pard.radius.Width=0.5;
            
            pard.t_proteinID.object=struct('String','ID','Style','text','Value',0);
            pard.t_proteinID.position=[6,1];
            pard.t_proteinID.Width=0.5;
            
            pard.proteinID.object=struct('String','00000','Style','edit','Value',0);
            pard.proteinID.position=[6,1.5];
            pard.proteinID.Width=0.5;
            
            pard.t_proteinName.object=struct('String','Name','Style','text','Value',0);
            pard.t_proteinName.position=[6,2];
            pard.t_proteinName.Width=0.5;
            
            pard.proteinName.object=struct('String','slaaa','Style','edit','Value',0);
            pard.proteinName.position=[6,2.5];
            pard.proteinName.Width=0.5;
            
            pard.t_velocityEPM.object=struct('String','Velocity','Style','text','Value',0);
            pard.t_velocityEPM.position=[7,1];
            pard.t_velocityEPM.Width=0.5;
            
            pard.velocityEPM.object=struct('String','9','Style','edit','Value',0);
            pard.velocityEPM.position=[7,1.5];
            pard.velocityEPM.Width=0.5;
            
            pard.t_zmt0.object=struct('String','Origin','Style','text','Value',0);
            pard.t_zmt0.position=[7,2];
            pard.t_zmt0.Width=0.5;
            
            pard.zmt0.object=struct('String','-135','Style','edit','Value',0);
            pard.zmt0.position=[7,2.5];
            pard.zmt0.Width=0.5;
            
            pard.t_timeParPath.object=struct('String','Time Parameters','Style','text');
            pard.t_timeParPath.position=[8,1];
            pard.t_timeParPath.Width=1;
            
            pard.timeParPath.object=struct('String','/','Style','edit');
            pard.timeParPath.position=[8,2:3.5];
            pard.timeParPath.Width=1;
            
            pard.file_button.object=struct('String','Choose file','Style','pushbutton','Callback',{{@load_callback,obj}});
            pard.file_button.position=[8,3];
            pard.file_button.Width=1;   
            
            pard.saveAs.object=struct('String','Save created object(s) as...','Style','text');
            pard.saveAs.position = [9,1];
            pard.saveAs.Width = 1.5;
            
            pard.save_button.object=struct('String','Save','Style','pushbutton' ,'Callback',{{@exportProteinES_callback,obj}});
            pard.save_button.position = [9,2.5];
            pard.save_button.Width = 0.5;
            %pard.guiobject.TooltipString=sprintf('you can define a tooltip string. \n This tip is displayed when you hover the mouse on the control');
            
            pard.t_load.object=struct('String','Load created object(s) as...','Style','text');
            pard.t_load.position = [10,1];
            pard.t_load.Width = 1.5;
            
            pard.load_button.object=struct('String','Load','Style','pushbutton' ,'Callback',{{@importProteinES_callback,obj}});
            pard.load_button.position = [10,2.5];
            pard.load_button.Width = 0.5;
            
            %provide a description and name in the field: plugininfo.
            pard.plugininfo.name='proteinES_ctrl';
            pard.plugininfo.description='write a description for your plugin';
            pard.plugininfo.type='ProcessorPlugin'; %type of plugin. Currently: ProcessorPlugin, WorkflowModule, WorkflowFitter, Renderer, LoaderPlugin, SaverPlugin, ROI_Analyze, ROI_Evaluate,WorkflowIntensity
             %define global paramters that are used (passed on in p)
            pard.inputParameters={'sr_pixrec'};
            
            %define which properties and uicontrols are written directly to the global parameter
            %structure.
            pard.outputParameters={'guiobject'};
            
            %define which uicontrols are synchronized with each other. Define what to synchronize:
            %String and/or Value. You can also define a function which is
            %called after the global parameter has been changed. this
            %function is only called, if the parameter is changed by a
            %different class.
            pard.syncParameters={{'globalParameterName','t_radius',{'String','Value'},{@aftersync_callback,obj}}};
            

        end
    end
end

function exportProteinES_callback(a, b, obj)
    proteinESObj = obj.getPar('proteinES');
    if class(proteinESObj)== "struct"
        [f,p]=uiputfile({'*.mat'},'Name the proteinES file','/');
        if ~f
            return
        end
        save([p f], 'proteinESObj');
    else
        warning('There is no generated proteinES object.')
    end
end

function importProteinES_callback(a, b, obj)
    
    [f,p]=uigetfile({'*.mat'},'Choose a proteinES file','/');
    if ~f
        return
    end
    load([p f], 'proteinESObj');
    obj.setPar('proteinES', proteinESObj);

end

function callbackfunction(uiobject,data,obj,extradata)
disp('callback')
end

function aftersync_callback(obj)
disp('aftersync')
end

function load_callback(a,b,obj)
f=obj.getSingleGuiParameter('timeParPath');
[f,p]=uigetfile({'*.txt'},'Choose timePar file',f);
if ~f
    return
end
obj.setGuiParameters(struct('timeParPath',[p f]));
end