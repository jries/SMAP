classdef Make3Dshape<interfaces.DialogProcessor
    % PLUGIN_TEMPLATE Summary of this plugin goes here
    % put a description of your plugin here.
        %replace Plugin_Template by filename   
    properties
        %define class properties if needed
    end
    methods
        function obj=Make3Dshape(varargin)   %replace by filename        
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.showresults=true; %set true, if results are shown by default
            obj.history=false; %if set true, every time the plugin is called, its parameters are added to the history. Makes sense only if plugin changes the data
        end       
        function initGui(obj)
            %is called after the GUI, defined in guidef, is made. Here you
            %can add additional GUI components, e.g. some that are not
            %defined by simple uicontrols (e.g. uitable), set additional
            %synchronizations and callbacks.
            obj.guihandles.previewimage=axes(obj.handle,'Units','pixels');
            obj.guihandles.previewimage.Position=obj.guihandles.previewdummy.Position;
            obj.guihandles.previewimage.Position(4)=5*obj.guihandles.previewimage.Position(4);
            obj.guihandles.previewimage.Position(3)=0.6*obj.guihandles.previewimage.Position(3);
            visualSurCom(mkSurCom(obj.getSingleGuiParameter('inCap'), obj.getSingleGuiParameter('inBottom'), obj.getSingleGuiParameter('root'), obj.getSingleGuiParameter('inDia')),obj.guihandles.previewimage)
            hold(obj.guihandles.previewimage, 'on')
            visualSurCom(mkSurCom(obj.getSingleGuiParameter('outCap'), obj.getSingleGuiParameter('outBottom'), obj.getSingleGuiParameter('root'), obj.getSingleGuiParameter('outDia')),obj.guihandles.previewimage)
            hold(obj.guihandles.previewimage, 'off')
            axis(obj.guihandles.previewimage,'on')
            %axis(obj.guihandles.previewimage,'equal')
            axis(obj.guihandles.previewimage,[-obj.getSingleGuiParameter('outDia')/2 obj.getSingleGuiParameter('outDia')/2 0 (obj.getSingleGuiParameter('outCap')+obj.getSingleGuiParameter('outBottom')+obj.getSingleGuiParameter('root'))])
          
        end
        function out=run(obj,p)
            l = pointsInDefSpace(p);

            %Here you implement the functionality
            %p contains:
            %global parameters defined in pard.inputParameters
            %uicontrol parameters, already parsed (e.g. edit with numerical
            %value returns a double, not a text. popupmenu etc return three
            %fields: String (all entries), Value (selected item), selection
            %(string of selected item).
            
            %out=[] possible, but out needs to be defined. Out is stored
            %and accessible later or by other plugins
            
            %Only if you modify the localization date it makes sense to
            %enable Undo. This can be done with the follwoing lines:
            obj.setPar('undoModule','PluginName');
            notify(obj.P,'backup4undo');

            %Get access to the localization data via obj.locData.getlocs (see description there):
            %example:
            locs=obj.locData.getloc({'xnm','frame','ynm'},'position','roi','layer',1);
            %locs.xnm, locs.frame, locs.ynm contain data from all
            %localizations in roi, or if no roi is defined in field of
            %view. But only those, displayed (and filtered) in layer 1.
            
            %create some output:
            axis1=obj.initaxis('axisname');
            plot(locs.xnm,locs.ynm,'.','Parent',axis1);
            
            %set a global parameter
            obj.setPar('locFields',fieldnames(obj.locData.loc));
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
            
            pard.t_inner.object=struct('String','Inner cylinder','Style','text');
            pard.t_inner.position=[1,2];
            pard.t_inner.Width=1;
            
            pard.t_outer.object=struct('String','Outer cylinder','Style','text');
            pard.t_outer.position=[1,3];
            pard.t_outer.Width=1;
            
            pard.t_cap.object=struct('String','Cap depth','Style','text');
            pard.t_cap.position=[2,1];
            pard.t_cap.Width=1.5;
            
            pard.t_bottom.object=struct('String','Bottom depth','Style','text');
            pard.t_bottom.position=[3,1];
            pard.t_bottom.Width=1.5;
            
            pard.t_diameter.object=struct('String','Diameter','Style','text');
            pard.t_diameter.position=[4,1];
            pard.t_diameter.Width=1.5;
            
            pard.t_root.object=struct('String','Root','Style','text');
            pard.t_root.position=[5,1];
            pard.t_root.Width=1.5;
            
            pard.inCap.object=struct('String',20,'Style','edit');
            pard.inCap.position=[2,2];
            pard.inCap.Width=.5;
            
            pard.outCap.object=struct('String',150,'Style','edit');
            pard.outCap.position=[2,3];
            pard.outCap.Width=.5;
            
            pard.inBottom.object=struct('String',80,'Style','edit');
            pard.inBottom.position=[3,2];
            pard.inBottom.Width=.5;
                     
            pard.outBottom.object=struct('String',0,'Style','edit');
            pard.outBottom.position=[3,3];
            pard.outBottom.Width=.5;
            
            pard.inDia.object=struct('String',50,'Style','edit');
            pard.inDia.position=[4,2];
            pard.inDia.Width=.5;
                     
            pard.outDia.object=struct('String',120,'Style','edit');
            pard.outDia.position=[4,3];
            pard.outDia.Width=.5;
            
            pard.root.object=struct('String',0,'Style','edit');
            pard.root.position=[5,2.5];
            pard.root.Width=.5;
            
            pard.preview.object=struct('String','Preview','Style','pushbutton','Callback',{{@callbackfunction,obj}});
            pard.preview.position=[6,2.5];
            pard.preview.Width=0.5;
            
            pard.t_numMol.object=struct('String','#Molecular','Style','text');
            pard.t_numMol.position=[7,1];
            pard.t_numMol.Width=0.5;
            
            pard.numMol.object=struct('String',150,'Style','edit');
            pard.numMol.position=[7,2];
            pard.numMol.Width=0.5;
            pard.numMol.TooltipString=sprintf('The mean of molecular number in the defined space.');
            
            pard.t_viewType.object=struct('String','View type','Style','text');
            pard.t_viewType.position=[7,1];
            pard.t_viewType.Width=0.5;
            
            pard.viewType.object=struct('String',{{'side','top'}},'Style','popupmenu');
            pard.viewType.position=[7,2];
            pard.viewType.Width=0.5;
            
            pard.t_size.object=struct('String','Size of the image','Style','text');
            pard.t_size.position=[7,3];
            pard.t_size.Width=1;
            
            pard.size.object=struct('String', 120, 'Style','edit');
            pard.size.position=[7,4];
            pard.size.Width=0.5;
            
            pard.t_imgPath.object=struct('String','as','Style','text');
            pard.t_imgPath.position=[9,1];
            pard.t_imgPath.Width=1;
            
            pard.t_folderPath.object=struct('String','Save to','Style','text');
            pard.t_folderPath.position=[8,1];
            pard.t_folderPath.Width=1;
            
            pard.folderPath.object=struct('String','/','Style','edit');
            pard.folderPath.position=[8,2:3];
            pard.folderPath.Width=1;
            
            pard.folder_button.object=struct('String','Choose folder','Style','pushbutton','Callback',{{@saveTo_callback,obj}});
            pard.folder_button.position=[8,3];
            pard.folder_button.Width=1;        
            
            pard.imgPath.object=struct('String', 'my_image.png', 'Style','edit');
            pard.imgPath.position=[9,2:3];
            pard.imgPath.Width=2;     
            
            pard.previewdummy.object=struct('String', ' ', 'Style','text');
            pard.previewdummy.position=[5,3.7];
            pard.previewdummy.Width=2;   
            
            %provide a description and name in the field: plugininfo.
            pard.plugininfo.name='Make3Dshape';
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
            pard.syncParameters={{'globalParameterName','inCap',{'String','Value'},{@aftersync_callback,obj}}};
            

        end
    end
end



function  unitSur = mkSurCom(capDepth, bottomDepth, root, diameter)
    %% make a surface component for making a cylinder
    if capDepth > 0
        % Create a parabola cap
        a = aFinder(capDepth, diameter);
        pY = 0:1:capDepth;
        pX1 = ((pY-capDepth)/a).^(1/2);
        rmBottom = 1;
    else
        pX1=[];
        rmBottom = 0;
    end
    
    if bottomDepth > 0
        % add a bottom column
        pY = 0:1:bottomDepth;
        pX2 = repelem(diameter/2, numel(pY));
    else
        pX2=[];
        rmBottom = 0;
    end
    
    pX = [pX2(1:(numel(pX2)-rmBottom)) pX1];
    
    % set the root
    unitSur = [];
    unitSur.main = pX;
    unitSur.root = root;
    unitSur.mainDepth = capDepth + bottomDepth;
end

function a = aFinder(depth, dia)
%% This is for the determination of coefficient "a" for a Parabola, using depth (y-intercept) and diameter (x-intercept) as inputs
    Y = 0:1:depth;
    a = -depth/(dia/2).^2;
end

function  fig = visualSurCom(unitSur,h)
    reflection = [-unitSur.main; 0:unitSur.mainDepth];
    surfaceCurve = [reflection(:,end:-1:1), [unitSur.main; 0:unitSur.mainDepth]];
    plot(h,surfaceCurve(1,:), surfaceCurve(2,:))
end

function  [X, Y, Z] = mkCy(unitSur)
    [X,Y,Z] = cylinder(unitSur.main, 40); % create a unit cylinder (the range of z is from 0 to 1)
    X = reshape(X,[numel(X),1]);
    Y = reshape(Y,[numel(Y),1]);
    Z = reshape(Z,[numel(Z),1]);
    
    Z = Z * unitSur.mainDepth;
    Z = Z + unitSur.root;
end

function callbackfunction(uiobject,data,obj)
    visualSurCom(mkSurCom(obj.getSingleGuiParameter('inCap'), obj.getSingleGuiParameter('inBottom'), obj.getSingleGuiParameter('root'), obj.getSingleGuiParameter('inDia')),obj.guihandles.previewimage)
    hold(obj.guihandles.previewimage, 'on')
    visualSurCom(mkSurCom(obj.getSingleGuiParameter('outCap'), obj.getSingleGuiParameter('outBottom'), obj.getSingleGuiParameter('root'), obj.getSingleGuiParameter('outDia')),obj.guihandles.previewimage)
    hold(obj.guihandles.previewimage, 'off')
    axis(obj.guihandles.previewimage,'on')
    %axis(obj.guihandles.previewimage,'equal')
    axis(obj.guihandles.previewimage,[-obj.getSingleGuiParameter('outDia')/2 obj.getSingleGuiParameter('outDia')/2 0 (obj.getSingleGuiParameter('outCap')+obj.getSingleGuiParameter('outBottom')+obj.getSingleGuiParameter('root'))])
disp('callback')
end

function saveTo_callback(a,b,obj)
f=obj.getSingleGuiParameter('folderPath');
f=uigetdir(f,'Choose folder for saving the result');
if ~f
    return
end
obj.setGuiParameters(struct('folderPath',f));
setvisibility(obj)
end

function aftersync_callback(obj)

disp('aftersync')
end