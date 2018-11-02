classdef ChangeCamPixelsize<interfaces.DialogProcessor
    methods
        function obj=ChangeCamPixelsize(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','CombineFiles');
            notify(obj.P,'backup4undo');

            pixold=obj.locData.files.file(p.file1.Value).info.cam_pixelsize_um;
            if length(pixold)==1
                pixold(2)=pixold(1);
            end       
            pixnew=[p.pixx p.pixy];
            factor=pixnew./pixold;
            
            ind=obj.locData.loc.filenumber==p.file1.Value;
            obj.locData.loc.xnm(ind)=obj.locData.loc.xnm(ind)*factor(1);
            obj.locData.loc.ynm(ind)=obj.locData.loc.ynm(ind)*factor(2);
            %update filelist
            obj.locData.files.file(p.file1.Value).info.cam_pixelsize_um=pixnew;
            for k=1:length(obj.locData.files.file(p.file1.Value).tif)
                obj.locData.files.file(p.file1.Value).tif(k).info.cam_pixelsize_um=pixnew;
            end
            obj.locData.regroup;
            initGuiAfterLoad(obj)
             
        end
        function initGui(obj)
           
            %synch filelist short to table
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function updatefiles(obj)
        end
    end
end
function fileselect(a,b,obj)
pixold=obj.locData.files.file(a.Value).info.cam_pixelsize_um;
pixnew=pixold;
if length(pixold)==1
    pixnew(2)=pixnew(1);
end

obj.setGuiParameters(struct('pixx',pixnew(1),'pixy',pixnew(2),'currentpix',pixold))
obj.guihandles.currentpix.String=num2str(pixold);
end
function factor_callback(a,b,obj)
p=obj.getGuiParameters;
Ynew=p.pixx*p.factorxy;
obj.setGuiParameters(struct('pixy',Ynew));
end
function pard=guidef(obj)

pard.file1.object=struct('Style','popupmenu','String','File','Callback',{{@fileselect,obj}});
pard.file1.position=[1,1];
pard.file1.object.TooltipString='choose localization file data set';
pard.file1.Width=3;

pard.t1.object=struct('Style','text','String','current pixelsize:');
pard.t1.position=[2,1];
pard.currentpix.object=struct('Style','text','String','selectFile');
pard.currentpix.position=[2,2];

pard.t2.object=struct('Style','text','String','New pixelsize X:');
pard.t2.position=[3,1];
pard.pixx.object=struct('Style','edit','String','??');
pard.pixx.position=[3,2];
pard.t3.object=struct('Style','text','String','Y:');
pard.t3.position=[3,3];
pard.t3.Width=0.25;
pard.pixy.object=struct('Style','edit','String','??');
pard.pixy.position=[3,3.25];

pard.t4.object=struct('Style','pushbutton','String','Factor Y=k*X:','Callback',{{@factor_callback,obj}});
pard.t4.position=[4,1];
pard.factorxy.object=struct('Style','edit','String','1');
pard.factorxy.position=[4,2];


pard.syncParameters={{'filelist_short','file1',{'String'}}};

pard.plugininfo.type='ProcessorPlugin';
end