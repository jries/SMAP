classdef Channeldrift4Pi<interfaces.WorkflowModule
    properties
        guidefh
    end
    methods
        function obj=Channeldrift4Pi(varargin)
                obj@interfaces.WorkflowModule(varargin{:});
                obj.inputChannels=1;
                obj.isstartmodule=true;
        end   
        function initGui(obj)

        end
        function prerun(obj,p)
        end
        function run(obj,data,p)
            wf=obj.parent;
            peakcombiner=wf.module('PeakCombiner');
            framecorrection=peakcombiner.guihandles.framecorrection.Value;

            if obj.getPar('loc_preview') 
                obj.output(data,1)
                return
            end

            % first fit of few frames
            loader=wf.module('Loader4PiMat');
            framestop=loader.guihandles.framestop.String;
            loader.guihandles.framestop.String=num2str(p.frameblock);

            locsaver=wf.module('LocSaver');
            locsaver.saveon=false;

            peakfinder=wf.module('PeakFinder');
            cutoff=peakfinder.guihandles.cutoffvalue.String;
            peakfinder.guihandles.cutoffvalue.String='3'; %to GUI?

            fitter=wf.module('MLE_4Pi');
            oldlink=fitter.guihandles.globaltable.Data;
            fitter.guihandles.globaltable.Data{1,1}=false;
            fitter.guihandles.globaltable.Data{2,1}=false;    
           
            peakcombiner.guihandles.framecorrection.Value=false;
            obj.setPar('loc_4Pichanneldrift',[]);
            
            obj.clearinitialize;
            obj.initialize; %calls prerun of following modules
            obj.output(data,1)

            channeldrift1=getchanneldrift4Pi(obj.locData.loc,p.frameblock+1);
            obj.setPar('loc_4Pichanneldrift',channeldrift1);
            
            %next iteration
            if framecorrection
                loader.guihandles.framestop.String=framestop;
            end
            peakcombiner.guihandles.framecorrection.Value=true;
            
           
            obj.clearinitialize;
            obj.initialize; %calls prerun of following modules
            obj.output(data,1)
            
            channeldrift=getchanneldrift4Pi(obj.locData.loc,p.frameblock,channeldrift1.dx,channeldrift1.dy);

            %add previous channel drift
            obj.setPar('loc_4Pichanneldrift',channeldrift);
%             peakcombiner.guihandles.framecorrection.Value=true;
            peakfinder.guihandles.cutoffvalue.String=cutoff; %to GUI?
            fitter.guihandles.globaltable.Data=oldlink;
            locsaver.saveon=true;
            loader.guihandles.framestop.String=framestop;

            obj.clearinitialize;
            obj.initialize; %calls prerun of following modules
            obj.output(data,1)

            peakcombiner.guihandles.framecorrection.Value=framecorrection;

        end
        function pard=guidef(obj)
            pard.frameblockt.object=struct('Style','text','String','window (frames)');
            pard.frameblockt.position=[1,1];
            pard.frameblockt.Width=2;
            pard.frameblock.object=struct('Style','edit','String','3000');
            pard.frameblock.position=[1,3];
            pard.frameblock.Width=1;
            pard.plugininfo.type='WorkflowModule'; 
            pard.plugininfo.description='runs the 4Pi fitting two times to determine channel drift';
        end
    end
end

