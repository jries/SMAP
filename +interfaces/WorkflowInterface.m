classdef WorkflowInterface<handle
    properties
        outputChannels=1;
        outputModules
        inputChannels=1;  
        isstartmodule=false;
        initialized=false;
        UID
    end
    properties(Access=private) 
        inputData
        syncmode='frame'; %ID or frame. used for synching two inputs.
        runparameters;
    end
    methods
        function obj=WorkflowInterface
            s=char(java.rmi.server.UID());
            s(s=='-')=[];s(s=='.')=[];
            obj.UID=s;
        end
        function setInputModule(obj,inputchannel,module,outputchannel)
            if nargin<4
                outputchannel=1;
            end
            obj.inputChannels=max(obj.inputChannels,inputchannel);
            
            inputinfo.inputChannels=obj.inputChannels;
            inputinfo.inputData=obj.inputData;
            inputinfo.syncmode=obj.syncmode;
            
            module.addNextModule(obj,outputchannel,inputchannel,inputinfo);
            
        end
        function setInputChannels(obj,inputChannels,syncmode)
            if inputChannels>1
                obj.inputChannels=inputChannels;
                obj.inputData=interfaces.SyncBuffer(obj.inputChannels);
            end
            if nargin>2
                obj.syncmode=syncmode;
            end
        end
        function addNextModule(obj,module,outputchannel,inputchannel,inputinfo)        
            idx=length(obj.outputModules)+1;
            for k=1:length(obj.outputModules)
                if obj.outputModules(k).module==module&&obj.outputModules(k).outputchannel==outputchannel&&obj.outputModules(k).inputchannel==inputchannel
                    idx=k;
                    break
                end
            end
            obj.outputModules(idx).module=module;
            obj.outputModules(idx).outputchannel=outputchannel;
            obj.outputModules(idx).inputchannel=inputchannel;
            obj.outputModules(idx).inputinfo=inputinfo;        
        end
        function output(obj,data,outputchannel)
            if nargin<3
                outputchannel=1;
            end
            om=obj.outputModules;
            for k=1:length(om)
%                 om=obj.outputModules(k);
                if om(k).outputchannel==outputchannel
                    om(k).module.input(data,om(k).inputchannel,om(k).inputinfo);
                end
            end
        end
        
        function input(obj,data,inputchannel,inputinfo)
            p=obj.runparameters;
            inputChannels=inputinfo.inputChannels;
            if inputChannels==1
                output=obj.run(data,p); %simple: only one package, no sync needed
                if ~isempty(output)
                    obj.output(output);
                end
            else
                syncmode=inputinfo.syncmode;
                tag=data.(syncmode);
                 inputData=obj.inputData;
                inputData.add(data,tag,inputchannel);
                indc=inputData.iscomplete;
                for k=1:length(indc) 
                    dat=inputData.get(indc(k));
                    output=obj.run(dat,p); 
                    if ~isempty(output)
                        obj.output(output);
                    end
                end

            end
        end
        
        function reset(obj)
            if obj.inputChannels>1
                obj.inputData=interfaces.SyncBuffer(obj.inputChannels);
            end
        end
        function initialize(obj,a,b)
            obj.reset;
            obj.updateObjectParameters;
            p=obj.getAllParameters;
            obj.runparameters=p;
            if ~obj.initialized
                obj.prerun(p);
                obj.initialized=true;
            end
            for k=1:length(obj.outputModules)
                 obj.outputModules(k).module.initialize;
            end
        end
        function clearinitialize(obj,a,b)
            obj.initialized=false;
            for k=1:length(obj.outputModules)
                 obj.outputModules(k).module.clearinitialize;
            end
        end
        function prerun(obj,p)
        end
        
    end
end
