classdef Workflow_Template<interfaces.WorkflowModule
    % WOEKDLOW_TEMPLATE Summary of this plugin goes here
    % put a description of your plugin here.
        %replace Workflow_Template by filename  
    properties
        %here you can define properties. Can be used as storage for
        %structures initialized in prerun
    end
    methods
        function obj=Workflow_Template(varargin) %replace Workflow_Template by filename 
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=2;  %set the  number of input Channels. Usually it is 1.
        end
        function pard=guidef(obj)
            pard=[];
%             see Plugin_Template.m 
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
            obj.setInputChannels(2,'frame'); %define how to synchronize the channels. Either use 'frame' or 'ID'. 
            
        end
        function prerun(obj,p)
            %called before the workflow is run. Here you can initialize
            %parameters, store them in properties etc.
           
        end
        function output=run(obj,datin,p)
            %p: structure with fields corresponding ot inputParameters and
            %uicontrols of GUI.
            
            %datin: structure which contains the data passed between
            %modules.
            % if obj.inputChannels>1 it is a cell of these structures with
            % the lenght of obj.inputChannels
            %data<interfaces.WorkflowData
            %first module should produce this. E.g.
            dat=interfaces.WorkflowData;
            dat.data=rand(100); %you can store any kind of data in data. E.g. Images or structures with localizations
            dat.frame=1; %important if this is used for synchronization. You can use instead ID.
            dat.numberInTag=1; %says there is only one block of data per tag.
            dat.eof=false; %last data block shoudl have eof=true. This tells all workflow modules that all data has been there.
            
            %if this is not the first module, you can now process
            %datin.data
            
            % often it makes sense to see if there is data, since some
            % moudles might output empty data for the last data block,
            % others might not.
            if ~isempty(datin{1}.data)
                datnew=datin{1}.data.*2; %more than 1 inputChannel
                datnew=datin.data.*2; 
            else
                datnew=[];
            end
                %to output data, create output data structure. Either manually
                %as before, then set datout.dat. Or simply copy input.
                datout=datin{1};
                datout=datin;

                datout.data=datnew;

                %there are two options for outputting data to the next module.
                %if you have a single output, you can define
                output=datout{1};

                %Alternatively, you can manually output all data:
                output=[]; 
                obj.output(datout,outputChannel); %you can have more than 1 output channel. 
                obj.output(datout); %outputs on channel 1.
            
            
            %some times you want to clean up after everything is done, eg.
            %save data. Here exampel for 1 channel
            if datin.eof
                %do something
                
                %don't forget to output the data, in case this is not done
                %before in the code.
                obj.output(datin);
            end
        end
    end
end

