classdef Arrayfeeder<interfaces.WorkflowModule
    % Arrayfeeder
  
    properties   
        imstack
    end
    methods
        function obj=Arrayfeeder(varargin)
            obj@interfaces.WorkflowModule(varargin{:})
            obj.inputChannels=1;
            obj.isstartmodule=true;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end

        function run(obj,data,p)
            for k=1:size(obj.imstack,3)
                datout=interfaces.WorkflowData;
                datout.data=obj.imstack(:,:,k);

                datout.frame=k;
                datout.ID=k;

                obj.output(datout)
            end
            datout.data=[];
            datout.frame=k+1;
            datout.ID=k+1;
            datout.eof=true;
            obj.output(datout)
        end
    end
end



function pard=guidef(obj)

pard.plugininfo.type='WorkflowModule'; 
t1='feeds image';
pard.plugininfo.description=t1;
end