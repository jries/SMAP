classdef WorkflowData %<handle
    properties %(Access=private)
        data %=workflow.myhandle;
        frame
        ID
        numberInTag=1;
        eof=false;
    end
    methods
%         function set(obj,data)
%             obj.data=data;
% %             obj.data.set(data);
%         end
%         function d=get(obj)
%             d=obj.data;
% %             d=obj.data.get;
%         end
%         function out=copy(obj)
%             out=obj;
% %             out=interfaces.WorkflowData;
% %             out.frame=obj.frame;
% %             out.ID=obj.ID;
% %             out.numberInTag=obj.numberInTag;
% %             out.eof=obj.eof;
% %             out.data=[]; 
%         end
%         function clear(obj)
%             obj.data=[];
%             obj.frame=[];
%             obj.ID=[];
%             obj.numberInTag=0;
%             obj.eof=false;
%             
%         end
    end
end