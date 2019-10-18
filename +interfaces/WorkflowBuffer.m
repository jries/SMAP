classdef WorkflowBuffer<handle
    properties(Access=private)
        dat={};
        currenttag=1;
        bufferincrease=10;
        complete=[];
    end
    methods
        function add(obj,data,tag)
            lod=length(obj.dat);
            if lod<tag
                    obj.dat{tag}=data;
%                 obj.complete(lod+1:tag+obj.bufferincrease)=false;
%                  obj.complete(lod+1:tag+obj.bufferincrease)=false;
            else
            
            ld=length(obj.dat{tag});
 
            obj.dat{tag}(ld+1)=data;
            end
            if data.numberInTag>0 && length(obj.dat{tag})>=data.numberInTag
                obj.complete(end+1)=tag;
%                 obj.complete(tag)=true;
            end
            if tag>obj.currenttag
                if data.numberInTag==0&&length(obj.dat{obj.currenttag})>=1
                    obj.complete(end+1)=obj.currenttag;
%                     obj.complete(obj.currenttag)=true;
                end
                obj.currenttag=tag;    
            end
            if data.eof
                 obj.complete(end+1)=tag;
%                  obj.complete(tag)=true;
            end
  
        end
        function out=iscomplete(obj,tag)
            if nargin<2
                out=obj.complete;
            else
%                 out=obj.complete(tag);
                out=any(obj.complete==tag);
            end
        end
        function data=get(obj,tag)
            data=obj.dat{tag}(1:end);
%             obj.dat{tag};
%              obj.dat{tag}=interfaces.WorkflowData;
%             for k=1:length(obj.dat{tag}.data)
            obj.dat{tag}.data=[];
            
%             end
%             obj.complete(tag)=false;  
            obj.complete(obj.complete==tag)=[];
        end
    end
end