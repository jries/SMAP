classdef SEsites<matlab.mixin.Copyable
    properties
        pos
        ID
        info
        annotation
        evaluation
        name
        image
        handles 
        sePar
        indList
    end
    methods
        function obj=SEsites(varargin)
            obj@matlab.mixin.Copyable;
            if nargin>0
                sein=varargin{1};
                fn=properties(sein);
                for k=1:length(fn)
                    obj.(fn{k})=sein.(fn{k});
                end
            end
        
            list.value=1;
            list.string={'empty'};
            obj.annotation.list1=list;
            obj.annotation.list2=list;
            obj.annotation.list3=list;
            obj.annotation.list4=list;
            obj.annotation.comments='';
            obj.annotation.use=true;

            line.pos=zeros(2);line.value=0;line.angle=0;line.length=0;
            obj.annotation.line1=line;
            obj.annotation.line2=line;
            obj.annotation.rotationpos=line;
        end
        function setlineangle(obj,number,angledeg,len)
            angle=angledeg/180*pi;
            if nargin<4
                len=200;
            end
            len=len/2;
%             line.length=len;line.angle=angle;line.value=len;
%             line.
            posh=[obj.pos(1)-len*cos(angle),obj.pos(2)-len*sin(angle);obj.pos(1)+len*cos(angle),obj.pos(2)+len*sin(angle)];
            obj.setlinepos(number,posh);
%             obj.annotation.rotationpos=line;
        end
        function setlinepos(obj,number,x1,y1,x2,y2)
            switch number
                case 0 
                    field='rotationpos';
                case 1
                    field='line1';
                case 2
                    field='line2';
            end
            if nargin>5
                obj.annotation.(field).pos=[x1 y1;x2,y2]/1000;
            else
                obj.annotation.(field).pos=x1/1000;
            end
            pos=obj.annotation.(field).pos;
            len=sqrt(sum((pos(1,:)-pos(2,:)).^2))*1000;
            angle=pos2angle(pos);
             obj.annotation.(field).length=len;obj.annotation.(field).angle=angle;obj.annotation.(field).value=len;
        end
            
    end
end