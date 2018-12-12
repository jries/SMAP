classdef LocTransformN0<handle
    %coordinates: right now it is totally confuse (coordinate vs particle)
    %coordinates(particle, xyz)
    properties
        transform2Reference
        transform2Target
        info={[]};
        transformZ2Reference
        transformZ2Target
        unit='nm'; %or pixel
%         cam_pixnm={[100 100],[100 100]};  %In future: allow pixel size for every channel? e.g. for multi-camera setup
        channels=2;
%         mirror
    end
    
    methods
        function obj=LocTransformN0(varargin)
            obj.transform2Reference{1}=affine2d(eye(3)); %initialize channel 1 to reference
            obj.transform2Target{1}=affine2d(eye(3));
            obj.transformZ2Reference{1}=affine3d(eye(4)); %initialize channel 1 to reference
            obj.transformZ2Target{1}=affine3d(eye(4));
            info{1}.xrange=[-inf inf];
            info{1}.yrange=[-inf inf];
        end
        function setTransform(obj,channel,varargin)
            % type: any of matlab transformations
            % mirror:  dimension along which to mirror: 0, 1, 2 , [1 2] or
            % 3
            %unit: nm pixel
            % also possibility to pass on structure
            if isempty(channel)
                channel=1:length(obj.info);
            end
            obj.channels=max(obj.channels,max(channel));
            
%             properties={'xrange','yrange','type','parameter','mirror','unit','channels','cam_pixnm'};
            propertiesch={'xrange','yrange','type','parameter','mirror','cam_pixnm'};
            if isstruct(varargin{1})   
                p=varargin{1};
                if max(channel)>length(obj.info)
                    obj.info{max(channel)}=[];
                end
                for c=1:length(channel)
                    obj.info{channel(c)}=copyfields(obj.info{channel(c)},p,propertiesch);
                end
                if isfield(p,'unit')
                    obj.unit=p.unit;
                end
                if isfield(p,'channels')
                    obj.channels=p.channels;
                end
%                 if isfield(p,'cam_pixnm')
%                     obj.cam_pixnm{channel}=p.cam_pixnm;
%                 end
            else
                  
            for k=1:2:length(varargin)
                if any(strcmp(propertiesch,varargin{k}))
                    for c=1:length(channel)
                        obj.info{channel(c)}.(varargin{k})=varargin{k+1};
                        
                    end
                else 
                    disp([varargin{k} ' is not a proper parameter']);
                end
            end
            ind=find(strcmp(varargin,'unit'));
            if ~isempty(ind)
                obj.unit=varargin{ind+1};
            end
            ind=find(strcmp(varargin,'channels'));
            if ~isempty(ind)
                obj.channels=varargin{ind+1};
            end
%             ind=find(strcmp(varargin,'cam_pixnm'));
%             if ~isempty(ind)
%                 for ch=channel
%                 obj.cam_pixnm{ch}=varargin{ind+1};
%                 end
%             end   
            end
        end
        function findTransform(obj,channel,coordreference,coordtarget,type,parameter)
            %XXX make sure reference and target are not exchanged 
            %matrix becomes singular for too large coordinates. Use um
            %instead of nm
            %use units as specified
            
            if nargin<5
                th=obj.info{channel};
                if isfield(th,'type') && ~isempty(th.type)
                    type=th.type;
                    if isfield(th,'parameter')
                        parameter=th.parameter;
                    end
                else
                    parameter=[];
                    for k=1:length(obj.info) %look for other channels if any type is specified, use the same
                         th=obj.info{k};
                        if isfield(th,'type') && ~isempty(th.type)
                            type=th.type;
                            if isfield(th,'parameter')
                                parameter=th.parameter;
                            end
                            break
                        end
                        
                    end
                    if isempty(parameter)
                        warning('no transformation type specified');
                        return
                    end
                end
            else
                obj.info{channel}.type=type;
                if nargin>5
                    obj.info{channel}.parameter=parameter;
                end
            end
            coordreference=double(coordreference/1000);
            coordtarget=double(coordtarget/1000);
            switch type
                case {'lwm','polynomial'}
                    obj.transform2Target{channel} = fitgeotrans(coordtarget(:,1:2),coordreference(:,1:2),type,parameter);
                    obj.transform2Reference{channel} = fitgeotrans( coordreference(:,1:2),coordtarget(:,1:2),type,parameter);
                otherwise             
                    obj.transform2Target{channel}= fitgeotrans(coordtarget(:,1:2),coordreference(:,1:2),type);
                    obj.transform2Reference{channel}= fitgeotrans(coordreference(:,1:2),coordtarget(:,1:2),type);
            end
            if size(coordreference,2)>2 %3D data set: also do z-transform
                coordtargettransformed=horzcat(obj.transformToReference(channel,coordtarget(:,1:2)),coordtarget(:,3));
                    %only affine3d possible      
                obj.transformZ2Target{channel}= findAffineTransformZ(coordtargettransformed,coordreference(:,3));
                obj.transformZ2Reference{channel}= findAffineTransformZ(coordreference,coordtarget(:,3));
            end
        end
        
       function co=convertcoordinates(obj, ci,unitref,unittar,channel)
           if nargin<5
               channel=1;
           end
           if ~strcmp(unitref,unittar)
                switch unittar
                    case {'pixels','pixel'}
                        cf=1./obj.info{channel}.cam_pixnm;
%                         cf=1./obj.cam_pixnm{channel};
                    case 'nm'
                        cf=obj.info{channel}.cam_pixnm;        
                end
    %             co(1,:)=ci(1,:)*cf(1);
    %             co(2,:)=ci(2,:)*cf(end);
                co(:,1)=ci(:,1)*cf(1);
                co(:,2)=ci(:,2)*cf(end);
            end
        end
        
        function co=transformToReference(obj,channel,ci,unit)
            if nargin>3 %unit specified
                ci=obj.convertcoordinates(ci,unit,obj.unit,channel);
            end
            ci=ci/1000;
            co=transformPointsInverse(obj.transform2Reference{channel},ci(:,1:2)); %inverse of inverse is forward          
            if size(ci,2)>2 %z coordinates present               
                 X=transformPointsInverse(obj.transformZ2Reference{channel},horzcat(co,ci(:,3)));
                 co(:,3)=X(:,3);
            end
             if nargin>3,co=obj.convertcoordinates(co,obj.unit,unit,1);end  %now in Reference channel (1), use this pixel size
             co=co*1000; %back to nm   
        end
        function co=transformToTarget(obj,channel,ci,unit)
            if nargin>3 %unit specified
                ci=obj.convertcoordinates(ci,unit,obj.unit,1); %start with reference coordinates ch 1
            end
            ci=ci/1000;
            co=transformPointsInverse(obj.transform2Target{channel},ci(:,1:2)); %inverse of inverse is forward          
            if size(ci,2)>2 && length(obj.transformZ2Target)>=channel %z coordinates present               
                 X=transformPointsInverse(obj.transformZ2Target{channel},horzcat(co,ci(:,3)));
                 co(:,3)=X(:,3);
            end
            if nargin>3,co=obj.convertcoordinates(co,obj.unit,unit,channel);end %now at target channel
             co=co*1000; %back to nm   
        end  
        function co=transformToTargetAll(obj,varargin)
            ci=varargin{1};
            numch=length(obj.transform2Target);
      
            si=size(ci);si(3)=numch;
            co=zeros(si);
            for ch=1:numch
                co(:,:,ch)=obj.transformToTarget(ch,varargin{:});
            end
        end

        function imout=transformImageToTarget(obj,channel,image,cam_pixnm,roi)
            imout=transformImage(obj.transform2Target{channel},image,cam_pixnm,roi);
        end

        function imout=transformImageToReference(obj,image,cam_pixnm,roi)
            imout=transformImage(obj.transform2Reference{channel},image,cam_pixnm,roi);
        end  
        
        function ind=getPart(obj,channel,coordinates,unit)
            % obj, channel, coordinates, unit
            if nargin>3 %unit specified
                coordinates=obj.convertcoordinates(coordinates,unit,obj.unit,channel);
            end
            th=obj.info{channel};
            if ~isfield(th,'xrange')||isempty(th.xrange)||~isfield(th,'yrange')||isempty(th.yrange)
                disp('no range specified, getPart in LocTransformN returns all coordinates');
                ind=true(size(coordinates,1),1);
            else
                ind=coordinates(:,1)>th.xrange(1) & coordinates(:,1)<=th.xrange(2) & coordinates(:,2)>th.yrange(1)& coordinates(:,2)<th.yrange(2); 
            end
        end


        function makeAffine2d(obj,channel,A)
            tform=affine2d(A);
            obj.transform2Target=tform;
            obj.transform2Reference=invert(tform);
        end
        
        function out=mirror(obj,channel,format)
            %out= logical (x y vs channels), channel specified only xy.
            %Format: return string

            if nargin <2 || isempty(channel)
                out=false(obj.channels,2);
                for k=2:obj.channels
                    out(k,:)=mirrorchannel(obj, k);
                end
            else
                out=mirrorchannel(obj, channel);
            end
            if nargin>2 && strcmp(format,'str')
                out=convert2str(out);
            end
          
  
        end
        function out=mirrorchannel(obj, channel)
             ci=[0,0;1,1];
            
            co=obj.transformToTarget(channel,ci);
            dc=sign(co(2,:)-co(1,:));
            out=[false false];
            if dc(1)<0
                out(1)=true;
            end
            if dc(2)<0
                out(2)=true;
            end
        end
        
    end
end


function out=convert2str(in)
s=size(in);
for k=1:s(1)
    out{k}='none';
    if in(k,1)
        out{k}='left-right';
    end
    if in(k,2)
        out{k}='up-down';
    end
    if in(k,1) && in(k,2)
        out{k}='both';
    end
end
if s(1)==1
    out=out{1};
end
end
