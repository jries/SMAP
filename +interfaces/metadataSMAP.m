classdef metadataSMAP<handle
    %metadata class for smap. mirrors the locData.files.file.info
    
    properties
        Width=512;
        Height=512;
%         roi_metadata=[];
        roi=[]; 
        
        camerainfo=[];
        allmetadata=[];
        exposure=1;
        emgain=1;
        EMon=false;
        conversion=1;
        offset=100;
%         pixsize=0.1;
        cam_pixelsize_um=0.1;
        timediff=30;
        comment='';
        numberOfFrames=0;
        pix2phot=[];
        basefile;
        assigned
        imagefile
        EMmirror=false; %used by Tiff loader workflow
    end
    
    methods
        function  obj=metadataSMAP
            fn=(properties(obj));
            for k=1:length(fn)
                obj.assigned.(fn{k})=true;
            end
%             obj.assigned.roi=true;
        end
        function roi=get.roi(obj)
            if isempty(obj.roi)
                roi=[0 0 obj.Width obj.Height];
%                 obj.assigned.roi=true;
            else
                roi=obj.roi;
                
            end
        end
        function cv=get.pix2phot(obj)
            if isempty(obj.pix2phot)
                if ~isnan(obj.EMon)&&obj.EMon
                    cv=obj.conversion/obj.emgain;
                else
                    cv=obj.conversion;
                end
                obj.pix2phot=cv;
%                 obj.assigned.pix2phot=true;
            else
                cv=obj.pix2phot;
            end
        end
        function set.roi(obj,roi)
            if length(roi)<4
                return
            end
            obj.roi=roi;
            obj.Width=round(roi(3));
            obj.Height=round(roi(4));
            obj.assigned.roi=true;
            obj.assigned.Width=true;
            obj.assigned.Height=true;
        end
    end
    
end

