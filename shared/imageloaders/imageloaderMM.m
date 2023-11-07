classdef imageloaderMM<interfaces.imageloaderSMAP
    %imageloaderMM image loader for micromanager  tiff stack files
    %   Detailed explanation goes here
    
    properties
        reader
        readoutimgtags={};
        imtags
        init
    end
    
    methods
        function obj=imageloaderMM(varargin)
            obj@interfaces.imageloaderSMAP(varargin{:});
        end
        function openi(obj,file)
            initMM(obj);
            try
                obj.reader.close;
            end
            if ~exist(file,'file')
                return
            end
            obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(file), false, [], false, false, true);
            obj.file=file;
%             obj.reader=bfGetReader(file);
            md=obj.getmetadata;
            [p,f]=fileparts(file);
            obj.metadata.basefile=[p ];
            initimagetags(obj)
            
        end
        function image=getimagei(obj,frame)
            image=readstack(obj,frame);
        end
        function prefit(obj)
            initimagetags(obj)
        end
        function closei(obj)
            obj.reader.close
            
%             clear(obj.reader)
        end
        
        function image=getimageonline(obj,number)
            image=obj.getimage(number);
            if isempty(image)&&obj.onlineAnalysis 
                    disp('wait')
                    obj.reader.close;
%                     delete(obj.reader)
                    pause(obj.waittime*2)
                    obj.reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(obj.file), false, [], false, false, true);
                    image=obj.getimage(number);
            end
        end
        
        function allmd=getmetadatatagsi(obj)
            img=obj.reader;
            imgmetadata=img.getImageTags(0,0,0,0);
            summarymetadata=img.getSummaryMetadata;
            
            allmd=gethashtable(imgmetadata);
            alls=gethashtable(summarymetadata);
            try
            comments=char(img.getDisplayAndComments.get('Comments'));
            allmd(end+1,:)={'Comments direct',comments};
            catch
            end
            %direct
            try
            troi=textscan(imgmetadata.get('ROI'),'%d','delimiter','-');
            %XXXXX
            roih=troi{:}';
%             roih(1)=512-roih(1)-roih(3);
            allmd(end+1,:)={'ROI direct',num2str(roih)};
            catch err
            end
            sl=0;fr=0;po=0;
            try
                sl=summarymetadata.get('Slices');
            end
            try
                fr=summarymetadata.get('Frames');
            end
            try
                po=summarymetadata.get('Positions');
            end            
            possibleframes=[img.lastAcquiredFrame,sl,fr,po];
            framesd=min(possibleframes(possibleframes>100));
            if isempty(framesd)
                framesd=max(possibleframes);
            end
%             framesd=max([img.lastAcquiredFrame,summarymetadata.get('Slices'),summarymetadata.get('Frames'),summarymetadata.get('Positions')]);
            allmd(end+1,:)={'frames direct',num2str(framesd)};
            
            allmd=vertcat(allmd,alls);
            obj.allmetadatatags=allmd;
        end
        
    end
    
end





function image=readstack(obj,imagenumber)
img=obj.reader.getImage(0,0,imagenumber-1,0);
imgpos={0,0,imagenumber-1,0};
if isempty(img)
    img=obj.reader.getImage(0,imagenumber-1,0,0);
    imgpos={0,imagenumber-1,0,0};
end
if isempty(img)
    img=obj.reader.getImage(0,0,0,imagenumber-1);
    imgpos={0,0,0,imagenumber-1};
end
if isempty(img)
    image=[];
    return
end
image=img.pix;
if isempty(image)
    return
end
% if numel(image)==obj.metadata.Width*obj.metadata.Height
    image=reshape(image,obj.metadata.Width,obj.metadata.Height)';
    if isa(image,'int16')
        image2=uint16(image);
        ind=image<0;
%         image2(ind)=image(ind)+2^16;
         image2(ind)=2^16-uint16(-image(ind));
        image=image2;
    end

if ~isempty(obj.readoutimgtags)
    imgmeta=obj.reader.getImageTags(imgpos{:});
    if obj.init
        obj.imtags=zeros(length(obj.readoutimgtags),obj.metadata.numberOfFrames);
%          for k=1:length(obj.readoutimgtags)
%                 if ischar(imgmeta.get(obj.readoutimgtags{k}))
%                     obj.imtags{k}=strings(1,obj.metadata.numberOfFrames);
%                 else
%                     obj.imtags{k}=zeros(1,obj.metadata.numberOfFrames);
%                 end
                obj.init=false;
%          end  
    end   
    
    for k=1:length(obj.readoutimgtags)
        try
            tag=imgmeta.get(obj.readoutimgtags{k});
        catch err
            tag=NaN;
            err
        end
        if ischar(tag)
            obj.imtags(k,imagenumber)=str2double(tag);
        else
            obj.imtags(k,imagenumber)=(tag);
        end
    end
end

% obj.imtags{:,imagenumber}
% else
%     image=[];
% end

%    if imagenumber<=obj.reader.getImageCount()
%        image=bfGetPlane(obj.reader,imagenumber);
%    else
%        image=[];
%    end
end

function initimagetags(obj)
        camset=obj.getPar('loc_cameraSettings');
%         obj.imtags={};
        if ~isempty(camset)&& myisfield(camset,'imagemetadata') && ~isempty(camset.imagemetadata)
            if iscell(camset.imagemetadata)
                obj.readoutimgtags=camset.imagemetadata;
            else
                obj.readoutimgtags=split(camset.imagemetadata,',');
            end
            obj.readoutimgtags=strtrim(obj.readoutimgtags);
            obj.init=true;
        else 
            obj.readoutimgtags={};
        end    
end
