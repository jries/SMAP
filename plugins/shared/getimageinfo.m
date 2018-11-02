function info=getimageinfo(file)
[p,f,ext]=fileparts(file);
if isempty(f)
    ext='.tif';
    allf=dir([p filesep '*' ext]);
    f=allf(1).name;
    file=[p filesep f];
end
 info.metafile=[];
 info.filename=file;
switch ext
    case '.tif' 
        warning('off','MATLAB:imagesci:tiffmexutils:libtiffWarning');
        omeind=strfind(file,'.ome.');
        isstack=false;
        if ~isempty(omeind)
            isstack=true;
            numberOfFrames=0;
             basefile=file(1:omeind-1);
                         %look if there is first file
            indu=strfind(basefile,'_');
            if ~isempty(indu) && ~isnan(str2double(basefile(indu(end)+1:end)))
                basefile2=basefile(1:indu(end)-1);
                if exist([basefile2 '.ome.tif'],'file')
                    basefile=basefile2;
                end
            end
            info.filename=[basefile '.ome.tif'];
                        %if several files: determine size of all of them
            allfs=dir([basefile '_*.ome.tif']);
            basepath=fileparts(basefile);
            for k=1:length(allfs)
                allfs(k).name=[basepath filesep allfs(k).name];
            end
            
            firstfile= dir(info.filename);
            firstfile.name=[basepath filesep firstfile.name];
            allfiles=vertcat(firstfile,allfs);
        else
            infoim=imfinfo(file,'tif');
            if length(infoim)>1
                isstack=true;
                numberOfFrames=length(infoim);
                basefile=file;
                info.filename=file;
            
            allfiles=dir(file);
            allfiles.name=[fileparts(file) filesep allfiles.name];
            allfiles(1).numberOfFrames=numberOfFrames;
            end
        end
        if  isstack
            info.basefile=basefile;
            
            

            metaname=[info.basefile '_metadata.txt'];
            
            if exist(metaname,'file')
                try
                info.metafile=metaname;
                fid=fopen(metaname);
                nread=100000;
                fseek(fid,-nread,'eof');
                minfo=fread(fid,[1,nread],'*char');
                fclose(fid);
                ind=strfindfast(minfo,'"FrameKey-',1,-1);
                
                ind2=strfindfast(minfo,'-',ind);
                str=minfo(ind:ind2-1);
                numberOfFrames=str2double(str)+1;
                 allfiles(1).numberOfFrames=numberOfFrames;
                 
                catch
                    
                end
            end
%                 metadata=minfoparsec(minfo,camcalib);
            if length(allfiles)>1||numberOfFrames==0
                numberOfFrames=0;
                for k=1:length(allfiles)              
                    th=Tiff(allfiles(k).name,'r');
%                     numframes=getImagesInTiff(allfiles(k).name);
                    desc=th.getTag('ImageDescription');
                    ind=strfindfast(desc,'SizeT=',1,1);
                    if ~isempty(ind)
                        sts=desc(ind+7:ind+28);
                        numframes=sscanf(sts,'%i');

                    else
                        ind=strfindfast(desc,'images=',1,1);
                        sts=desc(ind+7:ind+28);
                        numframes=sscanf(sts,'%i');
                    end
                    numberOfFrames=numberOfFrames+numframes;             
                    allfiles(k).numberOfFrames=numframes;
                end
            end
            info.allfiles=allfiles;
            if (numberOfFrames==0)
                info.numerOfFrames=inf;
            else
                info.numberOfFrames=numberOfFrames;
            end
            info.frames=info.numberOfFrames;
            th=Tiff(allfiles(1).name,'r');
            info.Width=th.getTag('ImageWidth');
            info.Height=th.getTag('ImageLength');
            

            info.format='stackTif';
            th.close;
%             info.tiffh=th;
            return
        
        end
        
        searchstr='img_0';
        fall=myfastdir(p, [searchstr '*.tif']);
        if isempty(fall)
%         if isempty(f)||~strcmp(f(1:length(searchstr)),searchstr) %single image tiffs 
            searchstr='';
        
        fall=myfastdir(p, [searchstr '*.tif']);
        end
        
        try
            info.files=string(fall);
        catch
            info.files=fall;
        end
        
         p=strrep(p,'\','/');
        info.path=p;
        inds=strfind(info.path,'/');
        if strcmp(p(inds(end)+1:end-1),'Pos')
            info.basefile=p(1:inds(end)-1);
        else
            info.basefile=p;
        end
        info.numberOfFrames=length(fall);
        if isempty(f)
            file=[p filesep fall{1}];
        end
        fi=imfinfo(file);
        info.Width=fi.Width;
        info.Height=fi.Height;

        metaname=[p filesep 'metadata.txt'];
        if exist(metaname,'file')
            info.metafile=metaname;
        end
        info.format='separateTif';
        try
        if length(fall)>2
        %pattern
            n1=fall{2};n2=fall{3};
            inds1=strfind(n1,'_');
            inds2=strfind(n2,'_');
            for k=1:length(inds1)-1
                if ~strcmp(n1(inds1(k):inds1(k+1)),n2(inds2(k):inds2(k+1)))
                    info.numberNameRange=[inds1(k)+1,inds1(k+1)-1];
                    break
                end
            end
        else
             inds1=strfind(fall{1},'_');
             info.numberNameRange=[inds1(end-1)+1,inds1(end)-1];
        end     
        catch
            info.numberNameRange=[];
        end
        info.frames=info.numberOfFrames;
        return
        
end

end