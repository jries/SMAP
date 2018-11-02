function mdinf=getmetadataMMnew(file)
info=getimageinfo(file);
mdinf(1,1:2)={'numberOfFrames info',info.numberOfFrames};
mdinf(end+1,1:2)={'Width info',info.Width};
mdinf(end+1,1:2)={'Height info',info.Height};
mdinf(end+1,1:2)={'format info',info.format};
if isfield(info,'numberNameRange')
    mdinf(end+1,1:2)={'number name range info',num2str(info.numberNameRange)};
end
metafile=info.metafile;

        fid=fopen(metafile);
        if fid>0

            
            
            minfo=fread(fid,[1,100000],'*char');
            nread=100000;
            status=fseek(fid,-nread,'eof');
            if status==0
                minfo2=fread(fid,[1,nread],'*char');
            else
                minfo2=minfo;
            end
            fclose(fid);
            tt=textscan(minfo,'%s','delimiter',',');
            t=tt{1};
            md=cell(length(t),2);
            indg=false(length(t),1);
            for k=1:length(t)
                x=textscan(t{k},'%s','delimiter',':' );
                xh=x{1};
                if length(xh)>1
                    indg(k)=true;
                    md(k,1)=strrep(xh(1),'"','');
                    md(k,2)=strrep(xh(2),'"','');
                end
                
            end
            md=md(indg,:);
            [~,ic]=unique(md(:,1));
            mdo=md(ic,:);   
            %remove frame keys
            mdo(strncmp(mdo(:,1),'FrameKey',8),:)=[];
            
            %manual
             ind=strfindfast(minfo,'"ROI": [',1);
             if isempty(ind)
                 ind=strfindfast(minfo,'"ROI": "',1);
             end
             troi=textscan(minfo(ind+10:ind+100),'%d','delimiter',',');
             mdo(end+1,:)={'ROI direct',num2str(troi{:}')};

             
            ind=strfindfast(minfo2,'"FrameKey-',1,-1);

            ind2=strfindfast(minfo2,'"',ind);
            str=minfo2(ind:ind2-1);
            if isempty(str)
                numberOfFrames=1;
            else
            numberOfFrames=max(cell2mat(textscan(str,'%f','delimiter','-')));
            end
%             ind3=strfind(str,'-');
%             str2=str(ind3(end)+1:end);
%             numberOfFrames=str2double(str2)+1;
            mdo(end+1,:)={'numberOfFrames frame key',num2str(numberOfFrames+1)};
            
            dd=myfastdir(fileparts(file), 'img_*.tif');
            numberoffiles=length(dd);
            
            mdo(end+1,:)={'frames direct',num2str(max(numberoffiles,numberOfFrames+1))};
            mdinf=vertcat(mdinf,mdo);
        end
             %sort
             [~,inds]=sortrows(mdinf,1);
             mdinf=mdinf(inds,:);
end