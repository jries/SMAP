function [md,allmd,alls]=getMetadataTifMM(ff, fields, fieldnumeric)
reader = javaObjectEDT('org.micromanager.acquisition.TaggedImageStorageMultipageTiff',fileparts(ff), false, [], false, false, true);

imgmetadata=reader.getImageTags(0,0,0,0);
summarymetadata=reader.getSummaryMetadata;

allmd=gethashtable(imgmetadata);
alls=gethashtable(summarymetadata);

if isempty(reader.getImageTags(0,1,0,0))
    impos=[0,0,1,0];
else
    impos=[0,1,0,0];
end
frames=summarymetadata.get('Frames');
fr=1;
imposc=num2cell(fr*impos);
imgmeta=reader.getImageTags(imposc{:});
if nargin<2 %ask for possible fields
    md=gethashtable(imgmeta);
    return
end

for fi=1:length(fields)
    val=imgmeta.get(fields{fi});
    if isnumeric(val)
        fieldnumerici(fi)=true;
        md{fi}=zeros(frames,1);
    else
        fieldnumerici(fi)=false;
        if  (nargin>=3 && fieldnumeric(fi))
            md{fi}=zeros(frames,1);
        else
            md{fi}{frames}=[];
        end
    end
end

if nargin<3
    fieldnumeric=fieldnumerici;
    convertval=false(length(fields),1);
else
    convertval=fieldnumeric~=fieldnumerici;

end

notcomplete=false;
for fr=1:frames
    imposc=num2cell(fr*impos);
    imgmeta=reader.getImageTags(imposc{:});
    if isempty(imgmeta)
        notcomplete=true;
        break
    end
    for fi=1:length(fields)

        val=imgmeta.get(fields{fi});
        if convertval(fi)
            val=str2double(val);
        end
        if fieldnumeric(fi)
            md{fi}(fr)=val;
        else
            md{fi}{fr}=val;
        end
    end
end
if notcomplete
    for fi=1:length(fields)
        md{fi}=md{fi}(1:fr-1);
    end
end

end




% 
% [md,allmd,alls]=getMetadataTifMM(ff);
% [fields,sortind]=sort(md(:,1));
% ld=listdlg('ListString',fields,'SelectionMode','multiple');
% usefields=fields(ld);
% md(sortind(ld),1)
% 
% isnumericval=false(1,length(ld));
% for k=1:length(ld)
%     val=md(sortind(ld(k)),2);
%     valn=str2double(val);
%     if ~isnan(valn)
%         isnumericval(k)=true;
%     end
% end
% 
% mdframe=getMetadataTifMM(ff,usefields,isnumericval);