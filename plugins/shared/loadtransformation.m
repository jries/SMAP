function transformation=loadtransformation(obj,file,dataset)
if isa(obj,'interfaces.LocalizationData')
    locData=obj;
else
    locData=obj.locData;
end
if nargin<3
    dataset=length(locData.files.file);
end
if exist(file,'file')
    l=load(file);
    if isfield(l,'transformation')
         transformation=l.transformation;
    else
        transformation=l.parameters_g.transformation;
    end
elseif isfield(locData.files.file(1),'transformation')
     if ~isempty(locData.files.file(dataset).transformation)
            transformation=locData.files.file(dataset).transformation;
     else
        for k=length(locData.files.file):-1:1
            if ~isempty(locData.files.file(k).transformation)
                transformation=locData.files.file(k).transformation;
                break
            end
        end
     end  
else
    transformation=[];
end