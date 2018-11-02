function PI = getPI(g, varargin)
% getPI(g), g is the GuiMainSMAP
% This is a function that Yu-Le created for extracting and ploting PIZStage-Position from SMLM files loaded
if ~isempty(varargin)
    pattern = varargin{1};
else
    pattern = 'PIZStage\-Position';
end

nFiles = g.locData.SE.numberOfFiles;
PI = zeros(1, nFiles);
for j = 1: nFiles
    path = g.locData.files.file(j).name;
    breakMark =regexp(path,'\\.');
    fileName = path(breakMark(end):end);
    fileName = regexprep(fileName, '\_Localization\_','\_bfp\_');
    fileName = regexprep(fileName, '\_sml\.mat','\_MMStack\_Pos0\.ome\.tif');
    path = regexprep(path, '\_Localization\_','\_bfp\_');
    path = regexprep(path, '\_sml\.mat','');
    path = [path fileName];

    info = imfinfo(path);
    if length(info.UnknownTags)==3
        [~,~,infoTag] = info.UnknownTags.Value;
    else
        infoTag = info.UnknownTags.Value;
    end
    
    PIZStage = regexprep(infoTag, ['.+' pattern], pattern);
    PIZStage = regexprep(PIZStage, '\,.*', '');
    PIZStage = regexprep(PIZStage, '.*\:', '');
    PIZStage = replace(PIZStage, '"', '');
    PI(j) = str2num(PIZStage);
end
    figure; plot(1:nFiles, PI)
end