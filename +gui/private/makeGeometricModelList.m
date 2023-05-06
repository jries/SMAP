function makeGeometricModelList
%     this function writes the m file 'modelList.m'
if ~isdeployed
    allModels = LocMoFit.getModelList();
    path_LocMoFit = which('LocMoFit');
    folder_LocMoFit = fileparts(path_LocMoFit);
    folder_LocMoFit = replace(folder_LocMoFit, '@LocMoFit', '');

    lListed = [];
    for k = length(allModels):-1:1
        tempModel = eval(allModels{k});
        lListed(k) = tempModel.listed;
    end

    allModels = allModels(logical(lListed));

    head = ['function modelObj = modelList(modelName)' newline...
        'if nargin>0' newline...
        'switch modelName' newline];
    body = [];
    bottom = ['end' newline...
        'else' newline...
        'modelObj = {' char(39) strjoin(allModels, [char(39) ',' char(39)]) char(39) '};' newline...
        'end' newline...
        'end'];
    for k = 1:length(allModels)
        oneModel = ['case ' char(39) allModels{k} char(39) newline...
            'modelObj = ' allModels{k} ';' newline];
        body = [body oneModel];
    end
    output = [head body bottom];
    writelines(output, [folder_LocMoFit 'modelList.m'])
end
end

function writelines(lines, filename)
    fid = fopen(filename,'wt');
    fprintf(fid, lines);
    fclose(fid);
end