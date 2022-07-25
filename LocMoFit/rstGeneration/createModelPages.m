%% Settings
global fid

path2Rst = 'C:\Users\ries\git\SMAP\LocMoFit\docs\modLibrary.rst';
dir_img = 'C:\Users\ries\git\SMAP\LocMoFit\docs\images\models\';
fid = fopen(path2Rst,'w','a');

%% Get all the listed models
model2Show = modelList;

h1_rst('Model library');
rst_explicit_markup('automodule','models');
linebreak_rst;

for k = 1:length(model2Show)
    oneModel = model2Show{k};
    rst_explicit_markup('autoclass',oneModel, 'members','','show-inheritance','');
    path2img = [dir_img oneModel '.png'];
    if exist(path2img)
        rst_explicit_markup('image',['./images/' oneModel '.PNG'], 'width', '500')
    end
    linebreak_rst;
end
fclose(fid);