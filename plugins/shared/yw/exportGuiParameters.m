function exportGuiParameters(obj, path)
    p = obj.getGuiParameters;
    if nargin==1
        [file,path] = uiputfile('*.mat','Save guiParameters to',[obj.pluginpath{end} '_guiPar.mat']);
        path = [path file];
    elseif nargin==2
    else
        warning('Invalid inputs.')
        return
    end
    save(path, 'p')
end