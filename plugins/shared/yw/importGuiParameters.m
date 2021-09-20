function importGuiParameters(obj, path)
    if nargin==1
        [file,path] = uigetfile('*.mat');
        path = [path file];
    elseif nargin==2
    else
        warning('Invalid inputs.')
        return
    end
    load(path);
    obj.setGuiParameters(p)
end