function pid=getpids(name,varargin)
[response,tasks] = system(['tasklist/fi "imagename eq "' name]);
result = regexp(tasks, ['(?<=' strrep(name, '.', '\.') '\s*)\d+'], 'match');
pid = str2double(result); % convert to numbers if needed
if nargin>1 && contains(varargin{1},'kill')
    if ispc
        cmd = ['Taskkill /IM ' name ' /F'];
    else
        cmd = ['kill' name];
    end
    system(cmd);
end
end