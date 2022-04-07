function [pid, tasknames]=getpids(name,varargin)
if ispc
[response,tasks] = system(['tasklist/fi "imagename eq "' name]);
result = regexp(tasks, ['(?<=' strrep(name, '.', '\.') '\s*)\d+'], 'match');
pid = str2double(result); % convert to numbers if needed
if nargin>1 && contains(varargin{1},'kill')
%     if ispc
        cmd = ['Taskkill /IM ' name ' /F'];
%     else
%         cmd = ['kill' name];
%     end
    system(cmd);
end
else %mac
    [response,tasks] = system(['ps -x']);
    inds=strfind(tasks, name);
    offset=200;
    pid=[];tasknames={};
    for k=1:length(inds)
        teststr=tasks(inds(k)-offset:inds(k));
        teststrall=tasks(inds(k)-offset:min(inds(k)+2*offset,length(tasks)));
        indh=strfind(teststr,10);
        teststr2=teststr(indh(end)+1:end);
        ind2=find(teststr2>32);
        teststr2=teststr2(ind2(1):end);
        ind2=find(teststr2==32);
        pidstr=teststr2(1:ind2(1)-1);
        indel=strfind(teststrall(indh(end)+1:end),10);
        tasknames{k}=teststrall(indh(end)+1:indh(end)+indel(1)-1);
        pid(k)=str2double(pidstr);
    end
end

end