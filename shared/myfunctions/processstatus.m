function [running,out]=processstatus(pid)
running=false;
if ismac
    [response,tasks] = system(['ps -x']);
    ll=strip(splitlines(tasks));
    ind=find(startsWith(ll,num2str(pid)));
    if isempty(ind)
        out='';
    else
        out=ll{ind};
        running=true;
    end
end

end