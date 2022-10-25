function [pid, status, results]=systemcallpython(envpath,command,runpath,logfile)
if nargin<3 || isempty(runpath)
    runpath=pwd;
end
cdir=pwd;
if ispc
    [p1,env]=(fileparts(envpath));
    condapath=fileparts(p1);
    pcall=['call "' condapath '\Scripts\activate.bat" ' env ' & cd "' runpath '" & ' command ];
else
    cd(runpath)
     pcall=[envpath '/bin/'  command ];
end
if nargin>3
%     if ispc
        pcall=[pcall ' 1> "' logfile '" 2>&1'];
%     else
%         pcall=[pcall ' > ' logfile];
%     end
end
pcall=strrep(strrep(pcall,'/',filesep),'\',filesep);
cstr=strsplit(command);
if ispc
    taskname=[cstr{1} '.exe'];
else
    taskname=[ cstr{1}];
end
[pidold,tasknames]=getpids(taskname);
[status, results]=system([pcall ' &'],'-echo');
[pidnew,tasknew]=getpids(taskname);
pid=setdiff(pidnew,pidold);
cd(cdir);
end

