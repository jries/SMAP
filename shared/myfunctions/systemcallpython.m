function [pid, status, results]=systemcallpython(envpath,command,runpath,logfile)
if nargin<3 || isempty(runpath)
    runpath=pwd;
end
if ispc
    cdir=pwd;
    [p1,env]=(fileparts(envpath));
    condapath=fileparts(p1);
    pcall=['call "' condapath '\Scripts\activate.bat" ' env ' & cd "' runpath '" & ' command ];
else
    cd (runpath)
     pcall=[envpath '/bin/ '  command ' &'];
end
if nargin>3
    pcall=[pcall ' 1> "' logfile '" 2>&1'];
end
pcall=strrep(strrep(pcall,'/',filesep),'\',filesep);
cstr=strsplit(command);
taskname=[cstr{1} '.exe'];
pidold=getpids(taskname);
[status, results]=system([pcall ' &'],'-echo');
pidnew=getpids(taskname);
pid=setdiff(pidnew,pidold);
cd(cdir);
end

