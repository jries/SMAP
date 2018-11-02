function fileout=getbasefile(filein)
path=strrep(fileparts(filein),'\','/');
inds=strfind(path,'/');
if strcmp(path(inds(end)+1:end-1),'Pos')
    fileout=path(1:inds(end)-1);
else
    fileout=path;
end
end