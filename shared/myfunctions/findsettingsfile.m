function fo=findsettingsfile(fin)
if exist(fin,'file')
    fo=fin;
else
    str=strrep(fin,'\','/');
    ind=strfind(str,'/settings/');
    if isempty(ind)
        fo='';
        return;
    end
    fo=[pwd fin(ind:end)];
    fo=strrep(strrep(fo,'/',filesep),'\',filesep);
end