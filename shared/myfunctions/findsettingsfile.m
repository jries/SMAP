function fo=findsettingsfile(fin,obj)
if exist(fin,'file')
    fo=fin;
else
    settingsdir=obj.getPar('SettingsDirectory');
    str=strrep(fin,'\','/');
    ind=strfind(str,'/settings/');
    if isempty(ind)
        fo=strrep(fin,'settings',settingsdir);
        return;
    end
   
    fo=strrep(fin(ind+1:end),'settings',settingsdir);
    fo=strrep(strrep(fo,'/',filesep),'\',filesep);
end