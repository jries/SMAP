function modifySettingsFile(file,varargin)
names={};
values={};
for k=1:2:length(varargin)
    names{end+1}=varargin{k};
    values{end+1}=varargin{k+1}; %need to be strings
end

fid=fopen(file);
if fid<0
    warning('modifySettingsFile: could not find fle')
    return
end
ind=1;
line=fgetl(fid);
while ischar(line)
    pline=parseline(line,names,values);
    textout{ind}=pline;
    ind=ind+1;
    line=fgetl(fid);
end
fclose(fid);
file2=file;
% file2(end-5)='X';
fid=fopen(file2,'w');
if fid<0
    warning('modifySettingsFile: could not open file for writing')
    return
end
for k=1:length(textout)
    fprintf(fid,'%s\n',textout{k});

end
fclose(fid);
end

function linout=parseline(line,names,values)
inde=strfind(line,'=');
    linout=line;
if isempty(inde)
    return;
end
inde=inde(1);
for k=1:length(names)
    if any(strfind(line(1:inde),names{k}))
        linout(inde+1:end)=[];
        linout=([linout values{k}]);
    end
end
    
end