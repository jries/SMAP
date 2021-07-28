function v=saverightversion(file,ls,version)
if nargin <3 || isempty(version)
    v='-v7';
else
    v=version;
end
try
    notwork=false;
save(file,'-struct','ls',v);
catch err
    disp(getReport(err, 'extended', 'hyperlinks', 'on'))
    notwork=true;
end
[msg,msgid]=lastwarn;
if notwork || strcmp(msgid,'MATLAB:save:sizeTooBigForMATFile')
    save(file,'-struct','ls','-v7.3');
    disp('File is now being saved as v7.3');
    v='-v7.3';
    lastwarn('cleared');
end


