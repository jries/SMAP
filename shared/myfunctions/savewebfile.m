function worked=savewebfile(fout,url)
worked=false;
try
    options = weboptions('Timeout', 1);
    docbin=webread(url,options);
catch err
    err
        disp(['could not download file ' url '. Please ensure you have internet connection'])
    return
end

try
fid=fopen(fout,'w');
fwrite(fid,docbin);
fclose(fid);
worked=true;
catch err
    err
    fid
    disp(['could not write file ' fout ' locally.']);
end


end


