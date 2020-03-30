function worked=savewebfile(fout,url)
worked=false;
try
    docbin=webread(url);
catch err
    err
        disp(['could not download file ' url '. Please ensure you have internet connection'])
    return
end

try
fid=fopen(fout,'w');
fwrite(fid,docbin);
fclose(fid);
catch err
    err
    fid
    disp(['could not write file ' fout ' locally.']);
end

worked=true;
end


