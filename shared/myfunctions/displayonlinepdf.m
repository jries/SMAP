function displayonlinepdf(fout,url,iszip)
try
    docbin=webread(url);
catch err
    err
    if exist(fout,'file')
        warndlg('could not load file, displaying local copy instead')
        myopenpdf(fout);
    else
        warndlg('could not load file. Please ensure you have internet connection')
    end
    return
end

try
    [pfad, file]=fileparts(fout);
    foutzip=[pfad filesep file  '.tar'];
    if iszip
        fout=foutzip;
    end
    
fid=fopen(fout,'w');
fwrite(fid,docbin);
fclose(fid);
catch err
    err
    warndlg('could not write file locally');
end

if iszip
    unzip(foutzip,pfad)
end
myopenpdf(fout);
end


