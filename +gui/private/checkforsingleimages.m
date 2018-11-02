function si=checkforsingleimages(file)
si=0;
il=imageLoader(file).info;
if (strcmp(il.format,'separateTif')||strcmp(il.format,'stackTif'))&&il.numberOfFrames>10
    button=questdlg('Raw Tif file detected. Localize these images?','Tif loader','Localize Images','Use for rendering','Cancel','Localize Images');
    if strcmpi(button,'Localize Images')
        si=1;
    elseif strcmpi(button,'Use for reconstruction')
        si=2;
    end
    
end
end