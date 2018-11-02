function il=getimageloader(obj,filename)
filename=strrep(filename,'\','/');
[~,~,ext]=fileparts(filename);

if isempty(ext)&&~strcmp(filename(end),'/') %directory
    filename=[filename '/'];
end
try
    il=imageloaderAll(filename,[],obj.P); %still exist?
    il.getimage(1);
catch
    maindir=obj.getGlobalSetting('DataDirectory');
    filenamef=findfilepath(filename,maindir);
    try
        il=imageloaderAll(filenamef,[],obj.P); %still exist?
        il.getimage(1);
    catch
        lastsml=obj.getPar('lastSMLFile');
        if ~isempty(lastsml)
        filenamef=findfilepathcomp(filename,fileparts(lastsml));
            try
                il=imageloaderAll(filenamef,[],obj.P); %still exist?
                il.getimage(1);
            catch
                [~,~,ext]=fileparts(filename);
                if isempty(ext)
                    filenamef=[filename '.tif'];
                end
                [f,path]=uigetfile(filenamef);
                if f
                    il=imageloaderAll([path f],[],obj.P); 
                else
                    il=[];
                end
            end
        end
    end
    %look for file in main directory
end
end

function filename=findfilepathcomp(filename,maindir)
% test: find common directory names
 filename=strrep(filename,'\','/');
 maindir=strrep(maindir,'\','/');
 [substr,loc1,loc2]=commonsubstring(filename,maindir);

filename=[maindir(1:loc2{1}) filename(loc1{1}:end)];
end

function filename=findfilepath(filename,maindir)
    d=dir(maindir);
    alldir={d([d.isdir]).name};
    ind=[1 strfind(filename,'/')];
    for k=1:length(ind)-1
        thisf=filename(ind(k)+1:ind(k+1)-1);
        if any(strcmp(alldir,thisf))&&~isempty(thisf)


            filename=[maindir '/' filename(ind(k)+1:end)];
            break
        end
    end

end