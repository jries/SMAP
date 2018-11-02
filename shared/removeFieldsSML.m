% removeFieldsSML

fieldstosave={'phot','frame','PSFxnm','locprecnm','xnm','ynm','bg','channel','filenumber','dummy'};

startdirectory='C:\Users\ries\Documents\Mund_RingPaper_Data\Radial\GFP\*_sml.mat';
h=selectManyFiles(startdirectory);
waitfor(h.handle)

filelist=h.filelist;
for k=1:length(filelist)
    disp(['load file: ' filelist{k}]);
    l=load(filelist{k});
    lout=l;
    fn=fieldnames(l.saveloc.loc);
    fieldsremove=setdiff(fn,fieldstosave);
    notfound=setdiff(fieldstosave,fn);
    if ~isempty(notfound)
        disp(['fields not found: ' notfound{1}]);
    end
    lout.saveloc.loc=rmfield(l.saveloc.loc,fieldsremove);
    newname=strrep(filelist{k},'_sml.mat','_s_sml.mat');
    if exist(newname,'file')
        newname=strrep(filelist{k},'.mat','_s_sml.mat');
    end
     disp(['save file: ' newname]);    
    v=saverightversion(newname,lout,'-v7');
    disp(v)
end