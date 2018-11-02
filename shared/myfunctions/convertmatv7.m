% convert -v7.3 files to v7 files
global path
[f,path]=uigetfile([path '*.mat'],'MultiSelect','on');
if ~iscell(f)
    f={f};
end
for k=1:length(f)
    
    fid=fopen([path f{k}]);
    x=fscanf(fid,'%c',22);
    
    fclose(fid);
    
    disp([x(1:21)  ': ' path f{k}] )
    if contains(x,'7.3')
        l=load([path f{k}]);
        v=saverightversion([path 'v7_' f{k}],l);
        disp([v ': ' path 'v7_' f{k}]);
    end
end

disp('done')