%% ask for user input
[file,path2setting] = uigetfile({'*_avg.mat','resulted averages'}, 'Select an avg.mat file...', '');

%% load files
load([path2setting file])

%% load the results into SMAP
for k = 1:length(finalAvg)
    g.locData.addfile(['Particle fusion iteration ' num2str(k-1)]);
    finalAvg{k}.frame = ones(size(finalAvg{k}.xnm));
    oneAvg = finalAvg{k};
    g.locData.addLocData(oneAvg)
end
initGuiAfterLoad(g)