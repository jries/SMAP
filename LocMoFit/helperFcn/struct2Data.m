function sortedData = struct2Data(parsArg)
%% STRUCT2DATA Convert the parsArg into the format that the uitable can handle
% A helper function for the SMLMModelFitGUI
typeInd = 6;
nameInd = 1;
lParOrder = {'x','y','z','zrot','xrot','yrot','weight','xscale','yscale','zscale'};
if isfield(parsArg, 'label')
    fn = {'name','value','fix','lb','ub','type','min','max','label'};
else
    fn = {'layer','value','fix','lb','ub','min','max'};
end
for k = length(parsArg.(fn{5})):-1:1
    for l = length(fn):-1:1
        uitableData{k,l} = parsArg.(fn{l})(k,:);
    end
end
% sort the rows based on the rules defined here
if ismember(fieldnames(parsArg), 'layer')
    lMPar = ismember(uitableData(:,typeInd),'mPar');
    sortedData = uitableData(lMPar,:);
    uitableData(lMPar,:)=[];
    for k = 1:length(lParOrder)
        lOnelPar = ismember(strtrim(uitableData(:,nameInd)), lParOrder{k});
        oneRow = uitableData(lOnelPar,:);
        if ~isempty(oneRow)
            sortedData(end+1,:) = oneRow;
        end
    end
else
    sortedData = uitableData;
end
end