function fo=findunusedfilename(fin)

% if ~exist(fin,'file')
%     fo=fin;
%     return;
% end

if ~isempty(strfind(fin,'_sml'))
    searchstr='_sml';
else
    [~,~,searchstr]=fileparts(fin);
end

ind=1;
fo=fin;
while exist(fo,'file')
    inds=strfind(fo,searchstr);
    indm=strfind(fo(1:inds-1),'_');
    insert=int2str(ind);
    if ~isempty(indm)&&indm(end)<inds(end)
        num=fo(indm(end)+1:inds(end)-1);
        numn=str2double(num);
        if ~isnan(numn)
            insert=int2str(numn+1);
            fo=fo([1:indm(end)-1 inds(end):end]);
        end
        ind=ind-1;
    end
        
    fo=strrep(fo,searchstr,['_' insert searchstr]);
    ind=ind+1;
end