function struct2uitable(htable, s,flip, formats)
if nargin<4
    cstring=false;
else
    cstring=true;
end
if nargin<3||isempty(flip)
    fliph=false;
else
    if strcmp(flip,'flip')
        fliph=true;
    else 
        fliph=false;
    end
end
if cstring
    for k=1:length(s)
        sh=s(k);
        fn=fieldnames(sh);
        for l=1:length(fn)
            if isnumeric(sh.(fn{l}))
                so(k).(fn{l})=num2str(sh.(fn{l}),formats);
            else
                so(k).(fn{l})=sh.(fn{l});
            end
        end
    end
  s=so;  
end

tab=struct2table(s);
% tab.names=tab.Properties.VariableNames;

if fliph
    data=table2cell(tab)';
    data=[tab.Properties.VariableNames' data];
%     horzcat(tab.Properties.VariableNames, data)
    htable.Data=data;
    
    htable.RowName={};
else
    htable.Data=table2cell(tab);
    htable.ColumnName=tab.Properties.VariableNames;
end


