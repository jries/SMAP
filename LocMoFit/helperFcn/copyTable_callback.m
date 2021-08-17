function copyTable_callback(a,b,t)
% Copy the contents of a table to the clipboard.
    header = t.Properties.VariableNames;
    mat2clip(cellstr(string([header;table2cell(t)])));
end