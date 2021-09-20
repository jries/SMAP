function rmRule_callback(a,b,obj, fn)
% Delete a row in the table
    htable = obj.guihandles.(fn);
    data = htable.Data;
    
    if isprop(obj,'name')
        name = obj.name;
    else
        name = obj.pluginpath{end};
    end
    
    selectedRow = obj.getPar(['selectedRowConvert_' name]);
    try
        % try here for unexpected row numbers caused by rows that has been
        % removed
        data(selectedRow,:) = [];
    catch 
    end
    htable.Data = data;
end