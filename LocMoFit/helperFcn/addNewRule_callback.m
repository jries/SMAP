function addNewRule_callback(a,b,obj,fn)
    % Create a new row in the table
    htable = obj.guihandles.(fn);
    
    % Indentify the number of columns
    numOfCol = length(htable.ColumnName);
    
    % Update the table data
    htable.Data = [htable.Data; cell(1,numOfCol)];
end

