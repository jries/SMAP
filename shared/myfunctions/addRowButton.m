function h = addRowButton(fig,htable)
    h = uicontrol(fig, 'Style','pushbutton','String','+','Callback',{@addRow,htable});
    
end

function addRow(a,b,htable)
    % Indentify the number of columns
    numOfCol = length(htable.ColumnName);
    
    % Update the table data
    htable.Data = [htable.Data; cell(1,numOfCol)];
end