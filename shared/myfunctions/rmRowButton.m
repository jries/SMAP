function h = rmRowButton(fig,htable)
    htable.CellSelectionCallback = @cellSelect;
    h = uicontrol(fig, 'Style','pushbutton','String','-','Callback',{@rmRow,htable});
    
end
function rmRow(a,b,htable)
    % Delete a row in the table
    data = htable.Data;
    selectedRow = htable.UserData.selectedRow;
    
    try
        % try here for unexpected row numbers caused by rows that has been
        % removed
        data(selectedRow,:) = [];
    catch 
    end
    htable.Data = data;
end

function cellSelect(a,b)
    try
        b.Source.UserData.selectedRow = b.Indices(1);
    catch 
    end
end