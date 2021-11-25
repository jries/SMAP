function columnWidth = getCurrentColumnWidth(hTable)
if ~isempty(hTable.Data)
    jScroll = findjobj(hTable);
end

if exist('jScroll','var')&&~isempty(jScroll)
    jTable = jScroll.getViewport.getView;
    columnWidth = {};
    for k = jTable.getColumnCount:-1:1
        columnWidth{k} = jTable.getColumnModel.getColumn(k-1).getWidth./1.5;
    end
else
    columnWidth = hTable.ColumnWidth;
end
end