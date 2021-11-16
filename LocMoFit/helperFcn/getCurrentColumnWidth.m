function columnWidth = getCurrentColumnWidth(hTable)
    jScroll = findjobj(hTable);
    jTable = jScroll.getViewport.getView;
    columnWidth = {};
    for k = jTable.getColumnCount:-1:1
       columnWidth{k} = jTable.getColumnModel.getColumn(k-1).getWidth./1.5;
    end
end