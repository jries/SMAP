function htable = createConvertTable(oldh, obj)
    % TODO: this function is under an integration to the SMLMModelFit
    % object. Please use the function there instead of this one.
    pos=oldh.Position;
    htable=uitable(oldh.Parent,'Data',{},'Position',[pos(1:2)+[5 10] 300 200]);
    
    if isa(obj,'ROIManager.Evaluate.LocMoFitGUI')
        colNames={'Source', 'Rule', 'Target_fit', 'Target_usr'};
        % check the loaded modules and hook all the LocMoFitGUI
        if length(obj.locData.SE.processors.eval.guihandles.modules.Data)>1
            loadedModuls = obj.locData.SE.processors.eval.guihandles.modules.Data(:,2);
            lLocMoFitGUI = contains(loadedModuls, 'LocMoFitGUI');
        else
            lLocMoFitGUI = 0;
        end
        if sum(lLocMoFitGUI)>0
            loadedLocMoFitGUI = loadedModuls(lLocMoFitGUI,:);
        else
            loadedLocMoFitGUI = [];
        end
        htable.CellEditCallback = {@convertTable_callback,4};
        colFormat = {[{'this step'}; loadedLocMoFitGUI]',[],{'none'},[]};
    else
        colNames={'Rule', 'Target_fit', 'Target_usr'};
        loadedLocMoFitGUI = [];
        htable.CellEditCallback = {@convertTable_callback,3};
        colFormat = {[],{'none'},[]};
    end
    htable.ColumnName = colNames;
    htable.ColumnFormat = colFormat;
    htable.ColumnEditable = true;
    htable.CellSelectionCallback = {@convertTableSelection_callback,obj};
    htable.RowName = [];
    delete(oldh);
end

function convertTable_callback(a,b,k) 
    if b.Indices(2) == k
        optionTarget = unique([a.ColumnFormat{3} {['usr_' b.NewData]}]);
        a.ColumnFormat{3} = optionTarget;
    end
end

function convertTableSelection_callback(a,b,obj) 
    try
    selectedRow = b.Indices(1);
        if isa(obj,'ROIManager.Evaluate.LocMoFitGUI')
            obj.setPar(['selectedRowConvert_' obj.name], selectedRow)
        elseif isa(obj,'ROIManager.Analyze.SMLMModelFit_dynamicRec_mCME')
            obj.setPar(['selectedRowConvert_' 'SMLMModelFit_dynamicRec_mCME'], selectedRow)
        elseif isa(obj,'matlab.ui.Figure')
            temp = findobj(obj,'Type','uicontrol','-and','String','selectedRowConvert');
            if isempty(temp)
                uicontrol(obj,'Style','text','String','selectedRowConvert','Value',selectedRow,'Visible','off')
            else
               temp.Value = selectedRow;
            end
        end
    catch
    end
end
