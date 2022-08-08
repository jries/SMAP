function helpLocMoFitGUI(a,b,obj,varargin)
% This function direct the click to the right web page.
%
%
    %% Find the current activated tab
    tabGroup_layer1 = obj.guihandles.tab1.Parent;
    tab_layer1_selected = tabGroup_layer1.SelectedTab;
    ID1 = tab_layer1_selected.Title;
    
    tabGroup_layer2 = findobj(tab_layer1_selected.Children,'-class','matlab.ui.container.TabGroup');
    if ~isempty(tabGroup_layer2)
        tab_layer2_selected = tabGroup_layer2.SelectedTab;
        ID2 = tab_layer2_selected.Title;
    else
        ID2 = [];
    end

    domain = 'LocMoFit/docs/_build/html/';
    domain = 'https://locmofit.readthedocs.io/en/latest/';
    
    if ~isempty(varargin)&&varargin{1}
        modelID_str = ID1(2:end);
        modelNames = obj.guihandles.(['modelname_CM_' modelID_str]);
        model_selected = modelNames.String{modelNames.Value};
        url = [domain, 'LocMoFit.modelLibrary.html#models.' model_selected];
    else
        if startsWith(ID1,'M')
            ID1 = 'M';
        end
        switch ID1
            case 'Settings'
               url = [domain, 'basics/GUIOverview.html#tab-settings'];
            case 'M'
               switch ID2
                   case 'Model'
                       rUrl = 'GUIOverview.html#sub-tab-model';
                   case 'Parameters'
                       rUrl = 'GUIOverview.html#sub-tab-parameters';
                   case 'Advance'
                       rUrl = 'GUIOverview.html#sub-tab-advance';
               end
               url = [domain 'basics/' rUrl];
            case 'Convert'
                url = [domain, 'basics/GUIOverview.html#tab-convert'];
            otherwise
        end
    end


    web(url)
end