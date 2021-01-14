function newf = flattenTabs(varargin)
    if isa(varargin{1},'matlab.ui.Figure')
        f = varargin{1};
        varargin(1) = [];
    else
        f = gcf;
    end
    nrow = varargin{1};
    
    tabGrp = findobj(f,'-class','matlab.ui.container.TabGroup');
    nTabs = size(tabGrp.Children,1);
    ncol = ceil(nTabs/nrow);
    
    newf = figure;
    
    for k = nTabs:-1:1
        [currentCol, currentRow] = ind2sub([nrow ncol],k);
        newFtabGrp{k} = uitabgroup(newf);
        c = copy(tabGrp.Children(k));
        c.Parent = newFtabGrp{k};
        newFtabGrp{k}.Position = [(currentCol-1)/ncol (nrow-currentRow)/nrow 1/ncol 1/nrow];
    end
end