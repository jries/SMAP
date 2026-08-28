function [ax,tab]=initaxis(handle,name,option)
if nargin<3
    option='';
end

if isa(handle,'matlab.ui.container.Tab')
    
    tab=handle;
    ax=tab.Children;
    tabgroup=tab.Parent;
%     state=tab.Parent.Parent.Visible;
    if ~strcmp(option,'keep')
    delete(handle.Children(:));
    handle.Title=name;
    ax=axes('Parent',handle);
    
    end
    
    
elseif isa(handle,'matlab.ui.container.TabGroup')
%     state=tab.Parent.Visible;
    children=handle.Children;
    found=0;
    for k=1:length(children)
        if strcmpi(name,children(k).Title)
            found=1;
            tab=children(k);
            if ~strcmp(option,'keep') %|| ~isa(children(k))
            delete(children(k).Children(:));
            
            ax=axes('Parent',children(k));
            else
                iax=find(arrayfun(@(o) isa(o, 'matlab.graphics.axis.Axes'), tab.Children));
                ax=tab.Children(iax(1));
                ileg=arrayfun(@(o) isa(o, 'matlab.graphics.illustration.Legend'), tab.Children);
                delete(tab.Children(ileg));
            end
            break
        end
    end
    if found==0
        tab=uitab(handle,'Title',name);
        ax=axes('Parent',tab);
    end
    tabgroup=handle;

elseif isa(handle,'Axes')
%     
    ax=handle;
    delete(ax.Children(:));
    tab=ax.Parent;
    tabgroup=tab.Parent;
elseif isa(handle,'matlab.ui.Figure')
    ax=axes(handle);
    return
end
tabgroup.SelectedTab=tab;
state=tab.Parent.Parent.Visible;
 % axes(ax); %%%%XXXXX taken out
 fig=tab.Parent.Parent;
 fig.Visible=state;
