function [outn,outind]=checknames(names,start)
%start: index or list of fieldnames
f=figure('ToolBar','none','MenuBar','none');
f.Position(3)=300;
ht=uitable('Parent',f,'Units','normalized','Position',[0,0.1,1,.9]);
uicontrol('Style','pushbutton','String','OK','Units','normalized','Position',[0.6,0.025,.2,.05],'Callback',{@okcallback,1})
uicontrol('Style','pushbutton','String','Cancel','Units','normalized','Position',[0.2,0.025,.2,.05],'Callback',{@okcallback,0})

if iscell(start)||ischar(start)
    ind=ismember(names,start);
else
    if islogical(start)
        ind=start;
    else
        ind=false(length(names),1);
        ind(start)=true;
    end
end
data(:,2)=names;
for k=1:length(names)
    data{k,1}=logical(ind(k));
end
    
ht.Data=data;
ht.ColumnEditable=[true false];
ht.ColumnWidth={30 220};
ht.ColumnName={};
waitfor(f)
    function okcallback(a,b,button)

        if button
            data=ht.Data;
            outind=[data{:,1}];
            outn=names(outind);
        else
            outind=false(length(names),1);
            outn={};
        end
        close(f)
    end

end