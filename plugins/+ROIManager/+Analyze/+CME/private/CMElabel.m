function txt=CMElabel(obj,event,label,image,imaxis)
pos=event.Position;
data=obj.Parent.Children;
idx=[];
for k=1:length(data)
    dh=data(k);
    if isa(dh,'matlab.graphics.chart.primitive.Line')
        idx=find(dh.XData==pos(1)&dh.YData==pos(2),1,'first');
        if ~isempty(idx)
            break
        end
    end
end

txt{1}=['X: ' num2str(pos(1))];
txt{2}=['Y: ' num2str(pos(2))];
if ~isempty(idx)
txt{3}=label{idx};
if nargin>4
    imagesc(1.3*image{idx}/max(image{idx}(:)),'Parent',imaxis)
    title(txt{3},'Parent',imaxis)
end
end
    