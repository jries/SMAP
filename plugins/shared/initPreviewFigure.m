function initPreviewFigure(obj)
f=obj.getPar('loc_outputfig');
if isempty(f) || ~isvalid(f)
    f=figure;
    obj.setPar('loc_outputfig',f);
end
figure(f);