function out=PolylineToInit(input,pos,angle)
%% Transforms 3D polyline to the XYZ coordinates of the initial model
%only for 3D polyline and 3D data
pars=input;
pars(:,1)=input(:,1)-pos(1,1);
pars(:,2)=input(:,2)-pos(1,2);
out=[];
[out.xnmR, out.ynmR]=rotcoord(pars(:,1),pars(:,2),angle);
out.znmR=pars(:,3);
%fieldvalues = @(pars)(cellfun(@(fieldName)(pars.(fieldName)),fieldnames(pars)));
%%Need to add error in case there is no polyline
end
