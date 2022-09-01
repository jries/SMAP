function out=PolylineToInit(input,pos,angle,Nctrlpoints)
%% Transforms 3D polyline to the XYZ coordinates of the initial model
%only for 3D polyline and 3D data
pars=input;
pars(:,1)=input(:,1)-pos(1,1);
pars(:,2)=input(:,2)-pos(1,2);
out=[];
[out.xnmR, out.ynmR]=rotcoord(pars(:,1),pars(:,2),angle);
out.znmR=pars(:,3);

if size(out.xnmR,1)>Nctrlpoints
    warning('Number of control points is smaller than the size of the polyline. Only first %i points will be used.',Nctrlpoints)
    temp.xnmR=out.xnmR(1:Nctrlpoints,:);
    temp.ynmR=out.ynmR(1:Nctrlpoints,:);
    temp.znmR=out.znmR(1:Nctrlpoints,:);
    out=temp;
elseif size(out.xnmR,1)<Nctrlpoints
    warning('Number of control points is lagger than the size of the polyline. %i equidictant points will be selected along the specified polyline by using interparc.',Nctrlpoints)
    temp=interparc(Nctrlpoints,out.xnmR,out.ynmR,out.znmR,'spline');
    out.xnmR=temp(:,1);
    out.ynmR=temp(:,2);
    out.znmR=temp(:,3);
end

%fieldvalues = @(pars)(cellfun(@(fieldName)(pars.(fieldName)),fieldnames(pars)));
%%Need to add error in case there is no polyline
end
