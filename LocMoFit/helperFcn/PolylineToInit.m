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
elseif size(out.xnmR,1)<Nctrlpoints  && Nctrlpoints~=2*size(out.xnmR,1)-1
    warning('Number of control points is lagger than the size of the polyline. %i equidictant points will be selected along the specified polyline by using interparc.',Nctrlpoints)
    temp=interparc(Nctrlpoints,out.xnmR,out.ynmR,out.znmR,'spline');
    out.xnmR=temp(:,1);
    out.ynmR=temp(:,2);
    out.znmR=temp(:,3);
elseif Nctrlpoints==2*size(out.xnmR,1)-1
    warning('Taking midpoints')
    oldpoints=1:2:2*size(out.xnmR,1)-1;
    %oldind=1:length(oldpoints);
    temp.xnmR=zeros(2*size(out.xnmR,1)-1,1);
    temp.ynmR=zeros(2*size(out.xnmR,1)-1,1);
    temp.znmR=zeros(2*size(out.xnmR,1)-1,1);
    c=1;
    for ptn=oldpoints
        if ptn==oldpoints(end)
            temp.xnmR(ptn)=out.xnmR(c);
            temp.ynmR(ptn)=out.ynmR(c);
            temp.znmR(ptn)=out.znmR(c);
        else
            temp.xnmR(ptn)=out.xnmR(c);
            temp.ynmR(ptn)=out.ynmR(c);
            temp.znmR(ptn)=out.znmR(c);
            temp.xnmR(ptn+1)=(out.xnmR(c)+out.xnmR(c+1))/2;
            temp.ynmR(ptn+1)=(out.ynmR(c)+out.ynmR(c+1))/2;
            temp.znmR(ptn+1)=(out.znmR(c)+out.znmR(c+1))/2;
        end
        c=c+1;        
    end
    out=temp;

end

%fieldvalues = @(pars)(cellfun(@(fieldName)(pars.(fieldName)),fieldnames(pars)));
%%Need to add error in case there is no polyline
end
