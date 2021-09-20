function [x,y,z,n]=dualRings(par,dx)
    if ~exist('par','var')
        x.type = 'continuous';
        x.par = getParameters();
        return
    end
    r = par.radius;
    d = par.ringDistance;

    numphi=round(2*pi*r/dx);
    theta = linspace(0,2*pi,numphi);
    theta(end)=[];

    x = r.*cos(theta); y = r.*sin(theta);
    z = repelem(d/2, length(theta));

    x = repmat(x,1,2); y = repmat(y,1,2);
    z = [z -z];
    n = ones(size(x));
end

function par = getParameters()
    par = [];
    par.radius = [];
    par.ringDistance= [];
end