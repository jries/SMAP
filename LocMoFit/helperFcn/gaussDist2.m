function likelihood = gaussDist2(refLocs,y,x,ysigmao,xsigmao,varargin)
    p = inputParser;
    addParameter(p, 'sigFactor',1);
    parse(p, varargin{:})
    % For one dimension: exp(-((x-mu/sigma)^2)/2)/(2*pi)^(1/2)*sigma
    f = p.Results.sigFactor;          % scale factor of the 3 sigma
    xsigma = xsigmao.*f;ysigma = ysigmao.*f;
    
    % the template bellow is for multiple sets of parameters at the same
    % time
    templateAll = zeros([1 size(refLocs.x,2) size(refLocs.x,1)]);
    template.x = templateAll;
    template.y = templateAll;
    template.n = templateAll;
    template.x(1,:,:)=refLocs.x';
    template.y(1,:,:)=refLocs.y';
    template.n(1,:,:)=refLocs.n'.*ones(size(refLocs.x,2),1);
    refLocs = template;
    termX = ((x - refLocs.x)./xsigma).^2;
    termY = ((y - refLocs.y)./ysigma).^2;
    svol = template.n.*(2*pi^-1./(xsigma.*ysigma));
    pairwiseVal = zeros(size(svol));
    lZero = svol(:) == 0;
    svol = svol(~lZero);
    termX = termX(~lZero);
    termY = termY(~lZero);
%     indKept = find(~lZero);
    pairwiseVal(~lZero) = sum(svol.*exp(-(termX+termY)/2),3);
    likelihood = sum(pairwiseVal,3);
end