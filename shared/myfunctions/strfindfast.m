function ind=strfindfast(bigs,finds,startind,direction)
%bigs: string, finds: searchstring, startind, direction=1,-1(from back)
if nargin<=2
ind1=1;
else
    ind1=startind;
end
if nargin<4
    direction=1;
    if isempty(startind)
        startind=length(bigs);
    end
end

if direction==-1
    toind=max(length(bigs)-1e7,1);
    bigs=bigs(end:-1:toind);
    finds=finds(end:-1:1);
end
    dind=100000;
    lenf=length(finds);
    lenbig=length(bigs);
    ind2=ind1+dind;
    ind=[];
    while  ind2<=lenbig+dind

        indf=strfind(bigs(ind1:min(ind2,lenbig)),finds);
        if ~isempty(indf)
            ind=ind1+indf(1)-1;
            break;
        end
        ind1=ind2-lenf;
        ind2=ind1+dind;
    end

if direction==-1
    ind=length(bigs)-ind+2+toind-1;
end