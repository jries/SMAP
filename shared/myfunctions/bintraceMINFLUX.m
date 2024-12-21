function [stdxo,photbo,sigminflux]=bintraceMINFLUX(xh, photh, L, axx)
nump=floor(log2(length(xh)));
stdx=zeros(nump,1);photb=stdx;
for k=1:nump
    stdx(k)=std(xh);
    photb(k)=mean(photh);
    xh=sumbintrace2(xh)/2;
    if length(xh)<10
        break
    end
    photh=sumbintrace2(photh);
end
photbo=photb(1:k);
stdxo=stdx(1:k);
if nargin>2 && ~isempty(L) %plot
    sigminflux=L./sqrt(8*photb(1:k));
end
if nargin>3 && ~isempty(axx) %plot
    hold(axx,'off')
    semilogx(axx,photb(1:k), stdx(1:k),'ko-')
    hold(axx,'on')
    sigsmlm=120./sqrt(photb(1:k));
    
   
    semilogx(axx,photb(1:k), sigsmlm,'m--')
    semilogx(axx,photb(1:k), sigminflux,'b-.')
    ylabel(axx,'std pos (nm)')
    xlabel(axx,'photons')
    legend(axx,'std','SMLM','MINFLUX')
end
end

function xb=sumbintrace2(x)
len=floor(length(x)/2)*2;
xb=x(1:2:len-1)+x(2:2:len);
end
