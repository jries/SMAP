function zas=stackas2z_so(sxa,sya,za,na,p)
wind=4;
windn=10;

wind=round(wind/p.dz*50);
windn=round(windn/p.dz*50);
midp=mean(p.fminmax);
numf=round(max(500/p.dz,length(sxa)/8));
range=find(abs(za-midp)<numf); %only consider vicinity of glass
sx=sxa(range);
sy=sya(range);
z=za(range);
n=na(range);


[~, indm]=max(n./max(0.5,sx)./max(0.5,sy));

winstd=round(500./p.dz);
range=max(1,indm-winstd):min(indm+winstd,length(n));

    [~,nindx]=min((sx(range)));
    [~,nindy]=min((sy(range)));

    nxind=round((nindx+nindy)/2)+range(1)-1;
    
    fstartn=max(nxind-windn,1);fstopn=min(nxind+windn,length(n));

    indx1=find(sx(fstartn:fstopn)>sy(fstartn:fstopn),1,'first')+fstartn-1;
    indx2=find(sy(fstartn:fstopn)>sx(fstartn:fstopn),1,'first')+fstartn-1;
    sxind=max(indx1,indx2)-1;
    fstart=max(sxind-wind,1);
    fstop=min(sxind+wind,length(n));
    

    sxfit=sx(fstart:fstop);syfit=sy(fstart:fstop);
    x=z(fstart:fstop);

    
    warning('off')
    [psx,S] = polyfit(x,sxfit,2);
    [psy,S] = polyfit(x,syfit,2);

    
    fmin=-psx(2)/2/psx(1);
    zsx=fmin;
    fminy=-psy(2)/2/psy(1);
    zsy=fminy;

%     S.normr
    nfit=n(fstart:fstop);
    

    p = polyfit(x,nfit,2);
%    warning('on')
    fmin=-p(2)/2/p(1);
    zn=fmin;
     warning('on')
    ro=roots(psy-psx);

    
    if S.df<3||S.normr>150||length(x)<5||length(ro)<2
        zas=NaN;
    else
        mp=z(fstop);
        d=abs(ro-mp);
        [~,mi]=min(d);
        zas=ro(mi);
        if ~(imag(zas)==0)
            zas=NaN;
        end

    end

