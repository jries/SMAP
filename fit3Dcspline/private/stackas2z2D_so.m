function zsx=stackas2z2D_so(sxa,za,na,p)
wind=4;
windn=10;

wind=wind/p.dz*50;
windn=windn/p.dz*50;
midp=mean(p.fminmax);
range=find(abs(za-midp)<500/p.dz); %only consider vicinity of glass
sx=sxa(range);
% sy=sya(range);
z=za(range);
n=na(range);


[~, indm]=max(n./max(0.5,sx).^2);

winstd=round(500./p.dz);
range=max(1,indm-winstd):min(indm+winstd,length(n));

    [~,nindx]=min((sx(range)));


    nxind=round(nindx)+range(1)-1;
    

    fstart=max(nxind-wind,1);
    fstop=min(nxind+wind,length(n));
    

    sxfit=sx(fstart:fstop);
    x=z(fstart:fstop);

    
    warning('off')
    [psx,S] = polyfit(x,sxfit,2);


    
    fmin=-psx(2)/2/psx(1);
    zsx=fmin;



   
