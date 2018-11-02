function [zas,zn]=stackas2z(sx,sy,z,n,p)
ploton=p.ploton;
wind=4;
windn=10;

wind=wind/p.dz*50;
windn=windn/p.dz*50;

%     [~,nxind]=max(n);
%     [~,nxind]=min((sx+sy)./n);

%more pre-processing
% ntest=n;
% sxmin=myquantilefast(sx,.1);symin=myquantilefast(sy,.1);
% ind=sx<sxmin&sy<symin;
% zn=mean(z(ind));
% winstd=round(std(z(ind)));
% indm=round(zn);
[~, indm]=max(n./sx./sy);
winstd=round(500./p.dz);
range=max(1,indm-winstd):min(indm+winstd,length(n));

    [~,nindx]=min((sx(range)));
    [~,nindy]=min((sy(range)));
%     [~,nindx]=min((sx(range))./n(range));
%     [~,nindy]=min((sy(range))./n(range));
    nxind=round((nindx+nindy)/2)+range(1)-1;
    
    fstartn=max(nxind-windn,1);fstopn=min(nxind+windn,length(n));
%     [sminx,sxindx]=min(sx(fstartn:fstopn));
%     [sminy,sxindy]=min(sy(fstartn:fstopn));
%     sxindx=sxindx+fstartn-1;
%     sxindy=sxindy+fstartn-1;
%     sxind=mean(sxindx,sxindy);
    
    indx1=find(sx(fstartn:fstopn)>sy(fstartn:fstopn),1,'first')+fstartn-1;
    indx2=find(sy(fstartn:fstopn)>sx(fstartn:fstopn),1,'first')+fstartn-1;
    sxind=max(indx1,indx2)-1;
    fstart=max(sxind-wind,1);
    fstop=min(sxind+wind,length(n));
    
%     fstart=max(min(sxindx-wind,sxindy-wind),1);
%     fstop=min(max(sxindx+wind,sxindy+wind),length(n));
    sxfit=sx(fstart:fstop);syfit=sy(fstart:fstop);
    x=z(fstart:fstop);
    if 1%sminx<200
    
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
%         mp=(fstart+fstop)/2;
        mp=z(fstop);
        d=abs(ro-mp);
        [~,mi]=min(d);
        zas=ro(mi);
        if ~(imag(zas)==0)
            zas=NaN;
        end
%         if fstart<=ro(1)<=fstop
%             zas(k)=ro(1);
%         else
%             zas(k)=ro(2);
%         end
%         if ro(2)<3 %outside frame
%             zas(k)=ro(1);
%         end
    end
    
    if ploton %&& ~isnan(zas)
%         ro
%         disp(zas)
%         disp(S.df)
%         disp(length(x))
%         disp(S.normr)
        try

        subplot(2,2,1)
        hold off
        plot(x,sxfit,'ro');
        hold on
        plot(x,syfit,'go');
        plot(x,polyval(psx,x),'r')
        plot(x,polyval(psy,x),'g')
        plot(zas,polyval(psy,zas),'ko')
        subplot(2,2,3)
        hold off
        plot(x,nfit,'o');
        hold on
        plot(x,polyval(p,x),'r')
        subplot(2,2,2)
        plot(z,sx,'o',z,sy,'o',x,polyval(psx,x),'k',x,polyval(psy,x),'k');
        subplot(2,2,4)
        plot(z,n,x,polyval(p,x));
%         S
%         plot(zn(k),1,'x')
        drawnow
        end
    end
    else
        zs=NaN;
        zn=NaN;
        zas=NaN;
    end


% figure(1)
% hold off
% plot(x,sfit,'o');
% hold on
% plot(x,polyval(p,x),'r')
% plot(fmin,1,'x')