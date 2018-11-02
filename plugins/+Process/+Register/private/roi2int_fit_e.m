function p=roi2int_fit_e(roi1,x,y,PSFxpix,PSFypix,roisize,bgroi)
%weights not implemented? Do htat!
sim=size(roi1);
if length(sim)==2
    sim(3)=1;
end
mp=round(sim+1)/2;
dn=min(min(mp)+1,round((roisize+1)/2));
[X,Y]=meshgrid(-dn:dn);
p=zeros(sim(3),1);

if nargin<7||isempty(bgroi)
p=zeros(sim(3),2);
for k=1:sim(3)

    gauss=exp((-(x(k)-X).^2)/2/PSFxpix(k)^2-((y(k)-Y).^2)/2/PSFypix(k)^2)/pi/PSFxpix(k)/PSFypix(k)/2;
    weights=sqrt(gauss);
    Xmat=horzcat(gauss(:), gauss(:)*0+1);
    roih=roi1(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
    p(k,:)=Xmat\roih(:);
    if 0%p(k,1)>2500
        p(k,:)
        figure(66)
        subplot(2,2,1)
        imagesc(-dn:dn,-dn:dn,roih);
        hold on
        plot(x(k),y(k),'+')
        hold off
        subplot(2,2,2);
        imagesc(-dn:dn,-dn:dn,gauss*p(k,1)+p(k,2))
        hold on
        plot(x(k),y(k),'+')
        hold off
        subplot(2,2,3);
        imagesc(-dn:dn,-dn:dn,gauss*p(k,1)+p(k,2)-roih)
        waitforbuttonpress
    end
end


else %fit bg


for k=1:sim(3)
    gauss=exp((-(x(k)-X).^2)/2/PSFxpix(k)^2-((y(k)-Y).^2)/2/PSFypix(k)^2)/pi/PSFxpix(k)/PSFypix(k)/2;
    weights=sqrt(gauss);
    Xmat=horzcat(gauss(:));
    bgh=bgroi(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k);
    roih=roi1(mp(1)-dn:mp(1)+dn,mp(2)-dn:mp(2)+dn,k)-bgh;
    p(k,1)=Xmat\roih(:);
    p(k,2)=sum(sum(bgh))/(2*dn+1)^2;
end
end
