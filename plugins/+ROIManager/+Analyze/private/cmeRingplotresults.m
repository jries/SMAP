function fh=cmeRingplotresults(p,results)
maxdrfac=p.maxdrro;
markersize=4;
N=results.N;
ro=results.ro;
rc=results.rc;
dr=results.dr;
sigma=results.sigma;
% rdensity=results.rdensity;
filenumber=results.filenumber;
filenumberrange=results.filenumberrange;

drrel=dr./ro;
% drrelmax1=drrel;drrelmax1(drrelmax1>1)=1;
drrelmax=drrel;drrelmax(drrelmax>maxdrfac)=maxdrfac;
drri=dr;drri(drri>ro)=ro(drri>ro);


fh=figure(90);
subplot(3,4,1)
hold off
pout=fithist(N,10,2); 
title(['max ',num2str(pout.max,'%3.0f'),', mean ',num2str(pout.mean,'%3.0f'),', median ',num2str(pout.median,'%3.0f'),', std ',num2str(pout.std,'%3.0f')]);

xlabel('N')
subplot(3,4,5)
plot(N,ro,'.','MarkerSize',markersize)
xlabel('N')
ylabel('ro')
% hold on
% plot(N,drri,'o')
% hold off
% legend('outer ring radius','ring thickness')
subplot(3,4,6)
plot(N,drrelmax,'.','MarkerSize',markersize)
xlabel('N')
ylabel('dr/rout')

subplot(3,4,7)
plot(ro,drrelmax,'.','MarkerSize',markersize)
xlabel('rout')
ylabel('dr/rout')

subplot(3,4,10)
% plot(mean(cat(3,rdensity{:}),3))
       hold off
       rdist=results.sumrdensityn(2)-results.sumrdensityn(1);
       rnh=(0:length(results.sumrdensity1)-1)*rdist;
if isfield(results, 'rdensity1')
    rall=cell2mat(results.rdensity1');
    stdr=std(rall,1);
    mr=mean(rall,1);

       fill([rnh(2:end),rnh(end:-1:2)],[mr(2:end)-stdr(2:end),mr(end:-1:2)+stdr(end:-1:2)],'c')
        hold on
    plot(rnh(2:end),mr(2:end)+stdr(2:end),rnh(2:end),mr(2:end)-stdr(2:end))
 

end


plot(rnh,results.sumrdensity1/results.numsites,rnh,results.sumrdensity2/results.numsites)
xlabel('radius (nm)')
ylabel('density')
title('radial density')


subplot(3,4,2)
% hori=histogram(ro);
% binedges=hori.BinEdges;
% [h,pd]=myhistfit(ro);

hold off
pout=fithist(ro,2,4); 
title(['max ',num2str(pout.max,'%3.0f'),', mean ',num2str(pout.mean,'%3.0f'),', median ',num2str(pout.median,'%3.0f'),', std ',num2str(pout.std,'%3.0f')]);

xlabel('outer radius')
% title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));

subplot(3,4,3)
hold off
pout=fithist(drrelmax,.05,4); 
title(['max ',num2str(pout.max,'%3.2f'),', ratio>1/<1: ' num2str(sum(drrelmax>1)/sum(drrelmax<1),'%3.2f')]);

% histogram(drrelmax);
% [~,pd]=myhistfit(drrelmax);
% h=histfit(drrelmax,20,'kernel');
% title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));
% [~,ind]=max(h(2).YData);
% maxX=h(2).XData(ind); 
% title(strcat('max ',num2str(maxX)));
xlabel('dr/rout');

subplot(3,4,9)
plot(filenumber,N,'+','MarkerSize',markersize)
hold on
plot(filenumberrange,results.Nfilemean,'ok','MarkerSize',markersize)
plot(filenumberrange,results.Nfilemedian,'*k','MarkerSize',markersize)
hold off
xlabel('filenumber');
ylabel('N');
legend('N','Nmean','Nmedian')
subplot(3,4,4)
q=myquantilefast(rc,0.995);
rca=rc(rc<q);
hold off
pout=fithist(rca,2,4); 
title(['max ',num2str(pout.max,'%3.0f'),', mean ',num2str(pout.mean,'%3.0f'),', median ',num2str(pout.median,'%3.0f'),', std ',num2str(pout.std,'%3.0f')]);

ylim([0 q])
% histogram(rc,hori.BinEdges)
% [h,pd]=myhistfit(ro,hori.BinEdges);
% [~,pd]=myhistfit(rc);
% title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));
xlabel('radius circfit')

subplot(3,4,8)
hold off
pout=fithist(sigma,1,4); 
title(['max ',num2str(pout.max,'%3.1f'),', mean ',num2str(pout.mean,'%3.1f'),', median ',num2str(pout.median,'%3.1f'),', std ',num2str(pout.std,'%3.1f')]);

% histogram(sigma);
% [h,pd]=myhistfit(sigma);
% title(strcat(num2str(pd.mu),'±',num2str(pd.sigma)));
xlabel('sigma imfit')

axim=subplot(3,4,12);
axis(axim,'equal')
axis(axim,'tight')
% xs=rc;ys=hori.BinEdges;
hdt=datacursormode(fh);
set(hdt,'UpdateFcn',{@CMElabel,results.sitenames,results.images,axim})

subplot(3,4,11);
pim=results.sumimage;
s=size(pim);
if length(s)==3
npim=squeeze(max(max(pim)));
pim(:,:,1)=pim(:,:,1)/npim(1);
pim(:,:,2)=pim(:,:,2)/npim(2);
pim(:,:,3)=pim(:,:,3)/npim(3);
end
imagesc(pim)
title('average image')
axis equal
axis tight
end

function pout=fithist(data,binw,fitrange)
if nargin<3
    fitrange=4;
end
mind=min(data);maxd=max(data);
binw=min(binw,(maxd-mind)/10);

range=floor(min(data)):binw:max(data);
h=histogram(data,range);

pout.median=median(data);
pout.max=getmax(h,fitrange);
pout.mean=mean(data);
pout.std=std(data);

end

function m=getmax(h,fitrange)
[~,ind]=max(h.Values(2:end-1));
ind=ind+1;
range=max(1,ind-fitrange):min(ind+fitrange,length(h.BinEdges)-1);
db=h.BinEdges(2)-h.BinEdges(1);
bincenters=h.BinEdges(range)+db/2;
values=h.Values(range);
fp=mygaussfit(bincenters,values,[max(values),h.BinEdges(ind)+db/2, db*fitrange,0],true);
m=fp(2);
hold on
plot(bincenters,mygaussforfit(fp,bincenters),'k')

end