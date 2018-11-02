function out=cluster_counting_2(locData,par)
% profile on
% global cluster hall hindex 

lengthscale=par.lengthscale;
% loc=locs;
ploton=par.showresults;
psfsigmafocus=par.segment_maxPSF;

locs=locData.getloc({'frame','xnm','ynm','phot','bg','PSFxnm','locprecnm','groupindex'},'layer',1,'position','roi','grouping','ungrouped');




pixrec=5;
gaussfac=1;
gaussmin=5; %clusters are starting points for meanshiftpeakfind
gaussmin=lengthscale/2;
locprecmax=20;
inlp=locs.locprecnm<locprecmax;
indInFocus=locs.PSFxnm<psfsigmafocus;
xmin=min(locs.xnm);ymin=min(locs.ynm);


badind=removesingles(struct('x',locs.xnm,'y',locs.ynm),lengthscale, 3,3);
%min neighbours, iterations
%maybe that is too much already? Although removing 
incluster=~badind;
indInFocus=indInFocus&inlp&(incluster);
% scale=median(locs.locprecnm);

posa.x=locs.xnm(incluster)-xmin;
posa.y=locs.ynm(incluster)-ymin;
posa.s=locs.locprecnm(incluster)*0.5;
xmax=max(posa.x);ymax=max(posa.y);

% if 0
posrem.x=locs.xnm(indInFocus)-xmin;
posrem.y=locs.ynm(indInFocus)-ymin;
posrem.s=max(locs.locprecnm(indInFocus)*gaussfac,10);
imsrrem=gaussrender(posrem,[0 xmax],[0 ymax],pixrec, pixrec);

% imsr=gaussrender(posa,[0 xmax],[0 ymax],pixrec, pixrec);
%spatial scale from locprec
maxima=maximumfindcall(single(imsrrem));
maxima(maxima(:,3)<0.00001,:)=[];
maximanm=maxima*pixrec; %in nanometers
% else

% maximaMS=meanshiftpeakfind(maximanm(:,2),maximanm(:,1),posa.x,posa.y,lengthscale);

% maximaMS=kmeansMeanshift(posa.x,posa.y,lengthscale);
% end

%remove all clusteres which are too dense. 
% indout=false(length(maxima),1);
% k=1;
% while k<=length(maxima)
%      r2=(maxima(k,1)-maxima(:,1)).^2+(maxima(:,2)-maxima(k,2)).^2;
%      if sum(r2<=4*scale^2)>1 %more than 1 maximum in 4 pixels radius (=20 nm) 
%         maxima(k,:)=[];
%         k=k-1;
% %         indout(k)=true;
%      end
%     k=k+1;
% end
% maxima(indout,:)=[];
% 


% [beadnum,numlocs,dist,numlocsred]=associatenearest(maximanm(:,2),maximanm(:,1),posa,scale*5);
[beadnum,numlocs,dist,numlocsred,inreduced]=associatenearest(maximanm(:,1),maximanm(:,2),posa,lengthscale);
mindistance2=lengthscale^2;
goodclusters=true(size(maximanm,1),1);
for k=1:size(maximanm,1)
    d=(maximanm(k,1)-maximanm(k+1:end,1)).^2+(maximanm(k,2)-maximanm(k+1:end,2)).^2;
    badind=find(d<mindistance2);
    if ~isempty(badind)
        goodclusters(k)=false;
        goodclusters(badind+k)=false;
    end
    
end

%identify clusters that are not close to ohter clusters, only use those
%values

if ploton
    initaxis(par.resultstabgroup,'peaks on SR');
%     par.resultstabs(1).Title='peaks on SR';
%     ax=axes('Parent',par.resultstabs(1));
%     axes(ax);
    assigned=dist<1.5*lengthscale;%+2*locs.locprecnm(incluster);
    
    beadnum(~assigned)=0;
%     assigned=dist<3*locs.locprecnm;
    impl=imsrrem/myquantile(imsrrem(:),0.999);
    impl(impl>1)=1;
    imagesc(impl);
    title('filtered reconstructed image, only in focus')    
    hold on
        plot(maxima(:,1),maxima(:,2),'k+')
%         plot(maximaMS(:,1)/pixrec,maximaMS(:,2)/pixrec,'wo')
    hold off
        col=jet(max(beadnum));
        
    initaxis(par.resultstabgroup,'clusters')  ; 
%     ax2=axes('Parent',par.resultstabs(2));
%     par.resultstabs(2).Title='clusters';
%     axes(ax2);
% colnum=beadnum(assigned)
colnum=ceil(rand(sum(assigned),1)*max(beadnum));
    scatter(posa.x(assigned),posa.y(assigned),12,col(colnum(beadnum(assigned)),:),'+')
    hold on 
    scatter(posa.x(~assigned),posa.y(~assigned),8,'k.');
    scatter(maximanm(:,1),maximanm(:,2),'kx')
%     scatter(maximaMS(:,1),maximaMS(:,2),'kx')
    hold off
     axis ij
end


%make cluster structure for further analysis
% x, y, frame, sigma, locprec, incluster
% stat.stdx, stdy, meansigma
inclusterf=find(incluster);
cluster=[];
   initaxis(par.resultstabgroup,'found clusters')
hold off
for k=1:max(beadnum)
    indbead=beadnum==k;
    cluster(k).x=locs.xnm(inclusterf(indbead));
    cluster(k).y=locs.ynm(inclusterf(indbead));
    cluster(k).bg=locs.bg(inclusterf(indbead));
    cluster(k).phot=locs.phot(inclusterf(indbead));
    cluster(k).inreduced=inreduced((indbead));
  
    
    cluster(k).frame=locs.frame(inclusterf(indbead));
    cluster(k).psf=locs.PSFxnm(inclusterf(indbead));
    cluster(k).groupindex=locs.groupindex(inclusterf(indbead));
    
    cluster(k).locprec=locs.locprecnm(inclusterf(indbead));
    
    cluster(k).locs=numlocsred(k);
    cluster(k).glocs=numel(unique(cluster(k).groupindex));
    cluster(k).glocsreduced=numel(unique(cluster(k).groupindex(cluster(k).inreduced)));
    
    cluster(k).issinglecluster=goodclusters(k);
    plot(cluster(k).x,cluster(k).y,'.')
    hold on
   
end



cluster(1).xmax=xmax;
cluster(1).ymax=ymax;
% cluster(1).name=
% [~,cluster(1).name]=fileparts(par.file);
% if isempty(hindex)
%     hindex=1
% end
% hall(hindex).cluster=cluster;

out=cluster;

%mean shift clustering
% x=[locs.xnm, locs.ynm]';
% bandWidth=40;
% plotFlag=0;
% % [clustCent,data2cluster,cluster2dataCell] = MeanShiftCluster(dataPts,bandWidth,plotFlag);
% [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(x,bandWidth);
% 
% figure(10),clf,hold on
% numClust = length(clustMembsCell);
% cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk'; cVec = [cVec cVec];cVec = [cVec cVec];cVec = [cVec cVec];
% for k = 1:min(numClust,length(cVec))
%     myMembers = clustMembsCell{k};
%     myClustCen = clustCent(:,k);
%     plot(x(1,myMembers),x(2,myMembers),[cVec(k) '.'])
%     plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k), 'MarkerSize',10)
% end
% title(['no shifting, numClust:' int2str(numClust)])

% parout=cluster_mmaple_evaluate(par);



