function out=cluster_segment_NPC(locs,par)
% profile on
% global cluster hall hindex 

diameterNPC=100;
lengthscale=par.lengthscale;
% loc=locs;
ploton=par.showresults;
psfsigmafocus=par.segment_maxPSF;
pixrec=5;
gaussfac=1;
gaussmin=5; %clusters are starting points for meanshiftpeakfind

indInFocus=locs.PSFxnm<psfsigmafocus;
xmin=min(locs.xnm);ymin=min(locs.ynm);

scale=median(locs.locprecnm);

posa.x=locs.xnm-xmin;
posa.y=locs.ynm-ymin;
posa.s=locs.locprecnm*0.5;
xmax=max(posa.x);ymax=max(posa.y);

% if 0
posrem.x=locs.xnm(indInFocus)-xmin;
posrem.y=locs.ynm(indInFocus)-ymin;
posrem.s=max(locs.locprecnm(indInFocus)*gaussfac,gaussmin);
% posrem.s=posrem.s*0+pixrec/2;
imsrrem=gaussrender(posrem,[0 xmax],[0 ymax],pixrec, pixrec);

% cutoff=2;
% imbw=double(imsrrem>cutoff);
% imsrrem(imsrrem>cutoff)=cutoff;
rim=10;
rRingO=double((diameterNPC/2+rim)/pixrec)
hfilterO=fspecial('disk',rRingO);
 rRingI=double(max(1,(diameterNPC/2-rim)/pixrec))
hfilterI=fspecial('disk',rRingI);

si=size(hfilterI);
n=(si(1)-1)/2;
so=size(hfilterO);
c=(so(1)+1)/2;
hfilterIb=0*hfilterO;
hfilterIb(c-n:c+n,c-n:c+n)=hfilterI;
hfilter=hfilterO/max(hfilterO(:))-hfilterIb/max(hfilterIb(:));
hfilter=hfilter/sum(hfilter(:));

% imfO=imfilter(imsrrem,hfilterO);
% imfI=imfilter(imsrrem,hfilterI);
imfD=imfilter(sqrt(imsrrem),hfilter);

hfilterGauss=fspecial('gauss',21,6);
hfilterGauss2=fspecial('gauss',50,2.5*rRingO);
% imfD=imfilter(imfD,hfilterGauss)-imfilter(imfD,hfilterGauss2);
figure(88);imagesc(imfD)
figure(89);imagesc(imsrrem)

imsr=gaussrender(posa,[0 xmax],[0 ymax],pixrec, pixrec);
%spatial scale from locprec
maxima=NMS2DBlockCcall(single(imsrrem),3);
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
[beadnum,numlocs,dist,numlocsred]=associatenearest(maximaMS(:,1),maximaMS(:,2),posa,scale*5);
if ploton
    recgui.initaxis(par.resultstabgroup,'peaks on SR')
%     par.resultstabs(1).Title='peaks on SR';
%     ax=axes('Parent',par.resultstabs(1));
%     axes(ax);
    assigned=dist<scale*5;
    impl=imsr/myquantile(imsr(:),0.999);
    impl(impl>1)=1;
    imagesc(impl);
    title('filtered reconstructed image, only in focus')    
    hold on
        plot(maxima(:,2),maxima(:,1),'k+')
        plot(maximaMS(:,1)/pixrec,maximaMS(:,2)/pixrec,'wo')
    hold off
        col=jet(max(beadnum));
        
    recgui.initaxis(par.resultstabgroup,'clusters')   
%     ax2=axes('Parent',par.resultstabs(2));
%     par.resultstabs(2).Title='clusters';
%     axes(ax2);
% colnum=beadnum(assigned)
colnum=ceil(rand(sum(assigned),1)*max(beadnum));
    scatter(posa.x(assigned),posa.y(assigned),12,col(colnum(beadnum(assigned)),:),'+')
    hold on 
    scatter(posa.x(~assigned),posa.y(~assigned),8,'k.');
    scatter(maximanm(:,2),maximanm(:,1),'ko')
    scatter(maximaMS(:,1),maximaMS(:,2),'kx')
    hold off
     axis ij
end


%make cluster structure for further analysis
% x, y, frame, sigma, locprec, incluster
% stat.stdx, stdy, meansigma
cluster=[];
for k=1:max(beadnum)
    indbead=beadnum==k;
    cluster(k).x=locs.xnm(indbead);
    cluster(k).y=locs.ynm(indbead);
    cluster(k).bg=locs.bg(indbead);
    cluster(k).phot=locs.phot(indbead);
  
    
    cluster(k).frame=locs.frame(indbead);
    cluster(k).psf=locs.PSFxnm(indbead);
    cluster(k).locprec=locs.locprecnm(indbead);
    
    cluster(k).locs=numlocsred(k);
   
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



