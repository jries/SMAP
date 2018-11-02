function histogram=cluster_mmaple_makehist(par,cluster)

numcluster=length(cluster);
for k=1:numcluster
    incluster=~removesingles(cluster(k),25, 2,1); %remove those which are far from cluster. Alternative: maximum distacne from center
    cluster(k).incluster= incluster;
%     incluster=cluster(k).incluster;
    cluster(k).meansigmapsf=mean(cluster(k).psf(incluster));
    cluster(k).meanbg=mean(cluster(k).bg(incluster));
    cluster(k).meanphot=mean(cluster(k).phot(incluster));
    cluster(k).meanlocp=mean(cluster(k).locprec(incluster));
    cluster(k).stdx=std(cluster(k).x(incluster));
    cluster(k).stdy=std(cluster(k).y(incluster));
    cluster(k).locsdf=min(sum(diff(cluster(k).frame(incluster))>par.dftime)+1,sum(incluster));
%     cluster(k).glocs=min(sum(diff(cluster(k).frame(incluster))>1)+1,sum(incluster));
    cluster(k).locsdfv=cluster(k).locsdf+0*cluster(k).x;
%     cluster(k).glocs=cluster(k).glocsreduced; %XXXXXXX
%     cluster(k).issinglecluster;
end



switch par.c_groupfield.Value
    case 1 
        gfield='locs';
    case 2
        gfield='glocs';
    case 3
        gfield='glocsreduced';
    case 4
        gfield='locsdf';
end
sigmapsf=[cluster(:).meansigmapsf]';
% clocs=[cluster(:).locsdf]';
clocs=[cluster(:).(gfield)]';
ig=clocs>0&sigmapsf>0;

ig=ig&[cluster(:).issinglecluster]';
%Manuel: remove
% ig=ig&sigmapsf>120&sigmapsf<140;

cluster(~ig)=[];
initaxis(par.resultstabgroup,'scatter');
%     par.resultstabs(3).Title='scatter';
%     ax3=axes('Parent',par.resultstabs(3));
%     axes(ax3)
subplot(2,3,1);
try
dscatter((([cluster.meansigmapsf]')),[cluster.(gfield)]')
hold on
n=min([cluster.meansigmapsf]):5:max([cluster.meansigmapsf]);
plot(n,bindata([cluster.meansigmapsf],[cluster.(gfield)],n,'median'))
plot(n,bindata([cluster.meansigmapsf],[cluster.(gfield)],n,'mean'))
% plot(n,bindata([cluster.meansigmapsf],[cluster.(gfield)],n,'geomean'))
hold off
xlabel('sigmapsf')
ylabel('locs')

subplot(2,3,2);
dscatter(sqrt(([cluster.stdx].*[cluster.stdy])'),[cluster.(gfield)]')
xlabel('cluster size')
ylabel('locs')
hold on
n=min(sqrt(([cluster.stdx].*[cluster.stdy])')):2:max(sqrt(([cluster.stdx].*[cluster.stdy])'));
plot(n,bindata(sqrt(([cluster.stdx].*[cluster.stdy])'),[cluster.(gfield)]',n,'median'))
hold off

subplot(2,3,3)
% dscatter((([cluster.meanbg]')),[cluster.(gfield)]')
% xlabel('bg')
% ylabel('locs')
dscatter((([cluster.meanlocp]')),[cluster.meansigmapsf]')
xlabel('locp')
ylabel('sigmapsf')
hold on
n=min([cluster.meanlocp]'):2:max([cluster.meanlocp]');
plot(n,bindata([cluster.meanlocp]',[cluster.meansigmapsf]',n,'median'))
hold off

subplot(2,3,4)
dscatter((([cluster.meanlocp]')),[cluster.(gfield)]')
xlabel('locp')
ylabel('locs')
hold on
n=min([cluster.meanlocp]'):2:max([cluster.meanlocp]');
plot(n,bindata([cluster.meanlocp]',[cluster.(gfield)]',n,'median'))
hold off


subplot(2,3,5)
dscatter(sqrt(([cluster.stdx].*[cluster.stdy])'),[cluster.meansigmapsf]')
xlabel('cluster size')
ylabel('sigmapsf')
hold on
n=min(sqrt(([cluster.stdx].*[cluster.stdy])')):2:max(sqrt(([cluster.stdx].*[cluster.stdy])'));
plot(n,bindata(sqrt(([cluster.stdx].*[cluster.stdy])'),[cluster.meansigmapsf]',n,'median'))
hold off

subplot(2,3,6)
dscatter([cluster.locs]',[cluster.(gfield)]')

xlabel('locs')
ylabel('locs selected')
catch
end
% %
% X=[[cluster(ig).locs]',[cluster(ig).meanlocp]',[cluster(ig).meanbg]',[cluster(ig).stdx]',[cluster(ig).stdy]',[cluster(ig).meansigmapsf]'];%,cluster(ig).meanbg',cluster(ig).stdx',cluster(ig).stdy',cluster(ig).meansigmapsf'];
% % [p,score]=pca(X);
% % figure(123)
% % plot(score(:,1),score(:,2),'.')
% idx=kmeans(X,5,'Replicates',10);
% figure(123)
% col=jet(max(idx));
% scatter(X(:,5),X(:,1),12,col(idx,:),'+')


% %filter
minsigma=par.c_PSFmin;
maxsigma=par.c_PSFmax;
sigmin=par.c_stdymin;
sigmax=par.c_stdymax;
% maxbg=150;

ig=[cluster(:).meansigmapsf]>minsigma & [cluster(:).meansigmapsf]<maxsigma;
ig=ig&sqrt([cluster(:).stdx].*[cluster(:).stdy])>sigmin & sqrt([cluster(:).stdx].*[cluster(:).stdy])<sigmax;

dlocs=[cluster(ig).(gfield)]';

initaxis(par.resultstabgroup,'histogram');
maxhist=myquantile(dlocs,0.998);
% figure(43)
[y,x]=hist(dlocs,0:maxhist);
xlim([1,100]);

plot(x,y)
title(length(dlocs));
% sum(y)
% length(cluster)

% y=y/sum(y(:));
histogram.h=y;
histogram.c=x;

end
