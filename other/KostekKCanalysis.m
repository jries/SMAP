% analyse data for Kostek

path='/Volumes/t2ries/users/Kostek/For Jonas/Distance data';
if 0
dirs=dir(path);
ind=1;
clear dists pname meanall stdall semall numdists
for k=1:length(dirs)
    if dirs(k).isdir && ~strcmp(dirs(k).name(1),'.')
        dirh=dirs(k).name;
        disp(dirh)
        dists{ind}=getdistances([path filesep dirh filesep]);
        pname{ind}=strrep(dirh,'_','-');
        ind=ind+1;
    end
end
else
    load([path filesep 'distances.mat'])
end

for k=1:length(pname)
    meanall(k)=mean(dists{k},'omitnan');
    stdall(k)=std(dists{k},'omitnan');
    numdists(k)=length(dists{k});
    semall(k)=stdall(k)/sqrt(numdists(k));
end


figure(88);
hold off
dd=0.4;
n=1:length(dists);
errorbar(meanall,n,stdall,'k+','horizontal','MarkerSize',12)
ax=gca;
ax.YTick=n;
ax.YTickLabel=pname;
ylim([0 n(end)+1])

hold on
for k=1:length(pname)
rectangle('Position',[meanall(k)-semall(k),n(k)-dd,semall(k)*2,2*dd])
end


%p-values
for k=1:length(pname)
    for l=k+1:length(pname)
        [ptrure(k,l) ,ptt(k,l)]=ttest2(dists{k},dists{l});
        if ptt(k,l)<0.001
            pstar(k,l)=3;
        elseif ptt(k,l)<0.01
            pstar(k,l)=2;
        elseif ptt(k,l)<0.05
            pstar(k,l)=1;
        else
            pstar(k,l)=0;
        end
    end
end
pstar
txtt=num2str(ptt,2)

function dists=getdistances(file)
dirh=dir([file filesep '*.mat']);
fh=[file dirh(1).name];
l=load(fh);
sites=l.saveloc.siteexplorer.sites;
dists=[];
for k=length(sites):-1:1
    si=sites(k);
    dc=(si.annotation.line2.pos(2,:)-si.annotation.line2.pos(1,:))*1000;
    dc=-dc; %invert axis
    ds=(si.annotation.line1.pos(2,:)-si.annotation.line1.pos(1,:))*1000;
    len(k)=norm(ds);
    lencd(k)=sum(dc.*ds)/norm(dc);
    
%     arel(k)=acos(lencd(k)/norm(ds)/1000);
    dxsite(k)=ds(1); dysite(k)=ds(2);
    arel(k)=(atan2(dc(2),dc(1))-atan2(ds(2),ds(1)));
    dists(k)=lencd(k);
end
if isempty(dists)
    
end
end