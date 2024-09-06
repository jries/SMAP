% FP blinking, later make plugin out of it
maxdTs= 5; %s
exposuret=g.locData.files.file(1).info.exposure/1000; %s
maxdTf=maxdTs/exposuret;
locd=g.locData.copy;
disp('regrouping with large dT, might take time')
locd.regroup([30 100],maxdTf,2);
locd.sort('groupindex');
%%
gi=locd.loc.groupindex;
maxg=gi(end);
ind1=1;
gih=gi(1);
newtrace=false(length(gi),1); newblink=false(length(gi),1);
ontime=zeros(length(gi),1);offtime=zeros(length(gi),1);
tracetime=zeros(length(gi),1);numberOfBlinks=zeros(length(gi),1);onfraction=zeros(length(gi),1);
dframe=diff(locd.loc.frame);
for k=1:length(gi) 
    if gi(k)>gih || k==length(gi) % new track. we lose last localization...
        indh=ind1:k-1;
        framesh=locd.loc.frame(indh);
        newtrace(ind1)=true;
        tracetime(ind1)=framesh(end)-framesh(1)+1;
        onfraction(ind1)=length(framesh)/tracetime(ind1);
        % numberOfBlinks(ind1)=sum(diff(framesh)>1)+1;
        dfh=[diff(framesh); -1];
        dfh1=[true;dfh>1];
        numberOfBlinks(ind1)=sum(dfh1);
        fdfh=[find(dfh1); length(framesh)+1];
        for l=1:length(fdfh)-1
            indbl=indh(fdfh(l));
            newblink(indbl)=true;

            ontime(indbl)=fdfh(l+1)-fdfh(l);
            offtime(indbl)=dfh(fdfh(l+1)-1);
            

        end

        findh=find(indh);
        gih=gi(k);
        ind1=k;
    end
end
%%
ysc='log';
ysc='linear';
figure(88)
subplot(2,3,1)
histogram(ontime(newblink),1:max(ontime),"DisplayStyle","stairs",Normalization="probability")
xlabel('on-times (frames)')
ax=gca; ax.YScale=ysc;ax.XLim(1)=0;
subplot(2,3,2)
offt=offtime(newblink);
histogram(offt(offt>0)*exposuret,"DisplayStyle","stairs",Normalization="probability")
xlabel('off-times (s)')
ax=gca; ax.YScale=ysc;ax.XLim(1)=0;
subplot(2,3,4)
histogram(tracetime(newtrace)*exposuret,0:exposuret:10,"DisplayStyle","stairs",Normalization="probability")
xlabel('duration of trace (s)')
ax=gca; ax.YScale=ysc;ax.XLim(1)=0;ax.XLim(2)=quantile(tracetime(newtrace)*exposuret,.9);
subplot(2,3,5)
histogram(onfraction(newtrace&numberOfBlinks>1),"DisplayStyle","stairs",Normalization="probability")
xlabel('on fraction')
subplot(2,3,6)
histogram(numberOfBlinks(newtrace)-1,0:quantile(numberOfBlinks(newtrace)-1,.99),"DisplayStyle","stairs",Normalization="probability")
xlabel('number of blinks')
ax=gca; ax.YScale=ysc;ax.XLim(1)=0;


%%
dT=[0 2 10 100 1000 10000];
% dT=0
locd2=g.locData.copy;
figure(223)
clf
clear mp dTtxt
for k=1:length(dT)
    locd2.regroup([30 100],dT(k),2);
    subplot(2,2,1)
    pmax=15000;
    histogram(locd2.grouploc.phot,0:100:pmax,"DisplayStyle","stairs",Normalization="probability")
    xlabel('photons')
    hold on
    xlim([0 pmax])
    subplot(2,2,2)
    histogram(locd2.grouploc.numberInGroup,1:50,"DisplayStyle","stairs",Normalization="probability")
    xlabel('on-time (frames)')
    hold on
    dTtxt{k}=num2str(dT(k));
    mp(k)=mean(locd2.grouploc.phot);
end
legend(dTtxt)
subplot(2,2,3)
plot(dT,mp)
ax=gca; ax.XScale="log";
xlabel('dT (frames)')
ylabel('mean photon numbers')
