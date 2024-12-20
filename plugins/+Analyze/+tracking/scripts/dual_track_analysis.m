% analyze dual-color tracks

%filters
minlenframes=5;
maxd=200; % nm max distance to associate localizations

cotracklength=4; %minimum data poits associated between tracks
cotrackfraction = 0.0; %minimum fraction of localizations associated

% test directed movement
aspectratio=1; %should be below this for directed motion.
lennmstartend=300; % minimum endpoint- startpoint
lennmmin=0; %minimum of largest extenstion standard deviation (not so useful)
minvelocity=200; %nm/s
maxvelocity=1500;

% visualizing co-tracks
onlyprogressivecotracks=true; %only when both partners are progressive
markintiffile=true;
tiffile=fout;
tiffile='/Users/ries/datalocal/2color_kinesin/25_50ms_561nm01_640nm02_600w52_676w37_1_MMStack_Default_combined.tif';

obj=g;
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','xpix','ypix','frame','track_id','channel','track_length','layer','filenumber'},'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','ungrouped');


exposuretime=obj.locData.files.file(locs.filenumber(1)).info.exposure;

% XXX if exposure is wront, overwrite here
% exposuretime= 10

% test directed movement, calculate statistics
usetracks=unique(locs.track_id(locs.track_id>0));
trackstat.lennmstartend=zeros(max(usetracks),1);
trackstat.stdshort=zeros(max(usetracks),1);
trackstat.stdlong=zeros(max(usetracks),1);
trackstat.lenframe=zeros(max(usetracks),1);
trackstat.partnertrackid=zeros(max(usetracks),1);
trackstat.velocity=zeros(max(usetracks),1);
locs.track_length_new=zeros(size(locs.xnm));

for k=1:length(usetracks)
    iduset=usetracks(k);
    tind=find(locs.track_id==iduset);
    trackstat.lenframe(iduset)=length(tind);
    locs.track_length_new(tind)=trackstat.lenframe(iduset);
    if trackstat.lenframe(iduset)<2
        continue
    end


    xh=locs.xnm(tind);yh=locs.ynm(tind);fh=locs.frame(tind);
    trackstat.lennmstartend(iduset)=sqrt((xh(end)-xh(1)).^2+(yh(end)-yh(1)).^2);
        
    c = cov(xh-mean(xh), yh-mean(yh));
    [a, ev] = eig(c);
    [ev,ind] = sort(diag(ev));
    [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
    trackstat.angle(iduset)=cart2pol(xa, ya);
    trackstat.stdshort(iduset)=real(sqrt(ev(1)));
    trackstat.stdlong(iduset)=real(sqrt(ev(2)));
    trackstat.velocity(iduset)=trackstat.lennmstartend(iduset)/(exposuretime/1000*trackstat.lenframe(iduset)); %nm/s
           
end

intrack=(locs.track_id>0);

longtracks=locs.track_length_new>minlenframes;


inchannel=locs.channel==1;
sum(inchannel)


indr=find((intrack&longtracks&inchannel));
locr.x=locs.xnm(indr);
locr.y=locs.ynm(indr);
locr.frame=locs.frame(indr);
locr.track_id=locs.track_id(indr);


inchannel2=locs.channel==2;
sum(inchannel2)
indt=find((intrack&longtracks&inchannel2));
loct.x=locs.xnm(indt);
loct.y=locs.ynm(indt);
loct.frame=locs.frame(indt);
loct.track_id=locs.track_id(indt);

[iAa,iBa,nA,nB,nseen]=matchlocsall(locr,loct,0,0,maxd);

trackstat.partnertrackid(locs.track_id(indr(iAa)))=(locs.track_id(indt(iBa)));
trackstat.partnertrackid(locs.track_id(indt(iBa)))=(locs.track_id(indr(iAa)));

% pairedr=false(size(locs.xnm));
% pairedt=pairedr;
partner=zeros(size(locs.xnm));
% pairedr(indr(iAa))=true;
% pairedt(indt(iBa))=true;
partner(indr(iAa))=indt(iBa);
partner(indt(iBa))=indr(iAa);

usetracks=unique(locs.track_id(intrack&longtracks));

roi=obj.locData.files.file(locs.filenumber(1)).info.roi;
pixelsize=obj.locData.files.file(locs.filenumber(1)).info.cam_pixelsize_um*1000;


%look for progressive movement
validstats=true(size(trackstat.velocity));
validstats = validstats & trackstat.stdshort./trackstat.stdlong<aspectratio;
validstats = validstats & trackstat.lennmstartend > lennmstartend;
validstats = validstats & trackstat.stdlong > lennmmin;
validstats = validstats & trackstat.velocity > minvelocity;
validstats = validstats & trackstat.velocity < maxvelocity;

trackstat.progressive=validstats;

comovement=false(size(trackstat.velocity));
progressivepartner=false(size(trackstat.velocity));
trackstat.channel=zeros(size(trackstat.velocity));


for k=1:length(usetracks)  
    idh=usetracks(k);
    indtr=locs.track_id==idh;

    %look at co-movement
    pid=0;
    comovement(idh)=trackstat.partnertrackid(idh)>0;
    if comovement(idh) %&& validstats
        partnerind=partner(indtr);
        partnerids=locs.track_id(partnerind(partnerind>0));
        [pid,npart]=mode(partnerids);
        trackstat.partnerids(idh)=pid;
        
        lenpartner=sum(locs.track_id==pid);
        minduallength=min(sum(indtr),lenpartner);

        comovement(idh)=comovement(idh) & npart>cotracklength;
        comovement(idh)=comovement(idh) & npart/minduallength>cotrackfraction;
        progressivepartner(idh)=trackstat.progressive(pid);
    end
    trackstat.channel(idh)=mode(locs.channel(indtr));
    
end

trackstat.comovement=comovement;
trackstat.progressivepartner=progressivepartner;


% plot tracks
figure(88)
hold off
ax=gca;

cols=[1 0 1
      0 1 1
      1 0 0 
      0 0 1
      .5 0.2 0.2
      0.2 0.2 .5
      .7 0 0.7
      0 0.7 .7];

colind=trackstat.channel+trackstat.comovement*2+(~trackstat.progressive)*4;   

for k=1:length(usetracks)  
    idh=usetracks(k);
    indtr=locs.track_id==idh;
    
    msize=3;
    lw=1;
    symb='-';

    if trackstat.comovement(idh)
        if trackstat.channel(idh)==2
                symb='x-';
                msize=3;
        else
                symb='+-';
                msize=7;
                lw=2;            
        end
    end
    hp=plot(ax,locs.xnm(indtr)/pixelsize(1)-roi(1),locs.ynm(indtr)/pixelsize(2)-roi(2),symb,'Color',cols(colind(idh),:),'LineWidth',lw,'Tag','test','MarkerSize',msize);
    hold(ax,"on")
    pidlabel=0*locs.track_id(indtr)+pid;
    dtRows = [dataTipTextRow("frame",locs.frame(indtr)),...
    dataTipTextRow("ID",locs.track_id(indtr)),...
    dataTipTextRow("partnerID",pidlabel)];
    hp.DataTipTemplate.DataTipRows(end+1:end+3) = dtRows;   
end


    axis(ax,'ij');
    axis(ax,'equal');

% Plot cotracks vs time
%%
%only both good
trackstat.coprogressive=(trackstat.color==1 & trackstat.comovement & trackstat.progressive & trackstat.progressivepartner);
goodpairs=find(trackstat.coprogressive);
figure(102)
clf
for k=1:length(goodpairs)
    subplot(5,5,k)
    hold off
    id1=locs.track_id==goodpairs(k);
    pid=trackstat.partnerids(goodpairs(k));
    id2=locs.track_id==pid;
    tmin=min(min(locs.frame(id1)),min(locs.frame(id2)));

    [x1,y1]=rotcoord(locs.xnm(id1)-mean(locs.xnm(id1)),locs.ynm(id1)-mean(locs.ynm(id1)),trackstat.angle(goodpairs(k)));
    [x2,y2]=rotcoord(locs.xnm(id2)-mean(locs.xnm(id2)),locs.ynm(id2)-mean(locs.ynm(id2)),trackstat.angle(pid));
    plot(locs.frame(id1),x1,'.-',locs.frame(id2),x2,'.-')
    hold on
    title(['frame: ' num2str(tmin) ', x: ' num2str(mean(locs.xpix(id1)),'%3.0f') ', y: ' num2str(mean(locs.ypix(id1)),'%3.0f')])  
end

% Calculate statistics
disp(sprintf('ch1\tch2\tch1prog\tch2prog\tdualcol\tch12progr'));
output=(sprintf([num2str(sum(trackstat.channel==1)), '\t' num2str(sum(trackstat.channel==2)),...
    '\t' num2str(sum(trackstat.channel==1 & trackstat.progressive)), ...
    '\t' num2str(sum(trackstat.channel==2 & trackstat.progressive)), ...
    '\t' num2str(sum(trackstat.partnertrackid>0)/2), ...
    '\t' num2str(sum(trackstat.coprogressive))]));
disp(output)
clipboard('copy',output)

%%
%add to movie
dx=3;
dimmark=3;
if markintiffile
    imstack=tiffreadVolume(tiffile);
    imstack=permute(imstack,[1 2 4 3]);
    immax=max(imstack(:));
    [lenx,leny]=size(imstack,[1,2]);
    idcoprogress=find(trackstat.coprogressive);
    for k=1:length(idcoprogress)
        ind=locs.track_id==idcoprogress(k);
        x=round(locs.xnm(ind)/pixelsize(1)-roi(1)); y=round(locs.ynm(ind)/pixelsize(2)-roi(2)); frame=locs.frame(ind);
        x=min(lenx-dx,x);x=max(dx+1,x);y=min(leny-dx,y);y=max(dx+1,y);
        for l=1:length(x)
            imstack(y(l)+dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
            imstack(y(l)-dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
            imstack(y(l)-dx:y(l)+dx, x(l)+dx,dimmark,frame(l))=immax;
            imstack(y(l)-dx:y(l)+dx, x(l)-dx,dimmark,frame(l))=immax;

        end

    end
fout2=strrep(tiffile,'_combined.tif','_mark.tif');
options.color=true;
saveastiff(squeeze(single(imstack)),fout2,options);
end

