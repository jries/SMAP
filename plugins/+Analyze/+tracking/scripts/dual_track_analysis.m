% analyze dual-color tracks

%filters
minlenframes=5;
maxd=350; %max distance to associate localizations

cotracklength=2; %minimum data poits associated between tracks
cotrackfraction = 0.0; %minimum fraction of localizations associated

% test directed movement
aspectratio=1; %should be below this for directed motion.
lennmstartend=150; % minimum endpoint- startpoint
lennmmin=200; %minimum of largest extenstion




obj=g;
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','frame','track_id','channel','track_length'},'layer',find(obj.getPar('sr_layerson')),'position','roi','grouping','ungrouped');

% test directed movement, calculate statistics
usetracks=unique(locs.track_id(locs.track_id>0));
trackstat.lennmstartend=zeros(max(usetracks),1);
trackstat.stdshort=zeros(max(usetracks),1);
trackstat.stdlong=zeros(max(usetracks),1);
trackstat.lenframe=zeros(max(usetracks),1);
trackstat.partnertrackid=zeros(max(usetracks),1);

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
    trackstat.stdshort(iduset)=real(sqrt(ev(1)));
    trackstat.stdlong(iduset)=real(sqrt(ev(2)));


               % [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
               % angle = cart2pol(xa, ya);          
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
figure(88)
hold off
ax=gca;

cols=[1 0 1
      0 1 1
      1 0 0 
      0 0 1
      .5 0 .5
      0 .5 .5];


for k=1:length(usetracks)
    idh=usetracks(k);
    validcotrack=false;
    validstats=true;
    indtr=locs.track_id==idh;
    colind=mode(locs.channel(indtr));
    
    lw=1;
    
    %look at stats
    validstats = validstats & trackstat.stdshort(idh)/trackstat.stdlong(idh)<aspectratio;
    validstats = validstats & trackstat.lennmstartend(idh) > lennmstartend;
    validstats = validstats & trackstat.stdlong(idh) > lennmmin;

    


    %look at co-movement
    pid=0;
    if trackstat.partnertrackid(idh)>0 %&& validstats
        partnerind=partner(indtr);
        partnerids=locs.track_id(partnerind(partnerind>0));
        [pid,npart]=mode(partnerids);
        if npart>cotracklength && npart/sum(indtr)>cotrackfraction
            validcotrack=true;
            lw=6;
            colind=colind+2;
            colind=colind+(~validstats)*2;
        end
    end

    

    % if any(pairedr(indtr)) 
    %     partnerind=partnerr(indtr);
    % elseif any(pairedt(indtr)) 
    %     partnerind=partnert(indtr);
    % end
    % 
    % if (any(pairedr(indtr)) || any(pairedt(indtr))) && validstats
    %     partnerids=locs.track_id(partnerind(partnerind>0));
    %     [pid,npart]=mode(partnerids);
    %     if npart>cotracklength && npart/sum(indtr)>cotrackfraction
    %         validcotrack=true;
    %         lw=6;
    %         colind=colind+2;
    %         colind=colind+(~validstats)*2;
    %     end
    % end
    
    hp=plot(ax,locs.xnm(indtr),locs.ynm(indtr),'.-','Color',cols(colind,:),'LineWidth',lw,'Tag','test');
    hold(ax,"on")
    % hp.DataTipTemplate.DataTipRows(1).Label = "X";
    % hp.DataTipTemplate.DataTipRows(2).Label = "Y"; 
    % hp.DataTipTemplate.DataTipRows(3).Label = "frame"; 
        pidlabel=0*locs.track_id(indtr)+pid;

        dtRows = [dataTipTextRow("frame",locs.frame(indtr)),...
        dataTipTextRow("ID",locs.track_id(indtr)),...
        dataTipTextRow("partnerID",pidlabel)];
        hp.DataTipTemplate.DataTipRows(end+1:end+3) = dtRows;
    
end




% Plot surviving tracks