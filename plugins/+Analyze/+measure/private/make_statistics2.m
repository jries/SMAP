function stat=make_statistics2(locs,p,ploton)
if nargin<3
    ploton=true;
end

if p.filter
    modetxt=p.layernames(logical(p.sr_layerson));
%     modetxt={'layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer'};
else
    modetxt={'ungroup','group'};
end

if isfield(locs{1},'znm')&&~isempty(locs{1}.znm)
    txt='znm';
   zexist=true;
else
    txt='PSFxnm';
        zexist=false;
end

% fn=locD.files.file.name;
if ploton
if p.overview
    figure(34);


    sf=0.75;
    ax1=subplot(3,2,1);ax1.Position(3)=ax1.Position(3)*sf;
   
    ax2=subplot(3,2,2);ax2.Position(3)=ax2.Position(3)*sf;
    ax3=subplot(3,2,3);ax3.Position(3)=ax3.Position(3)*sf;
    ax4=subplot(3,2,4);ax4.Position(3)=ax4.Position(3)*sf;
    ax5=subplot(3,2,5);ax5.Position(3)=ax5.Position(3)*sf;
   
    ax6=subplot(3,2,6);ax6.Position(3)=ax6.Position(3)*sf;ax7=[];ax8=[];
else
    ax1=initaxis(p.resultstabgroup,'Photons');
    ax1.Position(3)=0.55;
    ax2=initaxis(p.resultstabgroup,'locprec');
    ax2.Position(3)=0.55;
    ax3=initaxis(p.resultstabgroup,'lifetime');
    ax3.Position(3)=0.55;
    ax4=initaxis(p.resultstabgroup,'BG');
    ax4.Position(3)=0.55;
    ax5=initaxis(p.resultstabgroup,txt);
    ax5.Position(3)=0.55;
    ax6=initaxis(p.resultstabgroup,'frames');
    ax6.Position(3)=0.55;
    if zexist
        ax7=initaxis(p.resultstabgroup,['err ' txt]);
        ax7.Position(3)=0.55;
        ax8=initaxis(p.resultstabgroup,['err(' txt ')']);
        ax8.Position(3)=0.55;
    end


end
else
    ax1=[];ax2=[];ax3=[];ax4=[];ax5=[];ax6=[];ax7=[];ax8=[];
end
datrange=1:length(locs);

slegend={};
% look at frames
for k=datrange
    slegend{k}=[modetxt{k} num2str(datrange(k))];
%     slegend{end+1}='';
    frames=locs{k}.frame;
    if isempty(frames)
        continue
    end
    mf=max(frames);
    [hfr,n]=hist(frames,10);
    hfrc=hfr;
    hfrc(1:2)=[]; %ignore beginning
%     [~,ind]=max(diff(1./(hfr+.1*mean(hfr))));
%     ind=find(hfrc<(max(hfrc)-min(hfrc))/3+min(hfrc),1,'first');
    ind=find(hfrc<(max(hfrc))/3,1,'first');
    if isempty(ind)
        ind=length(n);
    else
    ind=ind+2;
    end
    mr=0.3;
    r2=round([n(ind)-mf*mr n(ind)+mf*mr]);
    frames2=frames(frames>=r2(1)&frames<=r2(2));
    [hfr2,n2]=hist(frames2,20);
    indco2=find(hfr2>(max(hfr2)-min(hfr2))/2+min(hfr2),1,'last');
    if isempty(indco2)
        indco2=length(n2);
    end
    falloffframe=n2(indco2);
    stat.frames.falloff(k)=falloffframe;
    [stat.frames.histogram(k).h,stat.frames.histogram(k).n]=hist(frames,100);
end

if ploton
%     axf=initaxis(p.resultstabgroup,'frames');
    axf=ax6;
    hold(axf, 'off')
    for k=datrange
        plothf(k)=plot(axf,stat.frames.histogram(k).n,stat.frames.histogram(k).h);
           hold(axf, 'on')
        plot(axf,ones(2,1)*stat.frames.falloff(k),[0,max(stat.frames.histogram(k).h)])
    end
    legend(plothf,slegend);
end
slf={'Frames'};

%photon stats
phot=getFieldAsVector(locs,'phot');
if isempty(phot{1})
    errdlg('no localizations in selected region')
    error('no localizations in selected region')
end
% if p.checkphot
%     for k=datrange
%         phot{k}(phot{k}<p.photrange(1))=[];
%         if length(p.photrange)>1
%              phot{k}(phot{k}>p.photrange(2))=[];
%         end
%     end
%     pr=p.photrange;
% else
    pr=0.99;
% end
[hphot,mmax]=plothist(phot,pr,[],0,ax1,modetxt,40);
sphot={'Photons'};
% phot1=1000;
% phot2=3000;
for k=datrange
    sphot{end+1}='';
    sphot{end+1}=[num2str(k) '.' modetxt{k} ];
    Nloc(k)=length(phot{k});
    meanphot(k)=mean(phot{k});
%     N1(k)=sum(phot{k}>phot1);
%     N2(k)=sum(phot{k}>phot2);
    inrange=phot{k}>p.photrange(1)&phot{k}<p.photrange(2);
    meanphotrange(k)=mean(phot{k}(inrange));
    
    sphot{end+1}=['N'  ' = ' num2str(Nloc(k)/1000,'%5.0f') 'k'];
    sphot{end+1}=['<P_all'  '> = ' num2str(meanphot(k),'%5.0f')];
    sphot{end+1}=['<P_range'  '> = ' num2str(meanphotrange(k),'%5.0f')];
%     sphot{end+1}=['r'  ' = ' num2str(N1(k)/N2(k),'%5.2f')];
%     dat(k)=fitexpphot(hphot{k},[],ploton);
%     dat(k)=meanexphere(phot{k}(inrange),hphot{k},p.photrange,ax1,mmax{k});
    dat(k)=meanexphere(phot{k},hphot{k},p.photrange,ax1,mmax{k});
    
    sphot{end+1}=(['Pexp'  ' = ' num2str(dat(k).mu,'%5.0f')]);   
end
stat.photons.Nloc=Nloc;
stat.photons.meanphot=meanphot;
stat.photons.mu=[dat(:).mu];

%locprec
locp=getFieldAsVector(locs,'locprecnm');
hlocp=plothist(locp,0.99,.25,0,ax2,modetxt);
slp={'locprec_x'};
for k=datrange
    slp{end+1}='';
    slp{end+1}=[num2str(k) '.' modetxt{k} ];
    loch=locp{k};
    loch(loch<=0)=[];
    loch(loch>10000)=[];
    px = mylognfit(loch);
    [~,ind]=max(hlocp{k}.h);
    smx=hlocp{k}.n(ind);
    
    %refine
    dwin=ceil(ind/4);
    rn=max(1,ind-dwin+1):min(ind+dwin,length(hlocp{k}.h));
    fpol=fit(hlocp{k}.n(rn)',hlocp{k}.h(rn)','poly3');
    smxf=fzero(@(x) 3*fpol.p1*x.^2+2*fpol.p2*x+fpol.p3,smx);
    if ~isempty(ax2)
        plot(hlocp{k}.n(rn),fpol(hlocp{k}.n(rn)),'g')
    end
    
    stat.locprec.max(k)=smxf;
    slp{end+1}=['max: ' num2str(smxf,3)];
%     hf=mylognpdf(hlocp{k}.n,px(1),px(2))*sum(hlocp{k}.h)*(hlocp{k}.n(2)-hlocp{k}.n(1));
%     if ploton
%     plot(hlocp{k}.n,hf/max(hlocp{k}.h),'k:')
%     end
    slp{end+1}=['median: ' num2str(median(locp{k}),3)];
    stat.locprec.median(k)=median(locp{k});
    %risng edge
    
    
    [~,indrise1]=find(hlocp{k}.h>0.2,1,'first');
    indrise1=min(4*indrise1,length(hlocp{k}.h));
    imaxx=max(hlocp{k}.h(1:indrise1));
    indrise=find(hlocp{k}.h>imaxx/2,1,'first');
    risingedge=hlocp{k}.n(indrise);
    stat.locprec.rising(k)=risingedge;
    slp{end+1}=['rising: ' num2str(risingedge,3)];
    
    geom=geomean(loch);
    slp{end+1}=['geomean: ' num2str(geom,3)];
end

%lifetime
lifetime=getFieldAsVector(locs,'numberInGroup');
% lifetimeall=lifetime;
% if isfield(p,'checklifetime')&&p.checklifetime
%     for k=datrange
%         lifetime{k}(lifetime{k}<p.lifetimerange(1))=[];
%         if length(p.lifetimerange)>1
%              lifetime{k}(lifetime{k}>p.lifetimerange(2))=[];
%         end
%     end
%     plr=p.lifetimerange;
% else
    plr=0.995;
% end


[hlifet,mmax]=plothist(lifetime,plr,1,0,ax3,modetxt);
 ax3.NextPlot='add';
slt={'lifetime'};
for k=datrange
    inrange=lifetime{k}>p.lifetimerange(1)&lifetime{k}<p.lifetimerange(2);
    slt{end+1}='';
    slt{end+1}=[num2str(k) '.' modetxt{k} ];
%     dat(k)=fitexpphot(hlifet{k},2,ploton);
    dat(k)=meanexphere(lifetime{k},hlifet{k},p.lifetimerange,ax3,mmax{k});
    slt{end+1}=(['texp'  ' = ' num2str(dat(k).mu,3)]);
    slt{end+1}=(['meanrange'  ' = ' num2str(mean(lifetime{k}(inrange)),3)]);
    slt{end+1}=(['meanall'  ' = ' num2str(mean(lifetime{k}),3)]);
    stat.lifetime.mu(k)=dat(k).mu;
end

%background
bg=getFieldAsVector(locs,'bg');
slb={'Background'};
if isempty(bg{1})
bg=getFieldAsVector(locs,'bg2');   
slb={'Background2'};
end
hbg=plothist(bg,0.95,1,0,ax4,modetxt);

for k=datrange
    slb{end+1}='';
    slb{end+1}=[num2str(k) '.' modetxt{k} ];
    mbg=mean(bg{k});
    slb{end+1}=['mean: ' num2str(mbg,'%5.0f')];
    [~,mind]=max(hbg{k}.h);
    
    maxbg=hbg{k}.n(mind);
        mx2=max(hbg{datrange(k)}.h(1:ceil(mind*0.55)));
    if mx2>0.3*maxbg %if by grouping second max is higher select the first max.
        maxbg=mx2;
    end
    
    slb{end+1}=['max: ' num2str(maxbg,'%5.0f')];
    stat.background.mean(k)=mbg;
    stat.background.max(k)=maxbg;
end

%z/sigma
if zexist
    v=getFieldAsVector(locs,'znm');
else
    v=getFieldAsVector(locs,'PSFxnm');
end
hz=plothist(v,.99,[],[],ax5,modetxt);
sls={txt};
for k=1:length(datrange)
    sls{end+1}='';
    sls{end+1}=[num2str(datrange(k)) '.' modetxt{datrange(k)} ];
    if length(hz)<k || isempty(hz{datrange(k)})
        continue
    end
    [~,ind]=max(hz{datrange(k)}.h);
    mx=hz{datrange(k)}.n(ind);   



    
    %refine
%     dwin=ceil(ind/4);
    dwin=8;
    rn=max(1,ind-dwin+1):min(ind+dwin,length(hz{datrange(k)}.h));
    fpol=fit(hz{datrange(k)}.n(rn)',hz{datrange(k)}.h(rn)','poly3');
    smxf=fzero(@(x) 3*fpol.p1*x.^2+2*fpol.p2*x+fpol.p3,mx);
    if ~isempty(ax2)
        plot(hz{datrange(k)}.n(rn),fpol(hz{datrange(k)}.n(rn)),'g')
    end
    
    
    sls{end+1}=['max: ' num2str(smxf,4)];
    stat.(txt).max(k)=smxf;
end


if ploton
fontsize=14;
if p.overview
    pos=[.25,0.0,-.1,0];
    fontsize=10;
    uicontrol('Parent',ax1.Parent,'style','text','String',sphot,'Units','normalized','Position',ax1.Position+pos,'FontSize',fontsize,'HorizontalAlignment','left')
    uicontrol('Parent',ax2.Parent,'style','text','String',slp,'Units','normalized','Position',ax2.Position+pos,'FontSize',fontsize,'HorizontalAlignment','left')
    uicontrol('Parent',ax3.Parent,'style','text','String',slt,'Units','normalized','Position',ax3.Position+pos,'FontSize',fontsize,'HorizontalAlignment','left')
    uicontrol('Parent',ax4.Parent,'style','text','String',slb,'Units','normalized','Position',ax4.Position+pos,'FontSize',fontsize,'HorizontalAlignment','left')
    uicontrol('Parent',ax5.Parent,'style','text','String',sls,'Units','normalized','Position',ax5.Position+pos,'FontSize',fontsize,'HorizontalAlignment','left')
    uicontrol('Parent',ax6.Parent,'style','text','String',slf,'Units','normalized','Position',ax6.Position+pos,'FontSize',fontsize,'HorizontalAlignment','left')
else
pos=[.7,0.025,.3,.95];

uicontrol('Parent',ax1.Parent,'style','text','String',sphot,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax2.Parent,'style','text','String',slp,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax3.Parent,'style','text','String',slt,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax4.Parent,'style','text','String',slb,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax5.Parent,'style','text','String',sls,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
uicontrol('Parent',ax6.Parent,'style','text','String',slf,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
end
end

if zexist
    v=getFieldAsVector(locs,'locprecznm');
    hz=plothist(v,.99,[],0,ax7,modetxt);
    slp={'locprecznm'};
    for k=datrange
        slp{end+1}='';
        slp{end+1}=[num2str(k) '.' modetxt{k} ];
        [~,ind]=max(hz{k}.h);
        mx=hz{k}.n(ind);    
        slp{end+1}=['max: ' num2str(mx,3)];
        stat.locprecznm.max(k)=mx;
    end    
    znm=getFieldAsVector(locs,'znm');
    rz=[-800 800];
    rsz=[0 100];
    him=myhist2(znm{1},v{1},10,1,rz,rsz);    
    if ploton &&~isempty(ax8)
        axes(ax8)
        imagesc(rz,rsz,him')
        axis xy
        xlabel('znm');ylabel('locprec z');
    end
    if ploton
        if p.overview
            uicontrol('Parent',ax6.Parent,'style','text','String',slp,'Units','normalized','Position',ax6.Position+pos,'FontSize',fontsize,'HorizontalAlignment','left')
        else
        uicontrol('Parent',ax6.Parent,'style','text','String',slp,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')
        end
    end   
end

if ploton && ~p.overview
   ax1.Parent.Parent.SelectedTab=ax1.Parent;
end

function [v,datrange]=getvals(locD,field,p,indin)
if p.filter %use filtered values
    for layer=1:length(p.sr_layerson)
        if p.sr_layerson(layer)
            if p.useroi                
                v{layer}=locD.getloc(field,'layer',layer,'position','roi','within',indin).(field);
            else
                v{layer}=locD.getloc(field,'layer',layer,'within',indin).(field);
            end
            
        else
            v{layer}=0;
        end
    end
    datrange=find(p.sr_layerson);
else %use all values, plot for unconnected and connected
    if p.useroi
        position='roi';
    else
        position='all';
    end
       struc=locD.getloc(field,'position',position,'grouping','ungrouped','within',indin);
       v{1}=struc.(field);
       struc=locD.getloc(field,'position',position,'grouping','grouped','within',indin);
       v{2}=struc.(field);
    datrange=1:2;
end

function [his,mmo]=plothist(v,quantile,dphot,hmin,ax,modetxt,qfac)
if nargin<7
    qfac=5;
end
his=[];
for k=1:length(v)
    if length(quantile)==1
    qq=myquantilefast(v{k},[1-quantile,quantile],30/(1-quantile));
    else
    qq=quantile;
    end
    q(k)=qq(2);q0(k)=qq(1);

    l(k)=length(v{k});
end

qmax=(max(q));
qmin=min(q0);
if qmax==qmin
    qmax=qmin+1;
end
lmax=max(l);
% qfac=log10(lmax)-1-1;
% qfac=5;

if nargin==2||isempty(dphot)
dphot=(10^ceil(log10(qmax/qfac)))/100;
end
if nargin<4||isempty(hmin)
    hmin=qmin;
end
nphot=hmin:dphot:qmax;
slegend={};
if ~isempty(ax)
    axes(ax)
    hold off
end
for k=1:length(v)
    
    if q(k)>0
%         sum(v{k})
        h=hist(v{k},nphot);
        [mmax,mi]=max(h(2:end-1)); 
        his{k}.h=h(2:end-1)/mmax;
        his{k}.n=nphot(2:end-1);
        if ~isempty(ax)
        plot(nphot(2:end-1),h(2:end-1)/mmax)
        hold on
        end
        slegend{end+1}=[modetxt{k}];
        mmo{k}=mmax;
    end
end
if ~isempty(ax)
legend(slegend,'Location','northeast')
end

function dat=fitexpphot(hin,fitstart,ploton)
h=double(hin.h);
xout=double(hin.n);
if length(h)>1
    [mmax,mi]=max(h(1:end-1)); 
    halft=find(h(mi:end)<mmax/2,1,'first')+mi;
if isempty(halft)
    halft=ceil(length(h)/2);
end
if nargin<2||isempty(fitstart)
    fitstart=ceil(mi*1.2);
end
fitr=fitstart:min(halft*5,length(h));

options=optimset('lsqcurvefit');
options.Display='off';
pf=lsqcurvefit(@expforfit,[1,xout(halft)],xout(fitr),h(fitr)/mmax,[],[],options);
if ploton
    plot(xout(fitr),expforfit(pf,xout(fitr)),'k--')
end
dat.mu=pf(2);
else
    dat.mu=0;  
end

function dat=meanexphere(v,hin,fitrange,ax,fac)
try
    h=double(hin.h);
   [mmax,mi]=max(h(1:end-1)); 
   xout=double(hin.n);
   maxpos=xout(mi);
if nargin<2||isempty(fitrange)
    
    
    if length(h)>1
        
        halft=find(h(mi:end)<mmax/2,1,'first')+mi;
    if isempty(halft)
        halft=ceil(length(h)/2);
    end
    fitstartind=ceil(mi*1.0);
    fitrange(1)=hin.n(fitstartind);
    fitrange(2)=hin.n(end);
    end
end
if length(fitrange)<2
    fitrange(2)=myquantilefast(v,.99,10000);
end
fitrange(1)=max(fitrange(1),maxpos);
% fitr=fitstart:min(length(h));
% % fitr=fitstart:min(halft*5,length(h));
dq=hin.n(2)-hin.n(1);
% rangev=[hin.n(fitr(1)) hin.n(fitr(end))];
% if ploton
%     ax=gca;
% else
%     ax=[];
% end
ax.NextPlot='add';
dat.mu=meanexp(v,dq,fitrange,ax,fac);
xlim(ax,[hin.n(1) hin.n(end)])
catch
    dat.mu=0;  
end

function out=expforfit(p,x)
        out=p(1)*exp(-x/p(2));

