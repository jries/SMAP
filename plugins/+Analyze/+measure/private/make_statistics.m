function stat=make_statistics(locD,p,indin,ploton)
if nargin<3||isempty(indin)
    indin=[];
end
if nargin<4
    ploton=true;
end
% ploton=false
if p.filter
    modetxt={'layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer','layer'};
else
    modetxt={'ungroup','group'};
end

if isfield(locD.loc,'znm')
    txt='znm';
   zexist=true;
else
    txt='PSFx';
        zexist=false;
end

fn=locD.files.file.name;
if ploton
if p.overview
    figure(34);


    sf=0.75;
    ax1=subplot(3,2,1);ax1.Position(3)=ax1.Position(3)*sf;
   
    ax2=subplot(3,2,2);ax2.Position(3)=ax2.Position(3)*sf;
    ax3=subplot(3,2,3);ax3.Position(3)=ax3.Position(3)*sf;
    ax4=subplot(3,2,4);ax4.Position(3)=ax4.Position(3)*sf;
    ax5=subplot(3,2,5);ax5.Position(3)=ax5.Position(3)*sf;
    ax6=subplot(3,2,6);ax6.Position(3)=ax6.Position(3)*sf;
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
    ax6=initaxis(p.resultstabgroup,['error ' txt]);
    ax6.Position(3)=0.55;
    if zexist
    ax7=initaxis(p.resultstabgroup,['error ' txt ' vs ' txt]);
    ax7.Position(3)=0.55;
    end


end
else
    ax1=[];ax2=[];ax3=[];ax4=[];ax5=[];ax6=[];ax7=[];
end
%photon stats
[phot,datrange]=getvals(locD,'phot',p,indin);
if isempty(phot)
    errdlg('no localizations in selected region')
    error('no localizations in selected region')
end

if p.checkphot
    for k=datrange
        phot{k}(phot{k}<p.photrange(1))=[];
        if length(p.photrange)>1
             phot{k}(phot{k}>p.photrange(2))=[];
        end
    end
end

% if ploton
% axes(ax1)
% hold off
% end
hphot=plothist(phot,0.95,[],0,ax1);
sphot={'Photons'};

phot1=1000;
phot2=3000;

for k=datrange
    sphot{end+1}='';
    sphot{end+1}=[num2str(k) '.' modetxt{k} ];
    Nloc(k)=length(phot{k});
    meanphot(k)=mean(phot{k});
    N1(k)=sum(phot{k}>phot1);
    N2(k)=sum(phot{k}>phot2);
    
    sphot{end+1}=['N'  ' = ' num2str(Nloc(k)/1000,'%5.0f') 'k'];
    sphot{end+1}=['<P'  '> = ' num2str(meanphot(k),'%5.0f')];
    sphot{end+1}=['r'  ' = ' num2str(N1(k)/N2(k),'%5.2f')];
    dat(k)=fitexpphot(hphot{k});
    sphot{end+1}=['\mu'  ' = ' num2str(dat(k).mu,'%5.0f')];
    
end

stat.photons.Nloc=Nloc;
stat.photons.meanphot=meanphot;
stat.photons.mu=[dat(:).mu];


%locprec
[locp,datrange]=getvals(locD,'locprecnm',p,indin);

hlocp=plothist(locp,0.99,.25,0,ax2);
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
    stat.locprec.max(k)=smx;
    slp{end+1}=['max: ' num2str(smx,3)];
    hf=mylognpdf(hlocp{k}.n,px(1),px(2))*sum(hlocp{k}.h)*(hlocp{k}.n(2)-hlocp{k}.n(1));
    plot(hlocp{k}.n,hf/max(hlocp{k}.h),'k:')
    slp{end+1}=['median: ' num2str(median(locp{k}),3)];
    stat.locprec.median(k)=median(locp{k});
    indrise=find(hlocp{k}.h>1/2,1,'first');
    risingedge=hlocp{k}.n(indrise);
    stat.locprec.rising(k)=risingedge;
    slp{end+1}=['rising: ' num2str(risingedge,3)];
end



%lifetime
[lifetime,datrange]=getvals(locD,'numberInGroup',p,indin);
% if ploton
% axes(ax3)
% % axp=recgui.initaxis(p.resultstabs(3),'lifet');
% % axp.Position(3)=0.55;
% hold off
% end
hlifet=plothist(lifetime,0.999,1,0,ax3);
slp={'lifetime'};
for k=datrange
    slp{end+1}='';
    slp{end+1}=[num2str(k) '.' modetxt{k} ];
    dat(k)=fitexpphot(hlifet{k},2);
    slp{end+1}=['\mu'  ' = ' num2str(dat(k).mu,3)];
    stat.lifetime.mu(k)=dat(k).mu;
end


%background
[bg,datrange]=getvals(locD,'bg',p,indin);
% axp=recgui.initaxis(p.resultstabs(4),'BG');
% axp.Position(3)=0.55;
% if ploton
% axes(ax4)
% hold off
% end
hbg=plothist(bg,0.95,1,0,ax4);
slp={'Background'};
for k=datrange
    slp{end+1}='';
    slp{end+1}=[num2str(k) '.' modetxt{k} ];
    mbg=mean(bg{k});
    slp{end+1}=['BG: ' num2str(mbg,'%5.0f')];
    stat.background.mean(k)=mbg;
end



%z/sigma
if isfield(locD.loc,'znm')
    [v,datrange]=getvals(locD,'znm',p,indin);
    txt='znm';

else
    [v,datrange]=getvals(locD,'PSFxnm',p,indin);
    txt='PSFx';
    zexist=false;
end



hz=plothist(v,.99,[],0,ax5);
slp={txt};
for k=1:length(datrange)
    slp{end+1}='';
    slp{end+1}=[num2str(datrange(k)) '.' modetxt{datrange(k)} ];
    [~,ind]=max(hz{datrange(k)}.h);
    mx=hz{datrange(k)}.n(ind);    
    slp{end+1}=['max: ' num2str(mx,3)];
end
 if ploton
text(ax1.XLim(2)*1.1,0.5,sphot,'FontSize',14,'Parent',ax1)
text(0,-0.2,fn,'FontSize',10,'Interpreter','none','Parent',ax1)

text(ax2.XLim(2)*1.1,0.5,slp,'FontSize',14,'Parent',ax2)
text(ax3.XLim(2)*1.1,0.5,slp,'FontSize',14,'Parent',ax3)

text(ax4.XLim(2)*1.1,0.5,slp,'FontSize',14)

text(ax5.XLim(2)*1.1,0.5,slp,'FontSize',14)
end

if zexist
    [v,datrange]=getvals(locD,'locprecznm',p,indin);
%     axes(ax6)
%     hold off
    hz=plothist(v,.99,[],0,ax6);
    slp={txt};
    for k=datrange
        slp{end+1}='';
        slp{end+1}=[num2str(k) '.' modetxt{k} ];
        [~,ind]=max(hz{k}.h);
        mx=hz{k}.n(ind);    
        slp{end+1}=['max: ' num2str(mx,3)];
    end
    
    
    
   
%     
    znm=getvals(locD,'znm',p,indin);
%     hold off
rz=[-800 800];
rsz=[0 100];
    him=myhist2(znm{1},v{1},10,1,rz,rsz);
    
    if ploton
    axes(ax7)
    imagesc(rz,rsz,him')
    axis xy
    xlabel('znm');ylabel('locprec z');
    end
%     for k=datrange
%         [~, znmm, znms]=bintrace(sort(znm{k},length(znm{k})/50);
%         [~, vm, vs]=bintrace(v{k},length(v{k})/50);
%         
%         plot(znmm,vm,'.')
        
%         hold on
%     end
if ploton
    text(ax6.XLim(2)*1.1,0.5,slp,'FontSize',14)
end
    
end

    


function [v,datrange]=getvals(locD,field,p,indin)
if p.filter %use filtered values
    for layer=1:length(p.sr_layerson)
        if p.sr_layerson(layer)
            if p.useroi                
                v{layer}=locD.getloc(field,'layer',layer,'position','roi','within',indin).(field);
%                 v{layer}=struc.(field);
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



function his=plothist(v,quantile,dphot,hmin,ax)

for k=1:length(v)
    qq=myquantile(v{k},[quantile,1-quantile]);
    q(k)=qq(1);q0(k)=qq(2);
    l(k)=length(v{k});
end
qmax=(max(q));
qmin=min(q0);
if qmax==qmin
    qmax=qmin+1;
end
lmax=max(l);

qfac=log10(lmax)-1.5;

if nargin==2||isempty(dphot)
dphot=(10^ceil(log10(qmax/qfac)))/100;
end
if nargin<4
    hmin=qmin;
end
nphot=hmin:dphot:qmax;
% if length(nphot)==1
%     nphot=[nphot; 2*nphot];
% end
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
        slegend{end+1}=num2str(k);
    end
end
if ~isempty(ax)
legend(slegend,'Location','northeast')
end

function dat=fitexpphot(hin,fitstart)
h=double(hin.h);
xout=double(hin.n);
if length(h)>1
[mmax,mi]=max(h(1:end-1)); 
halft=find(h(mi:end)<mmax/2,1,'first')+mi;
if isempty(halft)
 halft=ceil(length(h)/2);
end
if nargin<2
    fitstart=ceil(mi*1.2);
end
fitr=fitstart:min(halft*5,length(h));

options=optimset('lsqcurvefit');
options.Display='off';
pf=lsqcurvefit(@expforfit,[1,xout(halft)],xout(fitr),h(fitr)/mmax,[],[],options);
plot(xout(fitr),expforfit(pf,xout(fitr)),'k--')
dat.mu=pf(2);
else
    dat.mu=0;
    
end

function out=expforfit(p,x)
        out=p(1)*exp(-x/p(2));

