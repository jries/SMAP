function results=make_lineprofiles(locD,p)

global jervistemp;
jervistemp=[];
ax1=initaxis(p.resultstabgroup,'scatter');
hold off
[ax2,tab2]=initaxis(p.resultstabgroup,'xprofile');
ax2.Position(3)=0.55;ax2.Position(1)=0.05;
hold off
[ax3,tab3]=initaxis(p.resultstabgroup,'yprofile');
ax3.Position(3)=0.55;ax3.Position(1)=0.05;
hold off

if isfield(locD.loc,'znm')
ax4=initaxis(p.resultstabgroup,'zprofile');
ax4.Position(3)=0.55;ax4.Position(1)=0.05;
hold off
ax5=initaxis(p.resultstabgroup,'scatter xz');

hold off
end


t1={[13 'x profile']};t2={[13 'y profile']};t3={[13 'z profile']};

if p.setbinwidth
    binwidth=p.binwidth;
else
    binwidth=p.sr_pixrec;
end

for layer=1:length(p.sr_layerson)
    if p.sr_layerson(layer)
        [locs,~, hroi]=locD.getloc({'xnm','ynm','znm','locprecnm','locprecznm','xnmline','ynmline'},'layer',layer,'position','roi');
    linew=p.linewidth_roi/2;
    if isa(hroi,'imline')
        x=locs.xnmline;
        y=locs.ynmline;
        mx=0;
        my=0;
        pos=getPosition(hroi);
        linel=sqrt((pos(2,1)-pos(1,1))^2+(pos(2,2)-pos(1,2))^2)/2*1000;

    elseif isa(hroi,'impoint')
        x=locs.xnm;
        y=locs.ynm;
        pos=getPosition(hroi)*1000;
        mx=pos(1);
        my=pos(2);
        linel=linew;
    else
        x=locs.xnm;
        y=locs.ynm;
        pos=getPosition(hroi);
        mx=mean(pos(:,1));
        my=mean(pos(:,2));
    end
   
    z=double(locs.znm);
        locprecnm=double(locs.locprecnm);
    locprecznm=double(locs.locprecznm);
    if isempty(z)
        z=x*0;
        locprecznm=z;
    end

    
    if p.linelengthcheck

        indin=x>mx-p.linelength/2&x<mx+p.linelength/2;
        x=x(indin);y=y(indin);z=z(indin);locprecnm=locprecnm(indin);%locprecznm=locprecznm(indin);
        linel=p.linelength/2;
    end
    
    axes(ax1)
    plot(x,y,'.')
    axis equal tight
    hold on
    
%     t1{end+1}=['Layer ' 9 int2str(layer)];
%     t2{end+1}=['Layer ' 9 int2str(layer)];
%     t3{end+1}=['Layer ' 9 int2str(layer)];
    t1{end+1}=p.layernames{layer};
    t2{end+1}=p.layernames{layer};
    t3{end+1}=p.layernames{layer};  
    
    n=my-linew-binwidth:binwidth:my+linew+binwidth;
    profy=hist(y,n);profy([1 end])=[];n([1 end])=[];
    jervistemp.ny=n;
    jervistemp.profy=profy;
    axes(ax2)
    plot(n,profy);
    hold on
    fwhm=getFWHM(profy,n);
    t1{end+1}=['FWHM: ' 9 num2str(fwhm,3)];
    
    sigma=median(locprecnm);
    [fitp,fitprof,fittxt]=fitgeneral(profy,n,p,sigma);
    plot(n,fitprof,'k--')
    t1(end+1:end+length(fittxt))=fittxt;
    

    n=mx-linel-binwidth:binwidth:mx+linel+binwidth;
    jervistemp.binwidth=binwidth;
    profx=hist(x,n);profx([1 end])=[];n([1 end])=[];
    jervistemp.profx=profx;
    jervistemp.nx=n;
    axes(ax3)
    plot(n,profx);
    hold on
    
    fwhm=getFWHM(profx,n);
    t2{end+1}=['FWHM: ' 9 num2str(fwhm)];
     [fitp,fitprof,fittxt]=fitgeneral(profx,n,p,sigma);
    plot(n,fitprof,'k--')
    t2(end+1:end+length(fittxt))=fittxt;
    
    if ~isempty(locs.znm)
        minzh=max(-750,min(z));
        maxzh=min(750,max(z));
    n=minzh-3*binwidth:binwidth:maxzh+3*binwidth;
    profz=hist(z,n);profz([1 end])=[];n([1 end])=[];
    axes(ax4)
    plot(n,profz);
    hold on
    fwhm=getFWHM(profz,n);
    t3{end+1}=['FWHM: ' 9  num2str(fwhm)];
    [fitp,fitprof,fittxt]=fitgeneral(profz,n,p,fwhm/2.6);
    plot(n,fitprof,'k--')
    t3(end+1:end+length(fittxt))=fittxt;
    
    jervistemp.profz=profz;jervistemp.nz=n;
    
    
    axes(ax5)
    plot(x,z,'.')
    axis equal tight
    hold on
    end
    t1{end+1}='';
    t2{end+1}='';
    t3{end+1}='';
    end
end

fontsize=15;
pos=[.65,0.025,.4,.95];
uicontrol('Parent',ax2.Parent,'style','text','String',t1,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')

uicontrol('Parent',ax3.Parent,'style','text','String',t2,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')

if ~isempty(locs.znm)
    uicontrol('Parent',ax4.Parent,'style','text','String',t3,'Units','normalized','Position',pos,'FontSize',fontsize,'HorizontalAlignment','left')

end
tab2.Parent.SelectedTab=tab2;
assignin('base','jervistemp',jervistemp);
results=[t1 t2 t3];

function [fwhm,fwhmind]=getFWHM(profile,x)
    [mp, ip]=max(profile);
    i1=find(profile(1:ip)>mp/2,1,'first');
    i2=find(profile(ip:end)>mp/2,1,'last')+ip-1;
    if isempty(i2)||isempty(i1)
        fwhm=[];
        fwhmind=1;
    else
        if i1==i2
            i1=i1-1;i2=i2+1;
        end
        
    fwhm=x(i2)-x(i1);
    fwhmind=i2-i1;
    end


function [fitp,fitprof,fittext]=fitgeneral(profile,x,p,sigma)
if p.restrictsigma
    fif=@convoluteshape;
else
    fif=@convoluteshapefits;
end
switch p.fitmodel.Value
    case 1%Gauss
        [fitp,fitprof,fittext]=fitgauss(profile,x);
        fitp(5)=fitp(3);
    case 2%tophat
        [fitp,fitprof,fstart]=fif(x,profile,sigma,1);
        fittext=['Step L: ' 9 num2str(fitp(3),3) ];
    case 4%ring
        [fitp,fitprof,fstart]=fif(x,profile,sigma,3);
        fittext=['Ring R: ' 9 num2str(fitp(3),3) ];
    case 3%disk
        [fitp,fitprof,fstart]=fif(x,profile,sigma,2);
        fittext=['Disk R: ' 9 num2str(fitp(3),3)];   
    case 5 %double Gauss distance
%         [fitp,fitprof,fstart]=fif(x,profile,sigma,4);
        [fitp,fitprof,fittext]=fit2gauss(profile,x);
%         fittext=['Distance d: ' 9 num2str(fitp(4),3)];         

end
fittext={fittext};
if p.fitmodel.Value>1
if p.restrictsigma
    fittext(end+1)={['sigma: ' 9 num2str(sigma,4)]};
else
    fittext(end+1)={['sigma: ' 9 num2str(fitp(5),4)]}; 
end
end
    
function [fitp,fitprof,fittext]=fitgauss(profile,x)
[~,s]=getFWHM(profile,x);
s=s/2.6*(x(2)-x(1));
[mp, ip]=max(profile);
startp=[mp x(ip) s 0];
fitp=mygaussfit(x,profile,startp);
fitprof=mygaussforfit(fitp,x);
fittext=['sigma: ' 9 num2str(fitp(3),4)];

function [fitp,fitprof,fittext]=fit2gauss(profile,x)
%try with two Gaussian fits.
fp1=fit(x',profile','gauss1','Lower',[0 -inf 0],'Robust','Bisquare');
fp2=fit(x',profile'-(fp1(x)),'gauss1','Lower',[0 -inf 0],'Robust','Bisquare');

ft=fittype('a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c1)^2)+d');
% fp=fit(x',profile','gauss2','Lower',[0 -inf 0 0 -inf 0]);
% fitp=[fp.a1 fp.b1 fp.c1 fp.a2 fp.b2 fp.c2];
startp=[fp1.a1 fp2.a1 fp1.b1 fp2.b1 fp2.c1 0];

fp=fit(x',profile',ft,'Lower',[0 0 -inf -inf 0 0],'StartPoint',startp);
% [~,s]=getFWHM(profile,x);
% s=s/2.6*(x(2)-x(1));
% [mp, ip]=max(profile);
% startp=[mp x(ip) s 0];
% fitp=mygaussfit(x,profile,startp);
fitp=[fp.a1 fp.a2 fp.b1 fp.b2-fp.b1 fp.c1 fp.d ];
fitprof=fp(x);
fittext=['Distance d: ' 9 num2str(abs(fp.b2-fp.b1),4)];
