function results=make_lineprofiles(locD,p)

global jervistemp;
jervistemp=[];
ax1=initaxis(p.resultstabgroup,'scatter');
hold off
[ax2,tab2]=initaxis(p.resultstabgroup,'x profile');
ax2.Position(3)=0.55;ax2.Position(1)=0.08;
hold off
[ax3,tab3]=initaxis(p.resultstabgroup,'y profile');
ax3.Position(3)=0.55;ax3.Position(1)=0.08;
hold off

if isfield(locD.loc,'znm')
ax4=initaxis(p.resultstabgroup,'z profile');
ax4.Position(3)=0.55;ax4.Position(1)=0.08;
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
    plot(ax1,x,y,'.')
    axis(ax1,'equal')
     axis(ax1,'tight')
    xlabel(ax1,'Position along line ROI (nm)')
    ylabel(ax1,'Position perpendicular to line ROI (nm)')
    hold(ax1,'on')
    
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
    plot(ax2,n,profy);
    hold(ax2,'on')
    
    xlabel(ax2,'Position perpendicular to line ROI (nm)')
    ylabel(ax2,'counts')
    
    fwhm=getFWHM(profy,n);
    t1{end+1}=['FWHM: ' 9 num2str(fwhm,3)];
    
    sigma=median(locprecnm);
    [fitp,fitprof,fittxt]=fitgeneralprofile(profy,n,p,sigma);
    plot(ax2,n,fitprof,'k--')
    t1(end+1:end+length(fittxt))=fittxt;
    

    n=mx-linel-binwidth:binwidth:mx+linel+binwidth;
    jervistemp.binwidth=binwidth;
    profx=hist(x,n);profx([1 end])=[];n([1 end])=[];
    jervistemp.profx=profx;
    jervistemp.nx=n;
    axes(ax3)
    plot(ax3,n,profx);
    hold(ax3,'on')
    xlabel(ax3,'Position along line ROI (nm)')
    ylabel(ax3,'counts')
    
    
    fwhm=getFWHM(profx,n);
    t2{end+1}=['FWHM: ' 9 num2str(fwhm)];
     [fitp,fitprof,fittxt]=fitgeneralprofile(profx,n,p,sigma);
    plot(ax3,n,fitprof,'k--')
    t2(end+1:end+length(fittxt))=fittxt;
    
    if ~isempty(locs.znm)
        minzh=max(-750,min(z));
        maxzh=min(750,max(z));
    n=minzh-3*binwidth:binwidth:maxzh+3*binwidth;
    profz=hist(z,n);profz([1 end])=[];n([1 end])=[];
    axes(ax4)
    plot(ax4,n,profz);
    hold(ax4,'on')
    xlabel(ax4,'z (nm)')
    ylabel(ax4,'counts')
    fwhm=getFWHM(profz,n);
    t3{end+1}=['FWHM: ' 9  num2str(fwhm)];
    [fitp,fitprof,fittxt]=fitgeneralprofile(profz,n,p,fwhm/2.6);
    plot(ax4,n,fitprof,'k--')
    t3(end+1:end+length(fittxt))=fittxt;
    
    jervistemp.profz=profz;jervistemp.nz=n;
    
    
    axes(ax5)
    plot(x,z,'.')
    xlabel(ax5,'Position along line ROI (nm)')
    ylabel(ax5,'z (nm)')
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


