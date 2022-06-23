%plot 3D tracks, export movie.
%get ind of tracks
frametime=1; %real t ime

siten=8;
viewpoint=[-58 30];
xlimind=2;
ylimind=2;

% siten=20;
% viewpoint=[-138 30];
% xlimind=2;
% ylimind=1;
% 
% siten=32;
% viewpoint=[-228 33];
% xlimind=1;
% ylimind=1;

l1=g.locData.getloc({'groupindex','xnm','ynm','znm'},'Position',g.locData.SE.sites(siten));
gi(1)=mode(l1.groupindex);
% 
% l2=g.locData.getloc({'groupindex'},'Position',g.locData.SE.sites(2));
% gi(2)=mode(l2.groupindex);


waitframes=30;
xm=inf;xx=-inf; ym=inf;yx=-inf; zm=inf;zx=-inf; 
for tr=1:length(gi)
ind=g.locData.loc.groupindex==gi(tr);
xm=min(xm,min(g.locData.loc.xnm(ind)));xx=max(xx,max(g.locData.loc.xnm(ind)));
ym=min(ym,min(g.locData.loc.ynm(ind)));yx=max(yx,max(g.locData.loc.ynm(ind)));
zm=min(zm,min(g.locData.loc.znm(ind)));zx=max(zx,max(g.locData.loc.znm(ind)));
end


f=figure(78);
f.Color=[0 0 0]+0.;
ax=gca;
delete(ax.Children)
axis(ax,'equal');
% view(ax,27,30)

hl=plot3(ax,l1.xnm,l1.ynm,l1.znm,'k');
grid(ax,'on');
hold(ax,'on');
delete(hl)

ax.XTickLabel=[];ax.YTickLabel=[];ax.ZTickLabel=[];
axis(ax,'equal');
ax.XLim=[xm-40 xx+20];
ax.YLim=[ym-20 yx+20];
ax.ZLim=[zm-20 zx+20];
view(ax,viewpoint(1),viewpoint(2))

plot3(ax,ax.XLim, ax.YLim(ylimind)*[1 1], ax.ZLim(1)*[1 1],'k')
plot3(ax,ax.XLim(xlimind)*[1 1], ax.YLim, ax.ZLim(1)*[1 1],'k')
plot3(ax,ax.XLim(xlimind)*[1 1], ax.YLim(ylimind)*[1 1], ax.ZLim,'k')
% ax.XLim=[680 1150];
% ax.YLim=[2600 2950];
% ax.ZLim=[-400 -50];

tp=[640 2700 -50];
cxo=[1,0.2,0.2];
clear Fr
c{1}=[0,0,1];
cs{1}=[0.5,0.5,1];

c{2}=[0,.5,0];
cs{2}=[0,1,0];

toff=0;
froff=0;
for tr=1:length(gi)
ind=g.locData.loc.groupindex==gi(tr);
x=g.locData.loc.xnm(ind);
y=g.locData.loc.ynm(ind);
z=g.locData.loc.znm(ind);
time=g.locData.loc.time(ind);

ts=min(time):frametime:max(time);
for k=1:length(ts)
    indh=time<=ts(k);
    xh=x(indh);
    yh=y(indh);
    zh=z(indh);
    tpassed=ts(k)-ts(1)+toff;
    ht=text(ax,tp(1),tp(2),[num2str(tpassed,'%3.0f') ' ms'],'FontSize',23,'Color','w');
    hl=plot3(ax,xh,yh,zh,'Color',c{tr});

    hz=plot3(ax,xh,yh,0*zh+ax.ZLim(1),'Color',cs{tr});
    hy=plot3(ax,xh,0*yh+ax.YLim(ylimind),zh,'Color',cs{tr});
    hx=plot3(ax,0*xh+ax.XLim(xlimind),yh,zh,'Color',cs{tr});
    hd=plot3(ax,xh(end),yh(end),zh(end),'ro','MarkerFaceColor','r','MarkerSize',15);

    hzo=plot3(ax,xh(end),yh(end),1+ax.ZLim(1),'o','MarkerFaceColor',cxo,'MarkerSize',15);
    hyo=plot3(ax,xh(end),-1+ax.YLim(ylimind),zh(end),'o','MarkerFaceColor',cxo,'MarkerSize',15);
    hxo=plot3(ax,1+ax.XLim(xlimind),yh(end),zh(end),'o','MarkerFaceColor',cxo,'MarkerSize',15);


    grid(ax,'on')
    drawnow
    Fr(k+froff)=getframe(ax);
    if k<length(ts)
        delete(hd);delete(hx);delete(hy);delete(hz);delete(ht);delete(ht)
        delete(hxo);delete(hyo);delete(hzo);
    else
        delete(hd);delete(ht);
    end
end
Fr(k+froff:k+froff+waitframes)=Fr(k+froff);
toff=tpassed;
froff=k+waitframes;
end


smlfile=g.getPar('lastSMLFile');
if ~isempty(smlfile)
    pfad=fileparts(smlfile);
else
    pfad=fileparts(obj.locData.files.file(1).name);
end

[file,pfad]=uiputfile([pfad filesep '*.mp4']);
if file
    mysavemovie(Fr,[pfad  file],'FrameRate',30)
end 
