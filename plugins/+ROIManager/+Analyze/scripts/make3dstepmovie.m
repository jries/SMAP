%plot 3D tracks, export movie.
%get ind of tracks
frametime=10; %real t ime

l1=g.locData.getloc({'groupindex'},'Position',g.locData.SE.sites(1));
gi(1)=mode(l1.groupindex);

l2=g.locData.getloc({'groupindex'},'Position',g.locData.SE.sites(2));
gi(2)=mode(l2.groupindex);


waitframes=30;

f=figure(77);
f.Color=[0 0 0];
ax=gca;
delete(ax.Children)
axis(ax,'equal');
view(ax,27,30)

hl=plot3(ax,xh,yh,zh,'k');
grid(ax,'on');
hold(ax,'on');
delete(hl)

ax.XTickLabel=[];ax.YTickLabel=[];ax.ZTickLabel=[];
ax.XLim=[680 1150];
ax.YLim=[2600 2950];
ax.ZLim=[-400 -50];

tp=[640 2700 -50];
cxo=[1,0.2,0.2];
clear Fr
c{1}=[0,0,1];
cs{1}=[0.5,0.5,1];

c{2}=[0,.5,0];
cs{2}=[0,1,0];

toff=0;
froff=0;
for tr=1:2
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
    hy=plot3(ax,xh,0*yh+ax.YLim(2),zh,'Color',cs{tr});
    hx=plot3(ax,0*xh+ax.XLim(1),yh,zh,'Color',cs{tr});
    hd=plot3(ax,xh(end),yh(end),zh(end),'ro','MarkerFaceColor','r','MarkerSize',15);

    hzo=plot3(ax,xh(end),yh(end),1+ax.ZLim(1),'o','MarkerFaceColor',cxo,'MarkerSize',15);
    hyo=plot3(ax,xh(end),-1+ax.YLim(2),zh(end),'o','MarkerFaceColor',cxo,'MarkerSize',15);
    hxo=plot3(ax,1+ax.XLim(1),yh(end),zh(end),'o','MarkerFaceColor',cxo,'MarkerSize',15);


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
