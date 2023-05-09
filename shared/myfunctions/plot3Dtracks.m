function plot3Dtracks(ax,x,y,z,axoffset, pin)
if nargin<5
    axoffset=[0 0 0];
end
if nargin<6
    pin.plotsides=true;
end
if ~isfield(pin,'skipfirst')
    pin.skipfirst=0;
end
if ~isfield(pin,'use') || isempty(pin.use)
    usetracks=true(size(x));
else
    usetracks=pin.use;
end
%plot 3D tracks, 
% pin: 
% 
% f=figure(77);
% f.Color=[0 0 0];
% ax=gca;
delete(ax.Children)
axis(ax,'equal');
view(ax,52,21)

% hl=plot3(ax,xh,yh,zh,'k');
grid(ax,'on');
hold(ax,'on');
ax.Projection='perspective';
ax.Projection="orthographic";
% delete(hl)

% ax.XTickLabel=[];ax.YTickLabel=[];ax.ZTickLabel=[];
xmax=0; xmin=inf;ymax=0; ymin=inf;zmax=0; zmin=inf;
for k=1:length(x)
    xmax=max(xmax,max(x{k}));xmin=min(xmin,min(x{k}));
    ymax=max(ymax,max(y{k}));ymin=min(ymin,min(y{k}));
    zmax=max(zmax,max(z{k}));zmin=min(zmin,min(z{k}));
end
if axoffset(1)<=0
    xplane=xmin+axoffset(1);
else
    xplane=xmax+axoffset(1);
end
if axoffset(2)<=0
    yplane=ymin+axoffset(2);
else
    yplane=ymax+axoffset(2);
end
if axoffset(3)<=0
    zplane=zmin+axoffset(3);
else
    zplane=zmax+axoffset(3);
end

% ax.XLim=[680 1150];
% ax.YLim=[2600 2950];
% ax.ZLim=[-400 -50];

% tp=[640 2700 -50];
% cxo=[1,0.2,0.2];
% clear Fr
% c{1}=[0,0,1];
% cs{1}=[0.5,0.5,1];
% 
% c{2}=[0,.5,0];
% cs{2}=[0,1,0];
% m=round(length(x)/3);
m=0;
c=lines(length(x)+m);
c=turbo(length(x)+m);
c(1:m,:)=[];



c=c(pin.indsort,:);

%make unique
% c=c+(1:length(x))'/100/length(x);
% c(c>1)=1;
cs=min(1,(c+0.8)*0.5);
cs=c;

skipfirst=pin.skipfirst;

for tr=1:length(x)
    if ~usetracks(tr)
        continue
    end
% ind=g.locData.loc.groupindex==trackids(tr);
% x=g.locData.loc.xnm(ind);
% y=g.locData.loc.ynm(ind);
% z=g.locData.loc.znm(ind);
% time=g.locData.loc.time(ind);


    
    xh=x{tr}(skipfirst+1:end);
    yh=y{tr}(skipfirst+1:end);
    zh=z{tr}(skipfirst+1:end);

    if pin.plotsides
    lw=2;
    else
        lw=0.7;
    end
    hl=plot3(ax,xh,yh,zh,'Color',c(tr,:),'LineWidth',lw);
    if pin.plotsides
        hz=plot3(ax,xh,yh,zplane+0*xh,'Color',cs(tr,:),'LineWidth',0.55);
        hy=plot3(ax,xh,yplane+0*xh,zh,'Color',cs(tr,:),'LineWidth',0.55);
        hx=plot3(ax,xplane+0*xh,yh,zh,'Color',cs(tr,:),'LineWidth',0.55);
    end
%     hd=plot3(ax,xh(end),yh(end),zh(end),'ro','MarkerFaceColor','r','MarkerSize',15);

%     hzo=plot3(ax,xh(end),yh(end),1+ax.ZLim(1),'o','MarkerFaceColor',cxo,'MarkerSize',15);
%     hyo=plot3(ax,xh(end),-1+ax.YLim(2),zh(end),'o','MarkerFaceColor',cxo,'MarkerSize',15);
%     hxo=plot3(ax,1+ax.XLim(1),yh(end),zh(end),'o','MarkerFaceColor',cxo,'MarkerSize',15);




end
axis(ax,'tight')
plot3(ax,xplane*[1 1],yplane*[1 1],ax.ZLim,'k')
plot3(ax,xplane*[1 1],ax.YLim,zplane*[1 1],'k')
plot3(ax,ax.XLim,yplane*[1 1],zplane*[1 1],'k')

    grid(ax,'on')
