function f=corners2labeling(corners,p)
% p.fitrange, p.rings, p.corners, p.axis
if ~isfield(p,'corners')
    p.corners=8;
end
if ~isfield(p,'rings')
    p.rings=4;
end

hi=histcounts(corners,0:1:p.corners+1);
n=0:p.corners;
shi=sum(hi);
% hi=hi/shi;

corners=p.corners;rings=p.rings;
if isfield(p,'fitrange')
n1=find(n>=p.fitrange(1),1,'first');
n2=find(n<=p.fitrange(2),1,'last');
else
    n1=1;
    n2=corners+1;
end
range=n1:n2;

x=n(range)';
% clusterfromlabeling(x,corners,rings,.5)
ft=fittype('a*clusterfromlabeling(x,corners,rings,p)','problem',{'corners','rings'});
f=fit(x,hi(range)',ft,'problem',{corners, rings},'Lower',[0 0.01],'Upper',[inf .99],'Start',[shi .4]);
if isfield(p,'axis')
hold(p.axis,'off')

bar(n,hi)
hold(p.axis,'on')
plot(p.axis,n,f(n),'-g')
plot(p.axis,x,f(x),'-*r')
end
end