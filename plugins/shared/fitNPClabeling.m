function plabel=fitNPClabeling(corners,p)
% p.corners, rings,fitrange

    
 nb=0:p.corners;
 if length(corners)>9 %otherwise: histogram passed on
    hr=hist(corners,nb);
 else
     hr=corners;
 end
%  hold off
if ~isfield(p,'ploton') || p.ploton==true
 plot(nb(1:end),hr,'-x')
 hold on
 end
plabel=fithist(hr,nb,p);

% title(plabel)
end
function pf=fithist(hi,n,p)

shi=sum(hi);
% hi=hi/shi;

corners=p.corners;rings=p.rings;
n1=find(n>=p.fitrange(1),1,'first');
n2=find(n<=p.fitrange(2),1,'last');
range=n1:n2;
% x=(0:corners)';
x=n(range)';
% clusterfromlabeling(x,corners,rings,.5)
% ft=fittype('a*clusterfromlabeling(x,corners,rings,p)','problem',{'corners','rings'});
% f=fit(x,hi(range)',ft,'problem',{corners, rings},'Lower',[0 0.01],'Upper',[inf .99],'Start',[shi .4]);

ft=fittype('sqrt(a*clusterfromlabeling(x,corners,rings,p))','problem',{'corners','rings'});
f=fit(x,sqrt(hi(range))',ft,'problem',{corners, rings},'Lower',[0 0.01],'Upper',[inf .99],'Start',[shi .4]);

if ~isfield(p,'ploton') || p.ploton==true
plot(n,f(n).^2,'-g')
plot(x,f(x).^2,'-*r')
end

% if ~isfield(p,'ploton') || p.ploton==true
% plot(n,f(n),'-g')
% plot(x,f(x),'-*r')
% end
pf=f.p;
end