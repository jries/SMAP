function [xpos,fp,xerr]=mypeakfit(x,y)
% [fp,S]=polyfit(x,y,2);
fp=fit(double(x'),double(y'),'poly2','Upper',[0 inf inf]);
xpos=-fp.p2/2/fp.p1;
ci=confint(fp);
dc=(ci(2,:)-ci(1,:))/2;
dx=sqrt((1/fp.p1/2*dc(1,2)).^2+(fp.p2/2/fp.p1^2*dc(1,1)).^2);
xerr=dx;
end