function fitpar=mygaussfit(xdat,ydat,start,fixbg)
if nargin<4
    fixbg=false;
end
%a x sigma bg
% size(xdat)
% size(ydat)
% 
 oldopts=optimset('lsqcurvefit');
 oldopts=optimset(oldopts,'display','off');
newopts=optimset(oldopts,'MaxFunEvals',200,'MaxIter',200,'TolFun',1e-6,'Algorithm','levenberg-marquardt');
% newopts=optimset('MaxFunEvals',600,'MaxIter',600,'TolFun',1e-7);
% newopts=oldopts;
% fitpar=lsqcurvefit(@mygaussforfit,double(start),double(xdat),double(ydat));%,[0 0 0 0],[inf inf inf max(ydat)/100],newopts);
if fixbg
fitpar=lsqcurvefit(@mygaussforfit,double(start),double(xdat),double(ydat),[-inf -inf 0 0],[inf inf inf 1e-7],oldopts);
else
fitpar=lsqcurvefit(@mygaussforfit,double(start),double(xdat),double(ydat),[-inf -inf 0 -inf],[inf inf inf inf],oldopts);
end
% fitpar=lsqcurvefit(@mygaussforfit,double(start),double(xdat),double(ydat));