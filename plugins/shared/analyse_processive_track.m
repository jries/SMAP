function out=analyse_processive_track(x,y,time)

x=double(x);y=double(y);time=double(time-time(1));%angle=-double(angle)+pi;
angle=atan2(y(end)-y(1),x(end)-x(1));
fitfun=@(par,time) horzcat(par(1)+par(2)*cos(par(3))*time, par(4)+par(2)*sin(par(3))*time);
v0=sqrt((x(end)-x(1))^2+(y(end)-y(1))^2)/(time(end)-time(1));
startpar=[x(1),v0,angle,y(1)];
startxy=fitfun(startpar,time);

huber_loss = @(r, delta) (abs(r) <= delta) .* (0.5 * r.^2) + ...
                         (abs(r) > delta) .* (delta * (abs(r) - 0.5 * delta));
delta=1;
 % c=1;
% cauchy_loss = @(r) (c^2 / 2) * log(1 + (r / c).^2);

options = optimoptions('lsqnonlin', 'Display', 'off');

data=horzcat(x,y);
dataf=data; timef=time;

% % piecewise linear fits
    lend=length(x);
fr=0.6;
if lend>12/fr
    lsq_loss=@(r) r;

    ind1=(1:round(lend*fr))';
    objective_lsq = @(params) lsq_loss(fitfun(params, timef(ind1)) - dataf(ind1,:));
    [fitp1,r1,residuals1] = lsqnonlin(objective_lsq, startpar, [], [], options);
    
    ind2=(round(lend*(1-fr))+1:lend)';
    objective_lsq = @(params) lsq_loss(fitfun(params, timef(ind2)) - dataf(ind2,:));
    [fitp2,r2,residuals2] = lsqnonlin(objective_lsq, startpar, [], [], options);

    if r2/length(ind2)>r1/length(ind1)
        startpar=fitp1;
    else
        startpar=fitp2;
    end
    res=fitfun(startpar,timef)-dataf;
    r2d=sqrt(sum(res.^2,2));
    [m,stdrob,~,outlier]=robustMean(r2d(:));
    outlierb=false(size(timef));
    outlierb(outlier)=true;

    %fill in single zeros
    nn=vertcat(false, outlierb(3:end) & outlierb(1:end-2), false);
    outlierb=outlierb | nn;

    ind1=find(~outlierb,1,'first');
    ind2=find(~outlierb,1,'last');
    outlierb(ind1:ind2)=false;
    dataf(outlierb,:)=[]; timef(outlierb,:)=[];
    % remove first or last outliers. 
end
% %window: minimum 10 points, minimum 30%
% winlen=max(8, round(lend*0.3));
% % sliding window step: 3 points or 5%
% winstep=max(3, ceil(lend*0.05));
% numsteps=round((lend-winlen)/winstep);
% 
% r=zeros(numsteps,1);
% for k=1:numsteps
%     ind=((k-1)*winstep+1:min((k-1)*winstep+winlen,lend))';
%     objective_lsq = @(params) lsq_loss(fitfun(params, timef(ind)) - dataf(ind));
%     [fitp,r(k),residuals] = lsqnonlin(objective_lsq, starpar, [], [], options);
% end
% 

% fit from both ends, choose which is better


% end

% objective = @(params) cauchy_loss(fitfun(params, time) - data);

% Use lsqnonlin to solve the robust nonlinear least squares problem

for k=1:3
    objective = @(params) huber_loss(fitfun(params, timef) - dataf, delta);
    % objective = @(params) cauchy_loss(fitfun(params, timef) - dataf);
    [fitp,r,residuals] = lsqnonlin(objective, startpar, [], [], options);
    r2d=sqrt(sum(residuals.^2,2));
    [m,stdrob,~,outlier]=robustMean(r2d(:));
    if length(timef)-length(outlier) > 7 && ~isempty(outlier)
        dataf(outlier,:)=[]; timef(outlier,:)=[];
    % else 
        break
    end

end

out.x0=fitp(1);out.v=fitp(2);out.angle=fitp(3); out.y0=fitp(4);out.rmse=sqrt(mean(residuals(:).^2));
dr=diff(residuals)/sqrt(2);
[m,stdrob,~,outlier]=robustMean(dr(:));
out.gof=out.rmse/stdrob;
out.numpoints=length(timef);
out.fraction=length(timef)/length(time);
out.len=sqrt(sum((dataf(end,:)-dataf(1,:)).^2));


fitxy=fitfun(fitp,time);
    [xr,yr]=rotcoord(x,y,fitp(3));
    [xf,yf]=rotcoord(fitxy(:,1),fitxy(:,2),fitp(3));
out.plot.x=x; out.plot.y=y; out.plot.xfit=fitxy(:,1); out.plot.yfit=fitxy(:,2);
out.plot.xr=xr; out.plot.yr=yr; out.plot.xfitr=xf; out.plot.yfitr=yf;
out.plot.time=time;
out.filtered.time=timef;
out.filtered.x=dataf(:,1);out.filtered.y=dataf(:,2);
[out.filtered.xr,out.filtered.yr]=rotcoord(out.filtered.x,out.filtered.y,out.angle);
end
