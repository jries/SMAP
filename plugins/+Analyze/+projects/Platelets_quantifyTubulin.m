classdef Platelets_quantifyTubulin<interfaces.DialogProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a model of choice
    properties
        calibrationlines
    end
    methods
        function obj=Platelets_quantifyTubulin(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','sr_pos','sr_size'};
            obj.showresults=true;
        end
%         function initGui(obj)
%             lp=obj.guihandles.loclist.Position
%             obj.guihandles.loclist.delete;
%             h=uitable('Position',lp,'ColumnEditable',true,'ColumnName',[],'RowName',[]);
%             d=cell(100,1)
%         end
        function out=run(obj,p)           
            results=evalplatelet(obj,p);
            out.clipboard=results;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function addmt_callback(a,b,obj)

[locs,~,hroi]=obj.locData.getloc({'xnmline','ynmline','znm'},'layer',1,'position','roi');

pos=hroi.getPosition;
lenroi=sqrt(sum((pos(2,:)-pos(1,:)).^2));
numl=length(locs.xnmline)/lenroi;

cs=struct('pos',pos,'locsperum',numl,'string',num2str(numl,'%6.0f'));
if isempty(obj.calibrationlines)
    obj.calibrationlines=cs;
else
    obj.calibrationlines(end+1)=cs;
end

lls={obj.calibrationlines(:).string};

obj.guihandles.loclist.String=lls;
end

function remove_callback(a,b,obj,rem)
% p=obj.guihandles.loclist;
% if rem
% lls=p.String;
% lls(p.Value)=[];
% obj.calibrationlines(p.Value)=[];
% else
%     lls={''};
    obj.calibrationlines=[];
    obj.guihandles.loclist.String={''};
% end
% p.String=lls;
% p.Value=max(min(p.Value,length(lls)),1);

% obj.setGuiParameters(struct('loclist',p));
end

function results=evalplatelet(obj,p)

%determine locs in central platelet
cutoff=1;
if p.userroi
    poss='roi';
else
    poss='fov';
end
locs=obj.locData.getloc({'xnm','ynm','znm'},'layer',1,'position',poss);
pixrec=10;
rangex=[p.sr_pos(1)-p.sr_size(1) p.sr_pos(1)+p.sr_size(1)];
rangey=[p.sr_pos(2)-p.sr_size(2) p.sr_pos(2)+p.sr_size(2)];
im1=myhist2(locs.xnm,locs.ynm,pixrec,pixrec,rangex,rangey);
im1f=imfilter(im1,fspecial('Gauss',15,5));
imbs=im1f>cutoff;
sim=size(imbs);
components=bwconncomp(imbs);
rp=regionprops(components);
area=[rp(:).Area];
for l=1:length(rp)
    cx(l)=rp(l).Centroid(2);
    cy(l)=rp(l).Centroid(1);
end

d=sqrt((sim(1)/2-cx).^2+(sim(2)/2-cy).^2);
am=max(area);
d(area<am/2)=inf;
[~,indregion]=min(d);
bwf=false(sim);
bwf(components.PixelIdxList{indregion})=1;

if p.userroi
    d(indregion)=inf;
    [dmin,indregion]=min(d);
    if ~isinf(dmin)
        bwf(components.PixelIdxList{indregion})=1;
    end
end

ax=initaxis(p.resultstabgroup,'mask');
imagesc(ax,bwf')
axis(ax,'ij')

xr=floor((locs.xnm-rangex(1))/pixrec)+1;
yr=floor((locs.ynm-rangey(1))/pixrec)+1;
xr(xr>sim(1))=sim(1);
yr(yr>sim(2))=sim(2);

linind=sub2ind(sim,xr,yr);
inring=bwf(linind);
locsinring=sum(inring);

% determine average number of locs/mt/um
lls=p.loclist;
locsum=cellfun(@str2num,lls);
locsummean=mean(locsum);
locsummeanr=robustMean(locsum);


%total microtubule length
mttotlenum=locsinring/locsummeanr
results{1}=mttotlenum;

%do angle analysis
xh=locs.xnm(inring);
yh=locs.ynm(inring);

xhc=locs.xnm(inring)-cx(indregion)*pixrec-rangex(1);
yhc=locs.ynm(inring)-cy(indregion)*pixrec-rangey(1);

[th,rho]=cart2pol(xhc,yhc);
radius=mean(rho);
startp=[radius,0,0,0];
lb=[0,-inf,-inf,-1];
ub=[inf,inf,inf,1];
[fitp,r]=implicitfit(@circle_implicit,startp(1:3),xhc,yhc,0*xhc,lb(1:3),ub(1:3));
% [fitp,r]=implicitfit(@ellipse_implicit,startp,xhc,yhc,0*xhc,lb,ub);


cxn=fitp(2);cyn=fitp(3);
xhcn=xhc-cxn;
yhcn=yhc-cyn;

pts=[xhcn,yhcn];
 quad = quadratic_fit_points(pts);
  [ cen, axs, angle ] = quadratic_center(quad);

ax=initaxis(p.resultstabgroup,'fit');
plot(ax,xhcn,yhcn,'.')
circle(0,0,fitp(1),'Parent',ax)

hold on
plot(ax,0,0,'*')
h = quadratic_plot(quad);


axis(ax,'ij')

[th,rho]=cart2pol(xhcn,yhcn);

% [ths,inds]=sort(th);
% rhos=roh(inds);
dangle=pi/32;
angles=-pi:dangle:pi;
w=1./rho;
locsangle=zeros(size(angles));
for k=1:length(angles)-1
    inangle=th>angles(k)&th<angles(k+1);
    locsangle(k)=sum(w(inangle))/dangle;
end
mtangle=locsangle/locsummeanr*1000;
mtangle(end)=mtangle(1);

mtn=mtangle/max(mtangle)*fitp(1);
polar(ax,angles,mtn)
title(['total MT length: ' num2str(mttotlenum,'%3.0f'), 'um, R= ' num2str(fitp(1),'%3.0f') 'nm, mean MT = ',num2str(mean(mtangle(1:end-1)),'%3.1f'), ' ,turns: ' num2str(mttotlenum/2/pi/fitp(1)*1000,'%2.1f')])

ax2=initaxis(p.resultstabgroup,'mt(theta)');
plot(ax2,angles,mtangle)
xlabel(ax2,'theta')
ylabel(ax2,'number of microtubules')



if p.userroi
    outside=mtangle<max(mtangle)/8;
else
    outside=false(size(mtangle));
end
thetafraction= sum(~outside)/length(outside);

title(['total MT length (um): ' num2str(mttotlenum) ', fraction fitted: ' num2str(thetafraction)])

%assemble output

ovax=obj.getPar('ov_axes');
ovaxh=ovax.copy;
f=figure('Position',[50 100 500 750]);
subplot(3,2,1);
axh=gca;
ovaxh.Parent=f;
ovaxh.Position=axh.Position;
axh.delete;

imh=obj.getPar('sr_image');
subplot(3,2,2)
imagesc(imh.rangex,imh.rangey,imh.image)
axis equal
hold on
for k=1:length(obj.calibrationlines)
    plot(obj.calibrationlines(k).pos(:,1),obj.calibrationlines(k).pos(:,2),'m')
end
xlim(imh.rangex)
ylim(imh.rangey)

rangexnew=rangex-rangex(1)-(rangex(2)-rangex(1))/2;
rangeynew=rangey-rangey(1)-(rangey(2)-rangey(1))/2;

impl=myhist2(xhcn,yhcn,pixrec,pixrec,rangexnew,rangeynew);
subplot(3,2,3)
hold off
imagesc(rangexnew,rangeynew,impl')
axis ij
% plot(xhcn,yhcn,'.')
hold on
circle(0,0,fitp(1))

plot(0,0,'*')
polar(angles,mtn)
h = quadratic_plot(quad);

subplot(3,2,4)
plot(angles,mtangle)



xlim([-pi pi])
xlabel('theta')
ylabel('number of microtubules')
%text
str=obj.getPar('filelist_short').String;
sca='line calibrations (locs/um): ';
for k=1:length(obj.calibrationlines)
    sca=[sca obj.calibrationlines(k).string ', '];
end
str{end+1}=sca;
str{end+1}=['total MT length: ' num2str(mttotlenum,'%3.0f') ' um'];
str{end+1}=['R = ' num2str(fitp(1),'%3.0f') ' nm'] ;
str{end+1}=['mean MT per turn: ',num2str(mean(mtangle(1:end-1)),'%3.1f')];
str{end+1}=['turns from total MT and R: ' num2str(mttotlenum/2/pi/fitp(1)*1000,'%2.1f')];

str{end+1}=['Ellip. fit: axs= ' num2str(axs) ', angle= ' num2str(angle)];

if p.userroi
    str{end+1}=['fraction inside = ' num2str(thetafraction)];
    str{end+1}=['total MT length extrapolated: ' num2str(mttotlenum/thetafraction,'%3.0f') ' um'];
    str{end+1}=['mean MT per turn corrected: ',num2str(mean(mtangle(~outside)),'%3.1f')];
    str{end+1}=['turns from total MT and R corrected: ' num2str(mttotlenum/2/pi/fitp(1)*1000/thetafraction,'%2.1f')];
end

ht=uicontrol('Style','text','Units','normalized','Position',[0.05 0 .9 0.3]);
ht.String=str;
ht.HorizontalAlignment='left';
%save pdf
[po,fn]=fileparts(p.outputfile);
k=1;
while exist([po filesep fn num2str(k) '.pdf'],'file')
    k=k+1;
end
fno=[po filesep fn num2str(k) '.pdf'];
export_fig(gcf,fno,'-pdf');


results={locsummeanr,mttotlenum,fitp(1),axs(1),axs(2),mean(mtangle(1:end-1)),mean(mtangle(~outside)),thetafraction};
%calibration:locs / um, total MT length (um), radius, long axis, short
%axis, mean number of MT per turn (from MT(phi)), mean number of MT per turn averaged over ROI only, fraction of angles in
%ROI.
end

function outputfileb_callback(a,b,obj)
fn=obj.getSingleGuiParameter('outputfile');
if isempty(strfind(fn,'.pdf'))
    fn=[fn '*.pdf'];
end
[f,p]=uiputfile(fn);
if f
    obj.setGuiParameters(struct('outputfile',[p f]));
end
end
function loclist_callback(a,b,obj)
locsstored=[obj.calibrationlines(:).locsperum];
call=obj.calibrationlines;
used=false(size(call));
locslist=(obj.guihandles.loclist.String);
for k=1:length(locslist)
    locsum=str2double(locslist{k});
    ind=find((round(locsum)==round(locsstored)&~used));
    if ~isempty(ind) %already in list
        used(ind(1))=true;
    else %pasted: new data
        call(end+1)=struct('pos',zeros(2),'locsperum',locsum,'string',num2str(locsum));
        used(end+1)=true;
        locsstored(end+1)=locsum;
    end
    
end

call(~used)=[];
obj.guihandles.loclist.String={call(:).string};
obj.calibrationlines=call;
end

function pard=guidef(obj)
pard.addmt.object=struct('String','add single MT (locs/um)','Style','pushbutton','Callback',{{@addmt_callback,obj}});
pard.addmt.position=[1,1];
pard.addmt.Width=1.5;

pard.loclist.object=struct('String',{{''}},'Style','edit','Max',100,'Callback',{{@loclist_callback,obj}});
pard.loclist.position=[4,1];
pard.loclist.Width=1.5;
pard.loclist.Height=3;

% pard.remove.object=struct('String','remove','Style','pushbutton','Callback',{{@remove_callback,obj,1}});
% pard.remove.position=[5,1];
% pard.remove.Width=.75;

pard.clear.object=struct('String','clear','Style','pushbutton','Callback',{{@remove_callback,obj,0}});
pard.clear.position=[5,1.75];
pard.clear.Width=.75;

pard.userroi.object=struct('String','use ROI','Style','checkbox');
pard.userroi.position=[2,4];


pard.outputfile.object=struct('String','','Style','edit');
pard.outputfile.position=[7,1];
pard.outputfile.Width=3;

pard.outputfileb.object=struct('String','select','Style','pushbutton','Callback',{{@outputfileb_callback,obj}});
pard.outputfileb.position=[7,4];
pard.outputfileb.Width=1;

pard.plugininfo.name='Platelets_quantifyTubulin';
% pard.plugininfo.description=sprintf('Calculates profiles along a linear ROI and fits it with a model of choice. \n Flat: step function convolved with Gaussian (=Erf). \n Disk: Projection of a homogeneously filled disk, convolved with Gaussian. \n Ring: Projection of a ring, convolved with Gaussian. \n Distance: Two Gaussians in a distance d.');
pard.plugininfo.type='ProcessorPlugin';
end