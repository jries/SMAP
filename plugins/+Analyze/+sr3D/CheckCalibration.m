classdef CheckCalibration<interfaces.DialogProcessor
    properties
        zold
    end
    methods
        function obj=CheckCalibration(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
             obj.zold.changed=0;

        end

        function out=run(obj,p)
            locsall=obj.locData.getloc({'frame','xnm','ynm','znm','filenumber','phot'});
            maxfilen=max(locsall.filenumber);
            
            maxd=p.cam_pixelsize_nm(1);
            try
                p.refractiveIndexFactor=obj.locData.history{1}.children.fitparamters.fitterGUI.children.MLE_GPU_Yiming.refractive_index_mismatch;
            catch
                p.refractiveIndexFactor=1;
            end
            ax=obj.initaxis('filtered');
            locsall.beadnum=zeros(size(locsall.xnm));
              locsROI=obj.locData.getloc({'frame','xnm','ynm','znm','phot','filenumber'},'layer',1,'grouping','ungrouped','position','roi','removeFilter',{'filenumber'});
            for k=1:maxfilen
              
                fileroi=locsROI.filenumber==k;
                beadlocs=getBeadLocs(locsROI.xnm(fileroi),locsROI.ynm(fileroi),p);
                infile=locsall.filenumber==k;
                [bn,numlocs]=associatelocs(beadlocs.x,beadlocs.y,locsall.xnm(infile),locsall.ynm(infile),maxd);
                goodnum=bn>0;
                locsall.beadnum(infile)=bn+max(locsall.beadnum).*goodnum;
            end
            plotmanybeads(locsall,p)
%             plotzvsframe(locsROI,p)
            
            out=0;
        end
        
        function pard=guidef(obj)
            pard=guidef;
        end
        

    end
end

function plotmanybeads(locs,p)
mmax=max(locs.beadnum);
% indf1=locs.filenumber==1;
% indf2=locs.filenumber==2;
ax=initaxis(p.resultstabgroup,'zplots');
axis(ax)
hold off
for k=1:mmax
    indb=locs.beadnum==k;
    bf1=indb;
%     bf2=indb&indf2;
    locsb1.frame=locs.frame(bf1);
    locsb1.znm=locs.znm(bf1);
    locsb1.phot=locs.phot(bf1);
    
%     locsb2.frame=locs.frame(bf2);
%     locsb2.znm=locs.znm(bf2);
%     locsb2.phot=locs.phot(bf2);
    pf=plotzvsframe(locsb1,p);
    slope1(k)=pf(1);off1(k)=-pf(2)/pf(1);
    hold on
%     pf=plotzvsframe(locsb2,p);
%     slope2(k)=pf(1);off2(k)=-pf(2)/pf(1);
%     hold off
%     
    mx1(k)=mean(locs.xnm(bf1));
    my1(k)=mean(locs.ynm(bf1));
%     drawnow
%      waitforbuttonpress
end
initaxis(p.resultstabgroup,'slope vs z')
 indb=abs(slope1)>1.3|abs(off1)>10000|abs(slope1)<.8;
% indb=false(size(slope1));
indb=slope1==0;
plot(off1(~indb),slope1(~indb),'ro')


initaxis(p.resultstabgroup,'slope in image')
scatter3(mx1(~indb),my1(~indb),slope1(~indb),[],slope1(~indb))
xlabel('x position')
ylabel('y position')
zlabel('slope')
end


function pf=plotzvsframe(locs,p)

% diffz=10;
% refractiveIndexFactor=1.25;

z=locs.frame*p.dz*p.refractiveIndexFactor;
plot(z,locs.znm,'y.')
% 
% nd=7;
% dz=locs.znm;
% for k=1:nd
%     dz=diff(dz);
% end
fs=200/p.dz;
findmax=filter(1/fs*ones(fs,1),1,locs.phot./(abs(locs.znm)+p.dz));
 [~,ind]=max(findmax);
 
 hold on
% plot(z(ind),locs.znm(ind),'kd')
 
 
 
%   [~,ind]=min(abs(dz))
% nd=0;


dz1=diff(locs.znm);

dzf=filter(1/fs*ones(fs,1),1,dz1/p.refractiveIndexFactor/p.dz);
i2=find((dzf(ind:end))<0.3|dzf(ind:end)>3,1,'first')+ind-2;
i1=find((dzf(1:ind-1))<0.3|(dzf(1:ind-1))>3,1,'last')+1;

% i2=find((dz1(ind:end))<-p.dz*3,1,'first')+ind-2;
% i1=find((dz1(1:ind-1))<-p.dz*3,1,'last')+1;

if isempty(i1)
    i1=1;
end
if isempty(i2)
    i2=length(locs.znm);
end
if i2-i1<8
    pf=zeros(2,1);
    return
end
zlin=locs.znm(i1:i2);
flin=z(i1:i2);

pf=polyfit(flin,zlin,1);


plot(flin,zlin,'.')
plot(flin,polyval(pf,flin),'r')
hold off
title(['slope: ' num2str(pf(1))])
end

function pard=guidef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];

pard.text2.object=struct('String','check bead calibration','Style','text');
pard.text2.position=[1,1];

pard.dzt.object=struct('String','dz (nm)','Style','text');
pard.dzt.position=[2,1];

pard.dz.object=struct('String','10','Style','edit');
pard.dz.position=[2,2];

pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.type='ProcessorPlugin';
end
