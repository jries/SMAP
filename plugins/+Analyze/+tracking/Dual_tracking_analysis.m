classdef Dual_tracking_analysis<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    methods
        function obj=Dual_tracking_analysis(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=true;
        end
        function out=run(obj,p)
            
            out=runi(obj,p);
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function loadtiff(obj,a,b,field,sstr)
            oldf=obj.getSingleGuiParameter(field);
            if isempty(oldf)
                oldf=sstr;
            end

            [newf,newpath]=uigetfile(oldf);
            if newf
                obj.setGuiParameters(struct(field,[newpath newf]));
            end
        
        end
    end
end
function out=runi(obj,p)
out=[];
%filters
minlenframes=p.minlenframes;
maxd=p.maxd; % nm max distance to associate localizations

cotracklength=p.cotracklength; %minimum data poits associated between tracks
cotrackfraction = p.cotrackfraction; %minimum fraction of localizations associated

% test directed movement
aspectratio=p.aspectratio; %should be below this for directed motion.
lennmstartend=p.lennmstartend; % minimum endpoint- startpoint
lennmmin=0; %minimum of largest extenstion standard deviation (not so useful)
minvelocity=p.velocitymin; %nm/s
maxvelocity=p.velocitymax;

% visualizing co-tracks
% onlyprogressivecotracks=true; %only when both partners are progressive
% markintiffile=false;
% tiffile=fout;
% tiffile='/Users/ries/datalocal/2color_kinesin/25_50ms_561nm01_640nm02_600w52_676w37_1_MMStack_Default_combined.tif';
layers=find(obj.getPar('sr_layerson'));
% obj=g;
obj.locData.filter; %does this fix the bug?
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','xpix','ypix','frame','track_id','channel','track_length','layer','filenumber'},'layer',layers,'position','roi','grouping','ungrouped');

% unique(locs.channel)

exposuretime=obj.locData.files.file(locs.filenumber(1)).info.exposure;

% XXX if exposure is wront, overwrite here
% exposuretime= 10

% test directed movement, calculate statistics
usetracks=unique(locs.track_id(locs.track_id>0));
trackstat.lennmstartend=zeros(max(usetracks),1);
trackstat.stdshort=zeros(max(usetracks),1);
trackstat.stdlong=zeros(max(usetracks),1);
trackstat.lenframe=zeros(max(usetracks),1);
trackstat.partnertrackid=zeros(max(usetracks),1);
trackstat.velocity=zeros(max(usetracks),1);
locs.track_length_new=zeros(size(locs.xnm));
trackstat.partnerids=zeros(max(usetracks),1);

for k=1:length(usetracks)
    iduset=usetracks(k);
    tind=find(locs.track_id==iduset);
    trackstat.lenframe(iduset)=length(tind);
    locs.track_length_new(tind)=trackstat.lenframe(iduset);
    if trackstat.lenframe(iduset)<2
        continue
    end


    xh=locs.xnm(tind);yh=locs.ynm(tind);fh=locs.frame(tind);
    trackstat.lennmstartend(iduset)=sqrt((xh(end)-xh(1)).^2+(yh(end)-yh(1)).^2);
        
    c = cov(xh-mean(xh), yh-mean(yh));
    [a, ev] = eig(c);
    [ev,ind] = sort(diag(ev));
    [xa, ya] = deal(a(1,ind(end)), a(2,ind(end)));
    trackstat.angle(iduset)=cart2pol(xa, ya);
    trackstat.stdshort(iduset)=real(sqrt(ev(1)));
    trackstat.stdlong(iduset)=real(sqrt(ev(2)));
    trackstat.velocity(iduset)=trackstat.lennmstartend(iduset)/(exposuretime/1000*trackstat.lenframe(iduset)); %nm/s
           
end

intrack=(locs.track_id>0);

longtracks=locs.track_length_new>minlenframes;


inchannel=locs.channel==1;
sum(inchannel)


indr=find((intrack&longtracks&inchannel));
locr.x=locs.xnm(indr);
locr.y=locs.ynm(indr);
locr.frame=locs.frame(indr);
locr.track_id=locs.track_id(indr);


inchannel2=locs.channel==2;
sum(inchannel2)
indt=find((intrack&longtracks&inchannel2));
loct.x=locs.xnm(indt);
loct.y=locs.ynm(indt);
loct.frame=locs.frame(indt);
loct.track_id=locs.track_id(indt);

[iAa,iBa,nA,nB,nseen]=matchlocsall(locr,loct,0,0,maxd);

trackstat.partnertrackid(locs.track_id(indr(iAa)))=(locs.track_id(indt(iBa)));
trackstat.partnertrackid(locs.track_id(indt(iBa)))=(locs.track_id(indr(iAa)));

% pairedr=false(size(locs.xnm));
% pairedt=pairedr;
partner=zeros(size(locs.xnm));
% pairedr(indr(iAa))=true;
% pairedt(indt(iBa))=true;
partner(indr(iAa))=indt(iBa);
partner(indt(iBa))=indr(iAa);

usetracks=unique(locs.track_id(intrack&longtracks));

roi=obj.locData.files.file(locs.filenumber(1)).info.roi;
pixelsize=obj.locData.files.file(locs.filenumber(1)).info.cam_pixelsize_um*1000;


%look for progressive movement
validstats=true(size(trackstat.velocity));
validstats = validstats & trackstat.stdshort./trackstat.stdlong<aspectratio;
validstats = validstats & trackstat.lennmstartend > lennmstartend;
validstats = validstats & trackstat.stdlong > lennmmin;
validstats = validstats & trackstat.velocity > minvelocity;
validstats = validstats & trackstat.velocity < maxvelocity;

trackstat.progressive=validstats;

comovement=false(size(trackstat.velocity));
progressivepartner=false(size(trackstat.velocity));
trackstat.channel=zeros(size(trackstat.velocity));


for k=1:length(usetracks)  
    idh=usetracks(k);
    indtr=locs.track_id==idh;

    %look at co-movement
    pid=0;
    comovement(idh)=trackstat.partnertrackid(idh)>0;
    if comovement(idh) %&& validstats
        partnerind=partner(indtr);
        partnerids=locs.track_id(partnerind(partnerind>0));
        [pid,npart]=mode(partnerids);
        trackstat.partnerids(idh)=pid;
        
        lenpartner=sum(locs.track_id==pid);
        minduallength=min(sum(indtr),lenpartner);

        comovement(idh)=comovement(idh) & npart>cotracklength;
        comovement(idh)=comovement(idh) & npart/minduallength>cotrackfraction;
        progressivepartner(idh)=trackstat.progressive(pid);
    end
    trackstat.channel(idh)=mode(locs.channel(indtr));
    
end

trackstat.comovement=comovement;
trackstat.progressivepartner=progressivepartner;


% plot tracks

ax=obj.initaxis('xy');

cols=[1 0 1
      0 1 1
      1 0 0 
      0 0 1
      .5 0.2 0.2
      0.2 0.2 .5
      .7 0 0.7
      0 0.7 .7];

colind=trackstat.channel+trackstat.comovement*2+(~trackstat.progressive)*4;   

for k=1:length(usetracks)  
    idh=usetracks(k);
    indtr=locs.track_id==idh;
    
    msize=3;
    lw=1;
    symb='-';

    if trackstat.comovement(idh)
        if trackstat.channel(idh)==2
                symb='x-';
                msize=3;
        else
                symb='+-';
                msize=7;
                lw=2;            
        end
    end
    hp=plot(ax,locs.xnm(indtr)/pixelsize(1)-roi(1),locs.ynm(indtr)/pixelsize(2)-roi(2),symb,'Color',cols(colind(idh),:),'LineWidth',lw,'Tag','test','MarkerSize',msize);
    hold(ax,"on")
    pidlabel=0*locs.track_id(indtr)+trackstat.partnerids(idh);
    dtRows = [dataTipTextRow("frame",double(locs.frame(indtr))),...
    dataTipTextRow("ID",double(locs.track_id(indtr))),...
    dataTipTextRow("partnerID",double(pidlabel))];
    alldatatip=vertcat(hp.DataTipTemplate.DataTipRows,dtRows');
    %hp.DataTipTemplate.DataTipRows(end+1:end+3) = dtRows;   
    hp.DataTipTemplate.DataTipRows=alldatatip;
end


    axis(ax,'ij');
    axis(ax,'equal');

% Plot cotracks vs time
%%
%only both good
if contains(p.showtraces.selection,'progressive co-tracks')
    trackstat.coprogressive=(trackstat.channel==1 & trackstat.comovement & trackstat.progressive & trackstat.progressivepartner);
    goodpairs=find(trackstat.coprogressive);
    figure;
    % numrows=ceil(length(goodpairs)/5);
    f=0;
    for k=1:length(goodpairs)
        if 2*k-f>25
            f=f+30;
            figure
        end
        subplot(5,6,2*k-1-f)
        hold off
        id1=locs.track_id==goodpairs(k);
        pid=trackstat.partnerids(goodpairs(k));
        id2=locs.track_id==pid;
        tmin=min(min(locs.frame(id1)),min(locs.frame(id2)));
        tmax=max(max(locs.frame(id1)),min(locs.frame(id2)));
    
        [x1,y1]=rotcoord(locs.xnm(id1)-mean(locs.xnm(id1)),locs.ynm(id1)-mean(locs.ynm(id1)),trackstat.angle(goodpairs(k)));
        [x2,y2]=rotcoord(locs.xnm(id2)-mean(locs.xnm(id2)),locs.ynm(id2)-mean(locs.ynm(id2)),trackstat.angle(pid));
        plot(locs.frame(id1),x1,'.-',locs.frame(id2),x2,'.-')
        hold on
        title(['frame: ' num2str(tmin) ':', num2str(tmax),', x: ' num2str(mean(locs.xpix(id1)),'%3.0f') ', y: ' num2str(mean(locs.ypix(id1)),'%3.0f')])  
        xlabel('time(frame)')
        ylabel('xrot (nm)')
        subplot(5,6,2*k-f)
        plot(locs.xnm(id1)-mean(locs.xnm(id1)),locs.ynm(id1)-mean(locs.ynm(id1)),'.-',locs.xnm(id2)-mean(locs.xnm(id2)),locs.ynm(id2)-mean(locs.ynm(id2)),'.-')
        axis equal
        xlabel('x (nm)')
        ylabel('y (nm)')

    end
end

% Calculate statistics
disp(sprintf('ch1\tch2\tch1prog\tch2prog\tdualcol\tch12progr'));
output=(sprintf([num2str(sum(trackstat.channel==1)), '\t' num2str(sum(trackstat.channel==2)),...
    '\t' num2str(sum(trackstat.channel==1 & trackstat.progressive)), ...
    '\t' num2str(sum(trackstat.channel==2 & trackstat.progressive)), ...
    '\t' num2str(sum(trackstat.partnertrackid>0)/2), ...
    '\t' num2str(sum(trackstat.coprogressive))]));
disp(output)
clipboard('copy',output)

%%
%add to movie
dx=3;
dimmark=3;
if p.makemovie && ~isempty(p.tiffile) 
    if ~isempty(p.Tfile)
        tt=load(p.Tfile).transformation;
    else
        tt=obj.locData.files.file.transformation;
    end
    il=imageloaderMM; il.attachPar(obj.P);
    il.openi(p.tiffile);
    il.prefit;   
    numf=il.metadata.numberOfFrames+1;
    if tt.info{1}.xrange==tt.info{2}.xrange
        sx=ceil(il.metadata.Height/2);
        splitvert=true;
        imcomb=zeros(sx,il.metadata.Width,3,numf);
    else
        sx=ceil(il.metadata.Width/2);
        splitvert=false;  
        imcomb=zeros(il.metadata.Height,sx,3,numf);
    end

    
    
 
    for k=1:numf
        img=double(il.getimage(k));
        imgt=tt.transformImageToTarget(2,img,'pixel',il.metadata.roi);
        imrgb(:,:,1)=img;imrgb(:,:,2)=imgt;imrgb(:,:,3)=0;
        figure(99); image(imrgb/max(imrgb(:)));
        if splitvert
            imcomb(:,:,1,k)=img(1:sx,:);
            % imcomb(:,:,3,k)=img(:,1:sx);
            imcomb(:,:,2,k)=imgt(1:sx,:);             
        else
            imcomb(:,:,1,k)=img(:,1:sx);
            % imcomb(:,:,3,k)=img(:,1:sx);
            imcomb(:,:,2,k)=imgt(:,1:sx);    
        end
    end
    il.close;
    [fp,fn,ext]=fileparts(p.tiffile);
    fout=[fp filesep 'overlays' filesep 'combined_' fn ext];
    % fout=strrep(p.tiffile,'.ome.tif','_combined.tif');
    options.color=true;
    saveastiff(squeeze(single(imcomb)),fout,options);


    if p.addtracksmovie %&& ~isempty(p.tiffile)
        % imstack=tiffreadVolume(p.tiffile);
        % imstack=permute(imstack,[1 2 4 3]);
        imstack=imcomb;
        immax=max(imstack(:));
        [lenx,leny]=size(imstack,[1,2]);
        idcoprogress=find(trackstat.coprogressive);
        for k=1:length(idcoprogress)
            ind=locs.track_id==idcoprogress(k);
            x=round(locs.xnm(ind)/pixelsize(1)-roi(1)); y=round(locs.ynm(ind)/pixelsize(2)-roi(2)); frame=locs.frame(ind);
            x=min(lenx-dx,x);x=max(dx+1,x);y=min(leny-dx,y);y=max(dx+1,y);
            for l=1:length(x)
                imstack(y(l)+dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
                imstack(y(l)-dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
                imstack(y(l)-dx:y(l)+dx, x(l)+dx,dimmark,frame(l))=immax;
                imstack(y(l)-dx:y(l)+dx, x(l)-dx,dimmark,frame(l))=immax;
    
            end
    
        end
        fout2=[fp filesep 'overlays' filesep 'mark_' fn ext];
        % fout2=strrep(fout,'_combined.tif','_mark.tif');
        options.color=true;
        saveastiff(squeeze(single(imstack)),fout2,options);
    end
end
end
             
function pard=guidef(obj)

pard.minlenframest.object=struct('String','Min length track (frames)','Style','text');
pard.minlenframest.position=[1,1];
pard.minlenframest.Width=1.5;

pard.minlenframes.object=struct('String','5','Style','edit');
pard.minlenframes.position=[1,2.5];
pard.minlenframes.Width=.5;

pard.maxdt.object=struct('String','Max distance (nm)','Style','text');
pard.maxdt.position=[1,3];
pard.maxdt.Width=1.5;

pard.maxd.object=struct('String','200','Style','edit');
pard.maxd.position=[1,4.5];
pard.maxd.Width=.5;

pard.dmt.object=struct('String','directed movement:','Style','text');
pard.dmt.position=[2,1];

pard.aspectratiot.object=struct('String','aspect ratio <','Style','text');
pard.aspectratiot.position=[2,2];
pard.aspectratiot.Width=1.;
pard.aspectratio.object=struct('String','1','Style','edit');
pard.aspectratio.position=[2,2.7];
pard.aspectratio.Width=.5;


pard.lennmstartendt.object=struct('String','start-end (nm) >','Style','text');
pard.lennmstartendt.position=[2,3.5];
pard.lennmstartendt.Width=1.;
pard.lennmstartend.object=struct('String','300','Style','edit');
pard.lennmstartend.position=[2,4.5];
pard.lennmstartend.Width=.5;

pard.cotracklenght.object=struct('String','co track length (frames) >','Style','text');
pard.cotracklenght.position=[3,1];
pard.cotracklenght.Width=1.5;
pard.cotracklength.object=struct('String','4','Style','edit');
pard.cotracklength.position=[3,2.5];
pard.cotracklength.Width=.5;

pard.cotrackfractiont.object=struct('String','co track fraction >','Style','text');
pard.cotrackfractiont.position=[3,3];
pard.cotrackfractiont.Width=1.5;
pard.cotrackfraction.object=struct('String','0','Style','edit');
pard.cotrackfraction.position=[3,4.5];
pard.cotrackfraction.Width=.5;

pard.velocityt.object=struct('String','velocity (nm/s) min','Style','text');
pard.velocityt.position=[4,1];
pard.velocityt.Width=1.5;
pard.velocitymin.object=struct('String','200','Style','edit');
pard.velocitymin.position=[4,2.5];
pard.velocitymin.Width=.5;
pard.velocitymt.object=struct('String','max','Style','text');
pard.velocitymt.position=[4,4];
pard.velocitymt.Width=0.5;
pard.velocitymax.object=struct('String','1500','Style','edit');
pard.velocitymax.position=[4,4.5];
pard.velocitymax.Width=.5;

pard.showt.object=struct('String','Show:','Style','text');
pard.showt.position=[5,1];
pard.showt.Width=0.5;
pard.showtraces.object=struct('String',{{'progressive co-tracks','none'}},'Style','popupmenu');
pard.showtraces.position=[5,1.5];
pard.showtraces.Width=1.5;

pard.makemovie.object=struct('String','make movie','Style','checkbox');
pard.makemovie.position=[6,1];
pard.makemovie.Width=1.5;

pard.addtracksmovie.object=struct('String','add tracks to movie','Style','checkbox');
pard.addtracksmovie.position=[6,3];
pard.addtracksmovie.Width=1.3;

pard.tiffilet.object=struct('String','tif:','Style','text');
pard.tiffilet.position=[7,1];
pard.tiffilet.Width=0.5;

pard.tiffile.object=struct('String','','Style','edit');
pard.tiffile.position=[7,1.5];
pard.tiffile.Width=3.;

pard.tiffload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'tiffile','*.tif'}});
pard.tiffload.position=[7,4.5];
pard.tiffload.Width=0.5;

pard.Tfilet.object=struct('String','trafo:','Style','text');
pard.Tfilet.position=[8,1];
pard.Tfilet.Width=0.5;

pard.Tfile.object=struct('String','','Style','edit');
pard.Tfile.position=[8,1.5];
pard.Tfile.Width=3.;

pard.Tfload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'Tfile','*.mat'}});
pard.Tfload.position=[8,4.5];
pard.Tfload.Width=0.5;


pard.plugininfo.description=sprintf('co-tracking analysis');
pard.plugininfo.type='ProcessorPlugin';
end

