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
if isfield(obj.locData.pluginresults.Analyze.tracking,"Single_tracking_analysis")
    singletracks=obj.locData.pluginresults.Analyze.tracking.Single_tracking_analysis;
else
    error("Run Single_trackin_analysis first")
end
%filters
maxd=p.maxd; % nm max distance to associate localizations
layers=find(obj.getPar('sr_layerson'));
obj.locData.filter; %does this fix the bug?
[locs,indin]=obj.locData.getloc({'xnm','ynm','znm','xpix','ypix','frame','track_id','channel','track_length','layer','filenumber'},'layer',layers,'position','roi','grouping','ungrouped');

exposuretime=obj.locData.files.file(locs.filenumber(1)).info.exposure;
roi=obj.locData.files.file(locs.filenumber(1)).info.roi;
pixelsize=obj.locData.files.file(locs.filenumber(1)).info.cam_pixelsize_um*1000;

trackids=singletracks.trackids; %long tracks, id>0
usetracks=unique(trackids);
trackstat=singletracks.tracks;

trackstat.partnerids=zeros(max(usetracks),1);
trackstat.partnertrackid=zeros(max(usetracks),1);

% link tracks
track_singlechannel=ismember(locs.track_id,usetracks); %(locs.track_id>0, long);
inchannel1=locs.channel==1;
indr=find((track_singlechannel&inchannel1));
locr.x=locs.xnm(indr);
locr.y=locs.ynm(indr);
locr.frame=locs.frame(indr);
locr.track_id=locs.track_id(indr);

inchannel2=locs.channel==2;
indt=find((track_singlechannel&inchannel2));
loct.x=locs.xnm(indt);
loct.y=locs.ynm(indt);
loct.frame=locs.frame(indt);
loct.track_id=locs.track_id(indt);

partner=zeros(size(locs.xnm));

if ~isempty(indt) && ~isempty(indr)
    [iAa,iBa,nA,nB,nseen]=matchlocsall(locr,loct,0,0,maxd);
    
    trackstat.partnertrackid(locs.track_id(indr(iAa)))=(locs.track_id(indt(iBa)));
    trackstat.partnertrackid(locs.track_id(indt(iBa)))=(locs.track_id(indr(iAa)));
   
    partner(indr(iAa))=indt(iBa);
    partner(indt(iBa))=indr(iAa);
end

haspartner=false(size(trackstat.velocity));
% processivepartner=false(size(trackstat.velocity));
trackstat.channel=zeros(size(trackstat.velocity));
cotracklen=zeros(size(trackstat.velocity));
dx0=zeros(size(trackstat.velocity))+inf;
rmse=zeros(size(trackstat.velocity))+inf;
dangle=zeros(size(trackstat.velocity))+2*pi;
bothprocessive=false(size(trackstat.velocity));
%
%I should do this only for channel 1, right?
for k=1:length(usetracks)  
    idref=usetracks(k);
    indref=find(locs.track_id==idref);
    trackstat.channel(idref)=mode(locs.channel(indref));
    if ~trackstat.processive(idref)
            continue
    end
    
    %look at co-movement
    haspartner(idref)=trackstat.partnertrackid(idref)>0;
    if haspartner(idref)  %&& validstats
        indpartner=partner(indref);indpartner=indpartner(indpartner>0);
        partnerids=locs.track_id(indpartner);
        [idpartner,npart]=mode(partnerids);
        trackstat.partnerids(idref)=idpartner;
        if ~trackstat.processive(idpartner)
            continue
        end
        bothprocessive(idref)=true;
        cotracklen(idref)=npart;
        % here look at fits and compare:
        % 1. angles need to be similar
        dangle(idref)=mod(trackstat.angle(idref)-trackstat.angle(idpartner)+pi,2*pi)-pi;
        % 2. average distance between line: did not work, use x0
        [fboth, i1,i2]=intersect(locs.frame(indref),locs.frame(indpartner));
        l1x=trackstat.x0(idref);
        l1y=trackstat.y0(idref);
        l2x=trackstat.x0(idpartner);
        l2y=trackstat.y0(idpartner);  
        dx0(idref)=sqrt(sum((l1x-l2x).^2+(l1y-l2y).^2));
        % 3. raw data: calculate dx(t), dy(t), average distance      
        rmse(idref)=sqrt(mean((locs.xnm(indref(i1))-locs.xnm(indpartner(i2))).^2+(locs.ynm(indref(i1))-locs.ynm(indpartner(i2))).^2));
    end
end
comovement=bothprocessive & abs(dangle)<p.dangle/180*pi & dx0<p.dx0 & rmse < p.rmsedat & cotracklen > p.cotracklength;
trackstat.comovement=comovement;
trackstat.processivepartner=bothprocessive; %counts each 

% plot tracks
ax=obj.initaxis('xy');

% cols=[1 0 1
%       0 1 1
%       1 0 0 
%       0 0 1
%       .5 0.2 0.2
%       0.2 0.2 .5];

% colind=trackstat.channel+trackstat.processive*2+trackstat.comovement*2;  

for k=1:length(usetracks)  
    idref=usetracks(k);
    indref=locs.track_id==idref;  
    % msize=3;
    % lw=1;
    % symb='-';
    % if trackstat.comovement(idref)
    %     if trackstat.channel(idref)==2
    %             symb='x-';
    %             msize=3;
    %     else
    %             symb='+-';
    %             msize=7;
    %             lw=2;            
    %     end
    % end

    if trackstat.comovement(idref) %|| trackstat.processive(idref)
        if trackstat.channel(idref)==1
            hp=plot(ax,locs.xnm(indref)/pixelsize(1)-roi(1),locs.ynm(indref)/pixelsize(2)-roi(2),'+-','Color',[1 0 0],'LineWidth',3,'Tag','test','MarkerSize',7);
        elseif trackstat.channel(idref)==2
            hp=plot(ax,locs.xnm(indref)/pixelsize(1)-roi(1),locs.ynm(indref)/pixelsize(2)-roi(2),'+-','Color',[0 0 1],'LineWidth',2,'Tag','test','MarkerSize',3);
        end
        pidlabel=0*locs.track_id(indref)+trackstat.partnerids(idref);
        dtRows = [dataTipTextRow("frame",double(locs.frame(indref))),...
        dataTipTextRow("ID",double(locs.track_id(indref))),...
        dataTipTextRow("partnerID",double(pidlabel))];
        alldatatip=vertcat(hp.DataTipTemplate.DataTipRows,dtRows');
        %hp.DataTipTemplate.DataTipRows(end+1:end+3) = dtRows;   
        hp.DataTipTemplate.DataTipRows=alldatatip;
    elseif trackstat.processive(idref)
        if trackstat.channel(idref)==1
            plot(ax,locs.xnm(indref)/pixelsize(1)-roi(1),locs.ynm(indref)/pixelsize(2)-roi(2),'-','Color',[1 0 0])
        elseif trackstat.channel(idref)==2
            plot(ax,locs.xnm(indref)/pixelsize(1)-roi(1),locs.ynm(indref)/pixelsize(2)-roi(2),'-','Color',[0 0 1])
        end

    else
        if trackstat.channel(idref)==1
            plot(ax,locs.xnm(indref)/pixelsize(1)-roi(1),locs.ynm(indref)/pixelsize(2)-roi(2),'-','Color',[1 0 1])
        elseif trackstat.channel(idref)==2
            plot(ax,locs.xnm(indref)/pixelsize(1)-roi(1),locs.ynm(indref)/pixelsize(2)-roi(2),'-','Color',[0 1 1])
        end
    end
    hold(ax,"on")

end


    
    axis(ax,'equal');
    xlim(ax,[0 roi(3)])
    ylim(ax,[0 roi(4)/2])
axis(ax,'ij');

axr=obj.initaxis('rmse');
plot(axr,cotracklen(bothprocessive), rmse(bothprocessive),'bo',cotracklen(comovement), rmse(comovement),'r+')
xlabel(axr,'co track length (frames)')
ylabel(axr,'rmse dual color (nm)')

axf=obj.initaxis('fit');
plot(axf,dx0(bothprocessive), dangle(bothprocessive)*180/pi,'bo',dx0(comovement), dangle(comovement)*180/pi,'r+')
xlabel(axf,'distance fit x0 (nm)')
ylabel(axf,'distance fit angle (°)')
% Plot cotracks vs time
%only both good
if contains(p.showtraces.selection,'processive co-tracks')
   
    goodpairs=find(trackstat.comovement & trackstat.channel==1);
    figure;
    % numrows=ceil(length(goodpairs)/5);
    f=0;
    for k=1:length(goodpairs)
        if 2*k-f>30
            f=f+30;
            figure
        end
        subplot(5,6,2*k-1-f)
        hold off
        id1=locs.track_id==goodpairs(k);
        idpartner=trackstat.partnerids(goodpairs(k));
        id2=locs.track_id==idpartner;
        tmin=min(min(locs.frame(id1)),min(locs.frame(id2)));
        tmax=max(max(locs.frame(id1)),min(locs.frame(id2)));

        xc=mean(locs.xnm(id1));yc=mean(locs.ynm(id1)); %rotate around same center, ch1 is ref
        angleh=trackstat.angle(goodpairs(k));
    
        [x1,y1]=rotcoord(locs.xnm(id1)-xc,locs.ynm(id1)-yc,angleh);
        if x1(end)<x1(1)
            angleh=angleh+pi;
            [x1,y1]=rotcoord(locs.xnm(id1)-xc,locs.ynm(id1)-yc,angleh);
        end
        [x2,y2]=rotcoord(locs.xnm(id2)-xc,locs.ynm(id2)-yc,angleh);
        plot(locs.frame(id1),x1,'.-',locs.frame(id2),x2,'.-')
        hold on
        title(['f ' num2str(tmin) ':', num2str(tmax),', x ' num2str(xc/1000,'%3.0f') ', y ' num2str(yc/1000,'%3.0f'),', ID ' num2str(goodpairs(k)), ',' num2str(idpartner,'%3.0f')])  
        xlabel('time(frame)')
        ylabel('xrot (nm)')
        subplot(5,6,2*k-f)
        plot(locs.xnm(id1)-xc,locs.ynm(id1)-yc,'.-',locs.xnm(id2)-xc,locs.ynm(id2)-yc,'.-')
        title(['d\alpha ' num2str(dangle(goodpairs(k))*180/pi,'%1.0f'), ', dx ', num2str(dx0(goodpairs(k)),'%3.0f'), ', rmse ', num2str(rmse(goodpairs(k)),'%3.0f')]) 
        axis equal
        xlabel('x (nm)')
        ylabel('y (nm)')

    end
end




% Calculate statistics

% extracts filename from file path:
    % filePath = string(obj.getPar('lastSMLFile'));
fnum=find(obj.getPar('sr_layerson'),1,'first');
fl=obj.getPar('filelist_long').String;
filePath=string(fl{fnum});
    % Find all occurrences of the substring
    % slash_indices = strfind(filePath, '/');
    % % If the substring is found
    % if ~isempty(slash_indices)
    %     % Get the last occurrence index
    %     lastSlashIndex = slash_indices(end);
    %     % Extract the substring after the last occurrence
    %     fileName = extractAfter(filePath, lastSlashIndex);
    % else
    %     fileName = filePath;
    % end
[path,fileName]=fileparts(filePath); 
% new output:
output=(table(fileName,... 
    sum(trackstat.channel==1), ...
    sum(trackstat.channel==2), ...
    sum(trackstat.channel==1 & trackstat.processive), ...
    sum(trackstat.channel==2 & trackstat.processive), ...
    sum(trackstat.partnertrackid>0)/2, ...
    sum(trackstat.comovement & trackstat.channel==1),...
    'VariableNames', {'Filename','ch1', 'ch2', 'ch1_prog','ch2_prog', 'dual', 'dual_prog'}));

disp(output)
out.resultstable=output;
out.trackstat=trackstat;
% clipboard('copy',output)

%%
% %add to movie
% dx=3;
% dimmark=3;
% if p.makemovie && ~isempty(p.tiffile) 
%     if ~isempty(p.Tfile)
%         tt=load(p.Tfile).transformation;
%     else
%         tt=obj.locData.files.file.transformation;
%     end
%     il=imageloaderMM; il.attachPar(obj.P);
%     il.openi(p.tiffile);
%     il.prefit;   
%     numf=il.metadata.numberOfFrames+1;
%     if tt.info{1}.xrange==tt.info{2}.xrange
%         sx=ceil(il.metadata.Height/2);
%         splitvert=true;
%         imcomb=zeros(sx,il.metadata.Width,3,numf);
%     else
%         sx=ceil(il.metadata.Width/2);
%         splitvert=false;  
%         imcomb=zeros(il.metadata.Height,sx,3,numf);
%     end
% 
% 
% 
% 
%     for k=1:numf
%         img=double(il.getimage(k));
%         imgt=tt.transformImageToTarget(2,img,'pixel',il.metadata.roi);
%         imrgb(:,:,1)=img;imrgb(:,:,2)=imgt;imrgb(:,:,3)=0;
%         figure(99); image(imrgb/max(imrgb(:)));
%         if splitvert
%             imcomb(:,:,1,k)=img(1:sx,:);
%             % imcomb(:,:,3,k)=img(:,1:sx);
%             imcomb(:,:,2,k)=imgt(1:sx,:);     
% 
%             % imcomb(:,:,1,k)=img(sx:end,:);
%             % imcomb(:,:,2,k)=imgt(sx:end,:);  
%         else
%             imcomb(:,:,1,k)=img(:,1:sx);
%             % imcomb(:,:,3,k)=img(:,1:sx);
%             imcomb(:,:,2,k)=imgt(:,1:sx);    
%         end
%     end
%     il.close;
%     [fp,fn,ext]=fileparts(p.tiffile);
%     if ~exist([fp filesep 'overlays'],"dir")
%         mkdir([fp filesep 'overlays']);
%     end
%     fout=[fp filesep 'overlays' filesep 'combined_' fn ext];
%     % fout=strrep(p.tiffile,'.ome.tif','_combined.tif');
%     options.color=true;
%     saveastiff(squeeze(single(imcomb)),fout,options);
% 
% 
%     if p.addtracksmovie %&& ~isempty(p.tiffile)
%         % imstack=tiffreadVolume(p.tiffile);
%         % imstack=permute(imstack,[1 2 4 3]);
%         imstack=imcomb;
%         immax=max(imstack(:));
%         [lenx,leny]=size(imstack,[1,2]);
%         idcoprogress=find(trackstat.coprocessive);
%         for k=1:length(idcoprogress)
%             ind=locs.track_id==idcoprogress(k);
%             x=round(locs.xnm(ind)/pixelsize(1)-roi(1)); y=round(locs.ynm(ind)/pixelsize(2)-roi(2)); frame=locs.frame(ind);
%             x=min(lenx-dx,x);x=max(dx+1,x);y=min(leny-dx,y);y=max(dx+1,y);
%             for l=1:length(x)
%                 imstack(y(l)+dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
%                 imstack(y(l)-dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
%                 imstack(y(l)-dx:y(l)+dx, x(l)+dx,dimmark,frame(l))=immax;
%                 imstack(y(l)-dx:y(l)+dx, x(l)-dx,dimmark,frame(l))=immax;
% 
%             end
% 
%         end
%         fout2=[fp filesep 'overlays' filesep 'mark_' fn ext];
%         % fout2=strrep(fout,'_combined.tif','_mark.tif');
%         options.color=true;
%         saveastiff(squeeze(single(imstack)),fout2,options);
%     end
% end
end
             
function pard=guidef(obj)

% pard.minlenframest.object=struct('String','Min length track (frames)','Style','text');
% pard.minlenframest.position=[1,1];
% pard.minlenframest.Width=1.5;
% 
% pard.minlenframes.object=struct('String','5','Style','edit');
% pard.minlenframes.position=[1,2.5];
% pard.minlenframes.Width=.5;

pard.maxdt.object=struct('String','Max distance (nm)','Style','text');
pard.maxdt.position=[1,3];
pard.maxdt.Width=1.5;

pard.maxd.object=struct('String','200','Style','edit');
pard.maxd.position=[1,4.5];
pard.maxd.Width=.5;




pard.cotracklenght.object=struct('String','co track length (frames) >','Style','text');
pard.cotracklenght.position=[3,1];
pard.cotracklenght.Width=1.5;
pard.cotracklength.object=struct('String','4','Style','edit');
pard.cotracklength.position=[3,2.5];
pard.cotracklength.Width=.5;

pard.danglet.object=struct('String','diff angle (°) <','Style','text');
pard.danglet.position=[3,3];
pard.danglet.Width=1.5;
pard.dangle.object=struct('String','30','Style','edit');
pard.dangle.position=[3,4.5];
pard.dangle.Width=.5;

pard.rmsedatt.object=struct('String','rmse data <','Style','text');
pard.rmsedatt.position=[4,1];
pard.rmsedatt.Width=1.5;
pard.rmsedat.object=struct('String','200','Style','edit');
pard.rmsedat.position=[4,2.5];
pard.rmsedat.Width=.5;

pard.dx0t.object=struct('String','d x0 fit <','Style','text');
pard.dx0t.position=[4,3];
pard.dx0t.Width=1.5;
pard.dx0.object=struct('String','300','Style','edit');
pard.dx0.position=[4,4.5];
pard.dx0.Width=.5;



pard.showt.object=struct('String','Show:','Style','text');
pard.showt.position=[5,1];
pard.showt.Width=0.5;
pard.showtraces.object=struct('String',{{'none','processive co-tracks'}},'Style','popupmenu');
pard.showtraces.position=[5,1.5];
pard.showtraces.Width=1.5;

% pard.makemovie.object=struct('String','make movie','Style','checkbox');
% pard.makemovie.position=[6,1];
% pard.makemovie.Width=1.5;
% 
% pard.addtracksmovie.object=struct('String','add tracks to movie','Style','checkbox');
% pard.addtracksmovie.position=[6,3];
% pard.addtracksmovie.Width=1.3;
% 
% pard.tiffilet.object=struct('String','tif:','Style','text');
% pard.tiffilet.position=[7,1];
% pard.tiffilet.Width=0.5;
% 
% pard.tiffile.object=struct('String','','Style','edit');
% pard.tiffile.position=[7,1.5];
% pard.tiffile.Width=3.;
% 
% pard.tiffload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'tiffile','*.tif'}});
% pard.tiffload.position=[7,4.5];
% pard.tiffload.Width=0.5;
% 
% pard.Tfilet.object=struct('String','trafo:','Style','text');
% pard.Tfilet.position=[8,1];
% pard.Tfilet.Width=0.5;
% 
% pard.Tfile.object=struct('String','','Style','edit');
% pard.Tfile.position=[8,1.5];
% pard.Tfile.Width=3.;
% 
% pard.Tfload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'Tfile','*.mat'}});
% pard.Tfload.position=[8,4.5];
% pard.Tfload.Width=0.5;


pard.plugininfo.description=sprintf('co-tracking analysis');
pard.plugininfo.type='ProcessorPlugin';
end

