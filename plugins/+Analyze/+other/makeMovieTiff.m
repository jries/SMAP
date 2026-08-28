classdef makeMovieTiff<interfaces.DialogProcessor
    %Links molecules in consecutive frames for SPT analysis
    methods
        function obj=makeMovieTiff(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
            obj.showresults=false;
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
% if isfield(obj.locData.pluginresults.Analyze.tracking,"Single_tracking_analysis")
%     singletracks=obj.locData.pluginresults.Analyze.tracking.Single_tracking_analysis;
% else
%     error("Run Single_trackin_analysis first")
% end
% 
% if isfield(obj.locData.pluginresults.Analyze.tracking,"Dual_tracking_analysis")
%     dualtracks=obj.locData.pluginresults.Analyze.tracking.Dual_tracking_analysis;
% else
%     dualtracks=[];
% end

%
%add to movie
dx=3;
dimmark=3;
if isempty(p.tiffile) 
    error("load tiff file first")
end

if ~isempty(p.Tfile)
    tt=load(p.Tfile).transformation;
elseif p.cdT
    tt=obj.locData.files.file.transformation;
end

il=imageloaderMM; il.attachPar(obj.P);
il.openi(p.tiffile);
il.prefit;  
emmirror=il.metadata.EMon;
% il.metadata.EMmirror=true;
numf=il.metadata.numberOfFrames;

if p.cdT
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
        if emmirror
             img=img(:,end:-1:1);
        end
        imgt=tt.transformImageToTarget(2,img,'pixel',il.metadata.roi);
        imrgb(:,:,1)=img;imrgb(:,:,2)=imgt;imrgb(:,:,3)=0;
        % figure(99); image(imrgb/max(imrgb(:)));
        if splitvert
            imcomb(:,:,1,k)=img(1:sx,:);
            % imcomb(:,:,3,k)=img(:,1:sx);
            imcomb(:,:,2,k)=imgt(1:sx,:);     
    
            % imcomb(:,:,1,k)=img(sx:end,:);
            % imcomb(:,:,2,k)=imgt(sx:end,:);  
        else
            imcomb(:,:,1,k)=img(:,1:sx);
            % imcomb(:,:,3,k)=img(:,1:sx);
            imcomb(:,:,2,k)=imgt(:,1:sx);    
        end
    end
else
    imcomb=zeros(il.metadata.Height,il.metadata.Width,3,numf);
     for k=1:numf
         img=double(il.getimage(k));
         if emmirror
             img=img(:,end:-1:1);
         end
         imcomb(:,:,1,k)=img(:,:);
     end
end

il.close;

% 
% layers=find(obj.getPar('sr_layerson'));
% [locs,indin]=obj.locData.getloc({'xnm','ynm','znm','xpix','ypix','frame','track_id','track_length','layer','filenumber','channel'},'layer',layers,'position','roi','grouping','ungrouped');
% 

imstack=imcomb;
if p.alternating
    intch=mean(mean(imstack-min(imstack(:)),1),2);
    n1=1:2:size(imstack,4);
    n2=2:2:size(imstack,4);
    in1=sum(squeeze(intch(1,1,1,n1))/sum(intch(1,1,1,:))+squeeze(intch(1,1,2,n2))/sum(intch(1,1,2,:)));
    in2=sum(squeeze(intch(1,1,1,n2))/sum(intch(1,1,1,:))+squeeze(intch(1,1,2,n1))/sum(intch(1,1,2,:)));
    if in1>in2
        indch1=n1; indch2=n2;
    else
        indch1=n2; indch2=n1;
    end
    imstacka=zeros(size(imstack,1),size(imstack,2),size(imstack,3),size(imstack,4)/2);
    imstacka(:,:,1,:)=imstack(:,:,1,indch1);
    imstacka(:,:,2,:)=imstack(:,:,2,indch2);
    % imstacka(:,:,3,:)=imstack(:,:,1,indch1)+0*imstack(:,:,2,indch2);
    % imstacka(:,:,dimmark,:)=(imstack(:,:,dimmark,indch1)+imstack(:,:,dimmark,indch2))/2;
    imstack=imstacka;
end

% fout=[fp filesep 'overlays' filesep 'combined_' fn ext];
% % fout=strrep(p.tiffile,'.ome.tif','_combined.tif');
% options.color=true;
% saveastiff(squeeze(single(imcomb)),fout,options);


% % if ~isempty(dualtracks)
%     switch p.showtraces.selection
%         case 'all'
%             dualtracks.trackstat.
%         case 'processive'
% 
% 
%         case 'processive co-tracks'
%         case 'none'
%     end
% else
% switch p.showtraces.selection
%     case 'all'
%         ids=singletracks.trackids;
%     case 'processive'
%         ids=singletracks.processiveids;
%     case 'processive co-tracks'
%         ids=find(dualtracks.trackstat.comovement);      
%     case 'none'
%         ids=[];
% end
% end

% exposuretime=obj.locData.files.file(locs.filenumber(1)).info.exposure;
% roi=obj.locData.files.file(locs.filenumber(1)).info.roi;
% pixelsize=obj.locData.files.file(locs.filenumber(1)).info.cam_pixelsize_um*1000;

% imstack=imcomb;
immax=max(imstack(:));
[lenx,leny]=size(imstack,[1,2]);

% for k=1:length(ids)
%     ind=locs.track_id==ids(k);
%     x=round(locs.xnm(ind)/pixelsize(1)-roi(1)); y=round(locs.ynm(ind)/pixelsize(2)-roi(2)); frame=locs.frame(ind);
%     if ~p.addboxes
%         dx=0;
%     end
%         x=min(lenx-dx,x);x=max(dx+1,x);y=min(leny-dx,y);y=max(dx+1,y);
%     for l=1:length(x)
% 
%         if p.addboxes
%             imstack(y(l)+dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
%             imstack(y(l)-dx, x(l)-dx:x(l)+dx,dimmark,frame(l))=immax;
%             imstack(y(l)-dx:y(l)+dx, x(l)+dx,dimmark,frame(l))=immax;
%             imstack(y(l)-dx:y(l)+dx, x(l)-dx,dimmark,frame(l))=immax;
%         end
% 
%         if p.addpartial
%             if l>1
%                 imstack(:,:,dimmark,frame(l))=insertLineImage(x(1:l),y(1:l),imstack(:,:,dimmark,frame(l)),immax);
%             end
%         end
%         if p.addtracks
%             imstack(:,:,dimmark,frame(l))=insertLineImage(x,y,imstack(:,:,dimmark,frame(l)),immax);
%         end
% 
%     end
% 
% end
 


[fp,fn,ext]=fileparts(p.tiffile);
if ~exist([fp filesep 'overlays'],"dir")
    mkdir([fp filesep 'overlays']);
end
fout2=[fp filesep 'overlays' filesep 'c_' fn ext];
% fout2=strrep(fout,'_combined.tif','_mark.tif');
options.color=true;
saveastiff(squeeze(single(imstack)),fout2,options);

end

             
function pard=guidef(obj)


% pard.showt.object=struct('String','Show:','Style','text');
% pard.showt.position=[1,1];
% pard.showt.Width=0.5;
% pard.showtraces.object=struct('String',{{'none','all','processive','processive co-tracks'}},'Style','popupmenu');
% pard.showtraces.position=[1,1.5];
% pard.showtraces.Width=1.5;
% 
% pard.addboxes.object=struct('String','boxes','Style','checkbox');
% pard.addboxes.position=[2,1];
% pard.addboxes.Width=1;
% 
% pard.addpartial.object=struct('String','partial tracks','Style','checkbox');
% pard.addpartial.position=[2,2];
% pard.addpartial.Width=1;
% 
% pard.addtracks.object=struct('String','full tracks','Style','checkbox');
% pard.addtracks.position=[2,3];
% pard.addtracks.Width=1;

pard.alternating.object=struct('String','alternating','Style','checkbox');
pard.alternating.position=[2,4];
pard.alternating.Width=1;

pard.tiffilet.object=struct('String','tif:','Style','text');
pard.tiffilet.position=[3,1];
pard.tiffilet.Width=0.5;

pard.tiffile.object=struct('String','','Style','edit');
pard.tiffile.position=[3,1.5];
pard.tiffile.Width=3.;

pard.tiffload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'tiffile','*.tif'}});
pard.tiffload.position=[3,4.5];
pard.tiffload.Width=0.5;

pard.cdT.object=struct('String','dual-color trafo:','Style','checkbox');
pard.cdT.position=[4,1];
pard.cdT.Width=1;

pard.Tfile.object=struct('String','','Style','edit');
pard.Tfile.position=[4,2];
pard.Tfile.Width=2.5;

pard.Tfload.object=struct('String','load','Style','pushbutton','Callback',{{@obj.loadtiff,'Tfile','*.mat'}});
pard.Tfload.position=[4,4.5];
pard.Tfload.Width=0.5;


pard.plugininfo.description=sprintf('co-tracking analysis');
pard.plugininfo.type='ProcessorPlugin';
end

