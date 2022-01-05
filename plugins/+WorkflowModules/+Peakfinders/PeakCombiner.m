classdef PeakCombiner<interfaces.WorkflowModule
%     Combines candidate positions from two separate channels into one, and
%     then transforms it back to the second channel. Used for global
%     fitting of multiple simultaneous channels.
    properties (Access=private)
        transform
    end
    methods
        function obj=PeakCombiner(varargin)
            obj@interfaces.WorkflowModule(varargin{:});
             obj.inputParameters={'loc_loc_filter_sigma','EMon','loc_ROIsize','loc_fileinfo','loc_multifile'};
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
            initGui@interfaces.WorkflowModule(obj);
        end
        function prerun(obj,p)
            l=load(p.Tfile);
            if isfield(l,'transformation')
                obj.transform=l.transformation;  
            elseif isfield(l,'SXY')
                obj.transform=l.SXY(1).cspline.global.transformation;
            elseif isfield(l,'saveloc')
                obj.transform=l.saveloc.file.transformation;
            elseif isfield(l,'T')
                obj.transform.T=double(cat(3,eye(3,3),permute(l.T,[3 ,2 ,1]))); % xXXX create transform with that matrix.
                obj.transform.centercoord=l.imgcenter;
            else
                errordlg(['no transformation found in' p.Tfile])
            end
            obj.setPar('loc_globaltransform',obj.transform);
        end
        function dato=run(obj,data,p)
            if isempty(data.data)
                dato=data;
                return;
            end
            roi=p.loc_fileinfo.roi;
            maxima=data.data;%get;
            transform=obj.transform;
            
            if isa(transform,'interfaces.LocTransform')
                if isfield(transform.tinfo,'cam_pixelsize_nm')
                    pixelsize=transform.tinfo.cam_pixelsize_nm;
                else
                    pixelsize=p.loc_fileinfo.cam_pixelsize_um*1000;
                end
                if isfield(transform.tinfo,'units')&&strcmp(transform.tinfo.units,'pixels')
                    pixelsize=[1 1];
                end
                xnm=(maxima.xpix+roi(1))*pixelsize(1);
                ynm=(maxima.xpix+roi(2))*pixelsize(end);
                indref=obj.transform.getRef(xnm,ynm);
                xnmref=xnm(indref);
                ynmref=ynm(indref);
                Nref=maxima.phot(indref).^2; %rather sqrt??? Think about this now!
                Nt=maxima.phot(~indref).^2;
                [xt,yt]=transform.transformCoordinatesInv(xnm(~indref),ynm(~indref));
                [iA,iB,uiA,uiB]=matchlocs(xnmref,ynmref,xt,yt,[0 0],pixelsize(1)*4);
                try
                xallr=vertcat(xnmref(uiA), xt(uiB'),((xnmref(iA).*Nref(iA)+xt(iB).*Nt(iB))./(Nref(iA)+Nt(iB))));
                yallr=vertcat(ynmref(uiA), yt(uiB'),((ynmref(iA).*Nref(iA)+yt(iB).*Nt(iB))./(Nref(iA)+Nt(iB))));
                catch err
                    disp('error')
                end

                xallrpix=round(xallr/pixelsize(1));
                yallrpix=round(yallr/pixelsize(end));

                %for challenge: allow for localizations close to rim. Correct
                %both reference and target if too far outside
                %now implemented only for right-left
                if contains(p.Tfile,'Challenge')
                    dn=(p.loc_ROIsize-1)/2;
                    xallrpix=min(max(xallrpix,dn+1),roi(3)/2-dn);
                    yallrpix=min(max(yallrpix,dn+1),roi(4)-dn);
                end


                [xallt,yallt]=transform.transformCoordinatesFwd(xallrpix*pixelsize(1),yallrpix*pixelsize(end));

                xalltpix=round(xallt/pixelsize(1));
                yalltpix=round(yallt/pixelsize(end)); 



                indgood=xalltpix>roi(1)&yalltpix>roi(2)&xalltpix<=roi(1)+roi(3)&yalltpix<=roi(2)+roi(4);

                dx=xallt(indgood)/pixelsize(1)-xalltpix(indgood);
                dy=yallt(indgood)/pixelsize(end)-yalltpix(indgood);

                xallrpix=xallrpix(indgood)-roi(1);
                yallrpix=yallrpix(indgood)-roi(2);
                xalltpix=xalltpix(indgood)-roi(1);
                yalltpix=yalltpix(indgood)-roi(2);


                xout=vertcat(xallrpix',xalltpix');
                yout=vertcat(yallrpix',yalltpix');
                dxout=vertcat(0*dx',dx');
                dyout=vertcat(0*dy',dy');
                idout=vertcat(1:sum(indgood),1:sum(indgood));

                maxout.xpix=xout(:);
                maxout.ypix=yout(:);

                maxout.ID=idout(:);
                maxout.dx=dxout(:);
                maxout.dy=dyout(:);
                dato=data;
                dato.data=maxout;
            elseif isa(transform,'interfaces.LocTransformN')
                if myisfield(transform,'unit')&&strcmp(transform.unit,'pixel')
                else
                    disp('calibration should be in pixel units for 4pi');
                    %set pixelsize in transform here
                end
                xpix=(maxima.xpix+roi(1)); %still x,y inconsistency! solve
                ypix=(maxima.ypix+roi(2)); %put on camera chip
                cpix=[xpix,ypix];
                indref=transform.getPart(1,cpix);
                
                cref=cpix(indref,:);
                Nall=sqrt(maxima.phot);
                Nc=sqrt(maxima.phot(indref));
%                 xf=cpix(indref,1);yf=cpix(indref,2);intf=maxima.intensity(indref);
                ccombined=cref;
                if p.loc_multifile
                    multioffsetpix=[p.loc_fileinfo.Width 0];%imagesize/2 find out if x or y
                    transform.info{2}.xrange=transform.info{2}.xrange+p.loc_fileinfo.Width;
                else
                    multioffsetpix=0;
                end
                for k=2:transform.channels
                    indch=transform.getPart(k,cpix);
                    Nt=Nall(indch);
                    if p.framecorrection
                        ctarget=transform.transformToReferenceFramecorrection(k,cpix(indch,:)-multioffsetpix,data.frame);
                    else
                        ctarget=transform.transformToReference(k,cpix(indch,:)-multioffsetpix);
                    end
                    
                    [iA,iB,uiA,uiB]=matchlocs(ccombined(:,1),ccombined(:,2),ctarget(:,1),ctarget(:,2),[0 0],6);
                    if isempty(iA)
                        cnew=[];
                    else
                        cnew=(ccombined(iA,:).*Nc(iA)+ctarget(iB,:).*Nt(iB))./(Nc(iA)+Nt(iB));
                    end
                    ccombined=vertcat(ccombined(uiA,:),ctarget(uiB,:),cnew);
                    Nc=vertcat(Nc(uiA),Nt(uiB),(Nc(iA)+Nt(iB)));
                    
                end
            
%                 cr=round(ccombined);
                %offset between true and rounded position in channel 1
              
                if p.framecorrection
                    ct=transform.transformToTargetAllFramecorrection(ccombined,data.frame);
                else
                    ct=transform.transformToTargetAll(ccombined); %transfrom not rounded
                end
                ct(:,:,2)=ct(:,:,2)+multioffsetpix;
                ct(:,1,:)=ct(:,1,:)-roi(1); %bring back to ROI on camera
                ct(:,2,:)=ct(:,2,:)-roi(2);
                %test dc XXXXX
                ctt=ct;
            
                ctr=round(ctt);
                dc=ct-ctr;


                cout=permute(ctr,[2 3 1]);
                dcout=permute(dc,[2 3 1]);
                indout=repmat((1:size(ctr,1)),1,transform.channels);
                xh=cout(1,:,:);yh=cout(2,:,:);
                dxh=dcout(1,:,:);dyh=dcout(2,:,:);
                
                maxout.xpix=squeeze(xh(:));
                maxout.ypix=squeeze(yh(:));

                maxout.ID=indout(:);
                maxout.dx=squeeze(dxh(:));
                maxout.dy=squeeze(dyh(:));
                dato=data;
                dato.data=maxout;
            elseif isfield(obj.transform,'T') %new 4Pi

                xpix=(maxima(1).xpix+roi(1)); %still x,y inconsistency! solve
                ypix=(maxima(1).ypix+roi(2)); %put on camera chip
                cref=[xpix,ypix,ones(size(xpix))];
                
                Nc=(maxima(1).phot);
                ccombined=cref;
%                 ct(:,:,1)=cref(:,1:2);
                for k=2:length(maxima)    
                    Nt=maxima(k).phot;
                    T=obj.transform.T(:,:,k);
                    Tinv=inv(T);
                    cN=[maxima(k).xpix+roi(1),maxima(k).ypix+roi(2),ones(size(maxima(k).ypix))];
                    ctarget=transformT4Pi(cN,Tinv,obj.transform.centercoord);
%                     ctarget=(Tinv*cN')';

%                     if p.framecorrection
%                         ctarget=transform.transformToReferenceFramecorrection(k,cpix(indch,:)-multioffsetpix,data.frame);
%                     else
%                         ctarget=transform.transformToReference(k,cpix(indch,:)-multioffsetpix);
%                     end
                    
                    [iA,iB,uiA,uiB]=matchlocs(ccombined(:,1),ccombined(:,2),ctarget(:,1),ctarget(:,2),[0 0],6);
                    if isempty(iA)
                        cnew=[];
                    else
                        cnew=(ccombined(iA,:).*Nc(iA)+ctarget(iB,:).*Nt(iB))./(Nc(iA)+Nt(iB));
                    end
                    ccombined=vertcat(ccombined(uiA,:),ctarget(uiB,:),cnew);
                    Nc=vertcat(Nc(uiA),Nt(uiB),(Nc(iA)+Nt(iB)));
%                     ctar=(T*ccombined')';
%                     ct(end+1:end+size(ctar,1),:,k)=ctar(:,1:2);
                end

                for k=length(maxima):-1:1
                    T=obj.transform.T(:,:,k);
                    ct=transformT4Pi(ccombined,T,obj.transform.centercoord);
%                     ct=(T*ccombined')';
                    ct(:,1)=ct(:,1)-roi(1); %bring back to ROI on camera
                    ct(:,2)=ct(:,2)-roi(2);
                    ctt=ct(:,1:2);
                    ctr=round(ctt);
                    dc=ct(:,1:2)-ctr;

                    maxout(k).xpix=ct(:,1);
                    maxout(k).ypix=ct(:,2);
    
                    maxout(k).ID=k*ones(size(maxout(k).xpix));
                    maxout(k).dx=dc(:,1);
                    maxout(k).dy=dc(:,2);
                end
           

% %                 ct(:,:,2)=ct(:,:,2);
%                 ct(:,1,:)=ct(:,1,:)-roi(1); %bring back to ROI on camera
%                 ct(:,2,:)=ct(:,2,:)-roi(2);
%                 %test dc XXXXX
%                 ctt=ct;
%             
%                 ctr=round(ctt);
%                 dc=ct-ctr;
% 
% 
%                 cout=permute(ctr,[2 3 1]);
%                 dcout=permute(dc,[2 3 1]);
%                 indout=repmat((1:size(ctr,1)),1,length(maxima));
%                 xh=cout(1,:,:);yh=cout(2,:,:);
%                 dxh=dcout(1,:,:);dyh=dcout(2,:,:);
                
    
                dato=data;
                dato.data=maxout;

            else

                adslf
            end
        end
        function loadbutton(obj,a,b)
            fn=obj.guihandles.Tfile.String;
            path=fileparts(fn);
            if ~exist(path,'file')
                fn=[fileparts(obj.getPar('loc_outputfilename')) filesep '*.mat'];
            end
            [f,path]=uigetfile(fn,'Select transformation file _T.mat');
            if f
                Tload=load([path f],'transformation');
                if ~isfield(Tload,'transformation')
                    msgbox('could not find transformation in file. Load other file?')
                end
                
                obj.guihandles.Tfile.String=[path f];
                obj.setPar('transformationfile',[path f]);
            end      
        end  
    end
end

function ct=transformT4Pi(ccombined,T,centercoord)
ct=(T*(ccombined-centercoord)')'+centercoord;
end


function pard=guidef(obj)
pard.Tfile.object=struct('Style','edit','String','*_T.mat');
pard.Tfile.position=[1,1];
pard.Tfile.Width=3;

pard.loadbutton.object=struct('Style','pushbutton','String','load T','Callback',@obj.loadbutton);
pard.loadbutton.position=[1,4];

pard.framecorrection.object=struct('Style','checkbox','String','frame dependent transformation');
pard.framecorrection.position=[2,3];
pard.framecorrection.Width=2;


pard.syncParameters={{'transformationfile','Tfile',{'String'}}};

pard.plugininfo.type='WorkflowModule'; 
pard.plugininfo.description='Combines candidate positions from two separate channels into one, and then transforms it back to the second channel. Used for global fitting of multiple simultaneous channels.';
end