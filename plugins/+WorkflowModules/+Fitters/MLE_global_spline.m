classdef MLE_global_spline<interfaces.WorkflowFitter
%     Global fitter for multiple channels.
    properties
        fitpar
    end
    methods
        function obj=MLE_global_spline(varargin)
            obj@interfaces.WorkflowFitter(varargin{:})
            obj.inputChannels=1; 
             obj.setInputChannels(1,'frame');
             
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function fitinit(obj)
            obj.infofields={'x','y','ID','dx','dy'};
            obj.fitpar=getfitpar(obj);
            switch obj.fitpar.mode
                case '4pi'
                    fitterpath=[fileparts(obj.getPar('maindirectory')) filesep 'ries-private' filesep 'PSF4Pi'];
                    if ~isdeployed
                        addpath(fitterpath)
                    end
                    obj.fitpar.link=obj.fitpar.link([2 1 4 5 3 6]);
                    obj.fitpar.fitfunction=@mleFit_LM_4Pi;
                case 'Gauss'
                     obj.fitpar.fitfunction=@mleFit_LM_global_gauss;
                     disp('only implemented for symmetric Gauss: fittype=2');
                otherwise
                    obj.fitpar.fitfunction=@mleFit_LM_global; %later: include single channel, decide here
%                     GPUmleFit_LM_MultiChannel_Gauss
            end
             transform=obj.getPar('loc_globaltransform');
             
             if isa(transform,'interfaces.LocTransform')
                 obj.fitpar.mirrorud=contains(transform.tinfo.mirror.targetmirror,'up');
                 obj.fitpar.mirrorrl=contains(transform.tinfo.mirror.targetmirror,'right');
             elseif isa(transform,'interfaces.LocTransformN')
                 tinfo=transform.info;
                 for k=1:length(tinfo)
                     obj.fitpar.mirror{k}=tinfo{k}.mirror;
                 end
             else
                 disp('transform not recognized');
             end
            roisize=obj.getPar('loc_ROIsize');
            obj.numberInBlock=round(obj.fitpar.roisperfit*100/roisize^2/12)*12;
            
            if obj.fitpar.fitmode==5 ||obj.fitpar.fitmode==6              
                obj.fitpar.mirrorstack=false; %later: remove completely
                p=obj.getAllParameters;
                if p.overwritePixelsize
                    obj.setPar('overwrite_pixelsize',[p.pixelsizex p.pixelsizey])
                    cs=obj.getPar('loc_cameraSettings');
                    cs.cam_pixelsize_um=[p.pixelsizex p.pixelsizey];
                    obj.setPar('loc_cameraSettings',cs);
                else
                    obj.setPar('overwrite_pixelsize',[])
                end
                 obj.setPar('loc_iterations',p.iterations);
            end   
            obj.setPar('loc_numberOfChannels',2);
        end
        function nofound(obj,varargin)
            disp('fit function not working. Wrong Cuda version?')
        end

        function locs=fit(obj,imstack,stackinfo)
            if obj.fitpar.fitmode==3
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.zparhere=[obj.fitpar.zpar{X,Y}(:)];
            elseif obj.fitpar.fitmode==5 || obj.fitpar.fitmode==6
                X=stackinfo.X;Y=stackinfo.Y;
                obj.fitpar.splinefithere=[obj.fitpar.splinefit{X,Y}(:)];
            end
            if obj.fitpar.issCMOS
                varstack=getvarmap(obj.fitpar.varmap,stackinfo,size(imstack,1));
            else
                varstack=0;
            end
            switch obj.fitpar.mode
                case '4pi'
                     out=fitwrapper_4pi(imstack,obj.fitpar,stackinfo,varstack);
                     locs=fit2locs_4pi(out,stackinfo,obj.fitpar,imstack);
                otherwise
            out=fitwrapper_global(imstack,obj.fitpar,stackinfo,varstack);
            locs=fit2locs_global(out,stackinfo,obj.fitpar,imstack);
            end
        end

        function initGui(obj)
            initGui@interfaces.WorkflowFitter(obj);
            obj.guihandles.loadcal.Callback={@loadcall_callback,obj};

            %make global fit control
            hold=obj.guihandles.globaltable;
            gt=uitable(hold.Parent,'Position',hold.Position);
            
            gt.RowName={'x','y','z','N','Bg','ph'};
            gt.ColumnName={'l','x'};
            
            hh=obj.guiPar.fontsize+5;
            wh=hold.Position(3)-hh-50;
            gt.ColumnWidth={hh,wh};
            gt.Data={true,'1';true,'1';true,'1';false,'1';false,'1';true,'1'};
            gt.ColumnEditable=true;
            obj.guihandles.globaltable=gt;
            hold.Visible='off';
            hold.delete;

        end
            
    end
end

function loadscmos_callback(a,b,obj)
fs=obj.getSingleGuiParameter('scmosfile');
if isempty(fs)
    fs='*.*';
end
[file,pfad]=uigetfile(fs);
if file
    obj.setGuiParameters(struct('scmosfile',[pfad file]))
end

end


function locs=fit2locs_global(results,stackinfo,fitpar,image)
if isempty(results)
    locs=[];
    return
end
numl=size(results.P,1);
numpar=5;
numchannels=2;

v1=ones(numl,1,'single');
s=size(image);          
dn=ceil((s(1)-1)/2)*v1;

shiftx=0;%-0.5; %deviation from ground truth
shifty=0;%-0.5;
posx=stackinfo.xpix(results.indused)+shiftx;
posy=stackinfo.ypix(results.indused)+shifty;
frame=stackinfo.frame(results.indused);
P=results.P;
EMexcess=fitpar.EMexcessNoise;
CRLB=results.CRLB;
LogL=results.LogL;
           CRLB(isnan(CRLB))= 0; %XXXXXXXXX
           LogL(isnan(LogL))= 0; %XXXXXXXXX
           CRLB((CRLB)<0)= 0; %XXXXXXXXX
          % LogL((LogL)<0)= 0; %XXXXXXXXX
           
o=ones( numpar,1);fac=num2cell(o);z=zeros( numpar,1);off=num2cell(z);
faccrlb=fac;
% locs.xpix=P(:,2)-dn+posx;
if (fitpar.fitmode==5||fitpar.fitmode==6) && fitpar.mirrorstack
    fac{2}=-1;
    off{2}=dn+posx;
%     locs.xpix=dn-P(:,2)+posx;
else
    fac{2}=1;
    off{2}=-dn+posx;
%     locs.xpix=P(:,2)-dn+posx;
end
% locs.ypix=P(:,1)-dn+posy;
off{1}=-dn+posy;

switch fitpar.mode
    case 'Gauss'
        fac{3}=EMexcess;
        fac{4}=EMexcess;
        faccrlb{3}=EMexcess;
        faccrlb{4}=EMexcess;
        names={'ypix','xpix','phot','bg','PSFxpix'};
        namesav={'ypix','xpix','PSFxpix','phot'};
    
    case {'Spline','cspline'}
        fac{4}=EMexcess*fitpar.splinefithere.normf;
        fac{5}=EMexcess;
        faccrlb{4}=EMexcess;
        faccrlb{5}=EMexcess;
        fac{3}=-fitpar.dz*fitpar.refractive_index_mismatch; % negative factor to fix direction in SMAP
        off{3}=fitpar.z0*fitpar.dz*fitpar.refractive_index_mismatch;
        faccrlb{3}=fitpar.dz*fitpar.refractive_index_mismatch;
        names={'ypix','xpix','znm','phot','bg'};
        namesav={'ypix','xpix','znm','phot'};
        namesphot={'bg','phot'};
end

% locs.phot=P(:,4)*EMexcess;
% locs.bg=P(:,5)*EMexcess;

locs.frame=frame;

locs.logLikelihood=LogL;%/sum(fitpar.weightsch);

locs.peakfindx=posx;
locs.peakfindy=posy;
   
         sx=1*v1;
        locs.PSFxpix=sx;
        locs.PSFypix=sx;
% end

% fac and faccrlb should be of length(channels). Then use for each channel.
% If linked, sum them up??? no. Careful if to sum up the inverses.

%take care of weighting
locs.logLikelihood=locs.logLikelihood/sum(fitpar.weightsch);


linked=fitpar.link;
ind=1;
for k=1:length(names)
    if any(contains(namesphot,names{k})) 
        errfactor=1./sqrt(fitpar.weightsch);
        valfactor=1./fitpar.weightsch;
        valfactorlinked=min(valfactor);
        errfactorlinked=min(errfactor);
    else
        errfactor=sqrt(fitpar.weightsch);
        valfactor=ones(size(fitpar.weightsch));  
%         errfactorlinked=max(errfactor);
        errfactorlinked=(sqrt(mean(errfactor.^2)));
        valfactorlinked=1;
    end
    
    if linked(k)
        locs.(names{k})=(P(:,ind).*fac{k}(1)+off{k})*valfactorlinked;
        locs.([names{k} 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k}*errfactorlinked;
        ind=ind+1;
    else
        v=zeros(numl,numchannels,'single');
        ve=zeros(numl,numchannels,'single');
        for c=1:numchannels
            ch=num2str(c);
            fach=fac{k}(min(length(fac{k}),c)); %channel-dependent factor, e.g. for normalization of PSF
            locs.([names{k} ch])=(P(:,ind).*fach+off{k})*valfactor(c);
            locs.([names{k} ch 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k}*errfactor(c);
            v(:,c)=locs.([names{k} ch]);
            ve(:,c)=locs.([names{k} ch 'err']);
%             if any(contains(namesphot,names{k}))  % correct for intensities
%                 v(:,c)=v(:,c)/fitpar.weightsch(c);
%                 locs.([names{k} ch])=locs.([names{k} ch])/fitpar.weightsch(c);
%             end
            ind=ind+1;
        end
        if any(contains(namesav,names{k}))
            switch fitpar.mainchannel
                case 1
                    locs.([names{k} 'err'])=1./sqrt(1./ve(:,1).^2+1./ve(:,2).^2)*sqrt(2);
                    ve(:,1)=ve(:,1)/errfactor(1);
                    ve(:,2)=ve(:,2)/errfactor(2);
                    locs.(names{k})=sum(v./ve,2)./sum(1./ve,2); %average based on CRLB before weighting
                case 2
                    locs.(names{k})=v(:,1);
                    locs.([names{k} 'err'])=ve(:,1);
                case 3
                    locs.(names{k})=v(:,2);
                    locs.([names{k} 'err'])=ve(:,2);
            end
        end
    end
end
locs.iterations=P(:,end);
global testloc
testloc=locs;
end

% function errout=rescaleerr(errin,weight)
% errout=errin/sqrt(weight);
% end

function locs=fit2locs_4pi(results,stackinfo,fitpar,image)
if isempty(results)
    locs=[];
    return
end
numl=size(results.P,1);
numpar=6;
numchannels=4;

v1=ones(numl,1,'single');
s=size(image);          
dn=ceil((s(1)-1)/2)*v1;

shiftx=0;%-0.5; %deviation from ground truth
shifty=0;%-0.5;
posx=stackinfo.xpix(results.indused)+shiftx;
posy=stackinfo.ypix(results.indused)+shifty;
frame=stackinfo.frame(results.indused);
P=results.P;
EMexcess=fitpar.EMexcessNoise;
CRLB=results.CRLB;
LogL=results.LogL;
           CRLB(isnan(CRLB))= 0; %XXXXXXXXX
           LogL(isnan(LogL))= 0; %XXXXXXXXX
           CRLB((CRLB)<0)= 0; %XXXXXXXXX
          % LogL((LogL)<0)= 0; %XXXXXXXXX
           
o=ones( numpar,1);fac=num2cell(o);z=zeros( numpar,1);off=num2cell(z);
faccrlb=fac;
if (fitpar.fitmode==5||fitpar.fitmode==6) && fitpar.mirrorstack
    fac{2}=-1;
    off{2}=dn+posx;
else
    fac{2}=1;
    off{2}=-dn+posx;
end
off{1}=-dn+posy;

fac{4}=EMexcess;
fac{3}=EMexcess;
faccrlb{4}=EMexcess;
faccrlb{3}=EMexcess;
locs.frame=frame;

locs.logLikelihood=LogL;

locs.peakfindx=posx;
locs.peakfindy=posy;
sx=1*v1; 
locs.PSFxpix=sx;
locs.PSFypix=sx;

fac{5}=-fitpar.dz*fitpar.refractive_index_mismatch; %fix direction
off{5}=-fitpar.z0*fitpar.dz*fitpar.refractive_index_mismatch;
faccrlb{5}=fitpar.dz*fitpar.refractive_index_mismatch;

names={'ypix','xpix','phot','bg','zastig','phase'};
namesav={'ypix','xpix','zastig','phot'};
linked=fitpar.link(1:6);

ind=1;
for k=1:length(names)
    if linked(k)
        locs.(names{k})=P(:,ind).*fac{k}+off{k};
        locs.([names{k} 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k};
        ind=ind+1;
    else
        v=zeros(numl,numchannels,'single');
        ve=zeros(numl,numchannels,'single');
        for c=1:numchannels
            ch=num2str(c);
            locs.([names{k} ch])=P(:,ind).*fac{k}+off{k};
            if strcmp(names{k},'phase')
                locs.([names{k} ch])=mod(locs.([names{k} ch]),2*pi);
            end
            locs.([names{k} ch 'err'])=sqrt(CRLB(:,ind)).*faccrlb{k};
            v(:,c)=locs.([names{k} ch]);
            ve(:,c)=locs.([names{k} ch 'err']);
            ind=ind+1;
        end
        if any(contains(namesav,names{k}))
            switch fitpar.mainchannel
                case 1
                    locs.(names{k})=sum(v./ve,2)./sum(1./ve,2);
                    locs.([names{k} 'err'])=1./sqrt(1./ve(:,1).^2+1./ve(:,2).^2);
                case 2
                    locs.(names{k})=v(:,1);
                    locs.([names{k} 'err'])=ve(:,1);
                case 3
                    locs.(names{k})=v(:,2);
                    locs.([names{k} 'err'])=ve(:,2);
            end

        end
    end
end
locs.iterations=P(:,end);
end

function out=fitwrapper_global(imstack,fitpar,stackinfo,varstack)
numberOfChannels=2;
nfits=ceil(size(imstack,3)/numberOfChannels);
npar=5;


s=size(imstack);
if length(s)==2 
 s(3)=1;
end
if s(3)<numberOfChannels  %sorting: needs at least two 
    out=[];
 return
end
% fitpar=obj.fitpar;
EMexcess=fitpar.EMexcessNoise;
if isempty(EMexcess)
    EMexcess=1;
end
         
dT=zeros(npar,2,(nfits));
dT(1,1,:)=stackinfo.dy(1:numberOfChannels:end); %for nonrounding
dT(2,1,:)=stackinfo.dx(1:numberOfChannels:end);

dT(1,2,:)=stackinfo.dy(2:numberOfChannels:end);
dT(2,2,:)=stackinfo.dx(2:numberOfChannels:end);

disp('check MLE_globalspline')

imfit(:,:,:,1)=imstack(:,:,1:numberOfChannels:end);
if isfield(fitpar,'mirrorud') && fitpar.mirrorud
    imfit(:,:,:,2)=imstack(end:-1:1,:,2:numberOfChannels:end);
    dT(1,2,:)=-dT(1,2,:);
elseif isfield(fitpar,'mirror')  %now only mirror channel two!
    mirr=fitpar.mirror{2};
    switch mirr
        case 0 %no mirror
            imfit(:,:,:,2)=imstack(:,:,2:numberOfChannels:end);
        case 1 %righ-left mirror
            imfit(:,:,:,2)=imstack(:,end:-1:1,2:numberOfChannels:end);
            dT(2,2,:)=-dT(2,2,:);
        case 2 %up-down mirror
            imfit(:,:,:,2)=imstack(end:-1:1,:,2:numberOfChannels:end);
            dT(1,2,:)=-dT(1,2,:);
    end

    
else
    imfit(:,:,:,2)=imstack(:,:,2:numberOfChannels:end);
end

if fitpar.weightsch(1)~=1
    imfit(:,:,:,1)=imfit(:,:,:,1)*fitpar.weightsch(1);
end
if fitpar.weightsch(2)~=1
    imfit(:,:,:,2)=imfit(:,:,:,2)*fitpar.weightsch(2);
end

numframes=size(imfit,3); 
sharedA = repmat(int32(fitpar.link(1:5)'),[1 numframes]);
out.indused=1:numberOfChannels:numframes*numberOfChannels;   %XXXXXXX check. Wrong? results shoud be only displayed in one channel

switch fitpar.mode
    case 'Gauss'
        arguments{1}=imfit/EMexcess;
        arguments{2}=uint32(2);
        arguments{3}=uint32(sharedA);
        arguments{4}=uint32(fitpar.iterations);
        arguments{5}=single(fitpar.PSFx0);
        arguments{6}=single(dT); %XXXXXXXXXX
        %imstack, sharedflag, iterations, spline coefficients, channelshift,
        %fitmode, varmap
%         arguments{6}=fitpar.fitmode;
%         arguments{7}=fitpar.zstart/fitpar.dz;
    case {'Spline','cspline'}
        arguments{1}=imfit/EMexcess;
        arguments{2}=sharedA;
        arguments{3}=fitpar.iterations;
        arguments{4}=single(fitpar.splinefithere.coeff);
        arguments{5}=single(dT); %XXXXXXXXXX
        %imstack, sharedflag, iterations, spline coefficients, channelshift,
        %fitmode, varmap
%         arguments{6}=fitpar.fitmode;
        arguments{6}=fitpar.zstart/fitpar.dz;
    otherwise
        disp('fitmode not implemented for global fitting')
end

[P CRLB LogL]=fitpar.fitfunction(arguments{:});


%subtract dT for y
if fitpar.link(1)
    P(:,1)=P(:,1)+squeeze(dT(1,1,:));
    ind2=1;
else
    P(:,1)=P(:,1)+squeeze(dT(1,1,:));
    P(:,2)=P(:,2)+squeeze(dT(1,1,:));
    ind2=2;
end

%subtract dT for x
if fitpar.link(2)
    P(:,ind2+1)=P(:,ind2+1)+squeeze(dT(2,1,:));
else
    P(:,ind2+1)=P(:,ind2+1)+squeeze(dT(2,1,:));
    P(:,ind2+2)=P(:,ind2+2)+squeeze(dT(2,1,:));
end

% [PG,CRLBG, LLG] =  GPUmleFit_LM_MultiChannel_Gauss(d_data,2 (fitmode),uint32(sharedA),iterations,single(PSFstart),single(dT));


% [P, CRLB,LogL, res]=CPUmleFit_LM_MultiChannel_R(arguments{:});
% htot=P(:,8)
% 
% arguments{5}=varstack;
% arguments{6}=1;
% 
%     switch fitpar.fitmode
%         case {1,2,4} %fix
%             arguments{4}=fitpar.PSFx0;
%             arguments{1}=single(imstack/EMexcess);
% %         case 2 %free
%         case 3 %z
%             arguments{1}=single(imstack/EMexcess);
%             arguments{4}=single(fitpar.zparhere);
% %         case 4 %sx sy
%         case {5,6} %spline   
%             if fitpar.mirrorstack
%                 arguments{1}=single(imstack(:,end:-1:1,:)/EMexcess);
%             else
%                 arguments{1}=single(imstack/EMexcess);
%             end
%             arguments{4}=single(fitpar.splinefithere.cspline.coeff);
%     end

out.P=P;
out.CRLB=CRLB;
 out.LogL=LogL;
end
 
function out=fitwrapper_4pi(imstack,fitpar,stackinfo,varstack)
numberOfChannels=size(fitpar.splinefithere.coeff.phaseshifts,2);
nfits=ceil(size(imstack,3)/numberOfChannels);
npar=6;
numframes=size(imstack,3); 

s=size(imstack);
if length(s)==2 
 s(3)=1;
end
if s(3)<numberOfChannels  %sorting: needs at least two 
    out=[];
 return
end
EMexcess=fitpar.EMexcessNoise;
if isempty(EMexcess)
    EMexcess=1;
end

if isfield(fitpar.splinefithere.coeff,'normf')
    normf=fitpar.splinefithere.coeff.normf;
else
    normf=ones(4,1);
end
        
dT=zeros(npar,numberOfChannels,(nfits),'single');
imfit=zeros(s(1),s(2),ceil(s(3)/numberOfChannels),numberOfChannels,'single');
for k=1:numberOfChannels
    dT(1,k,:)=stackinfo.dy(k:numberOfChannels:end);
    dT(2,k,:)=stackinfo.dx(k:numberOfChannels:end);
    if fitpar.mirror{k}==1 || fitpar.mirror{k}==3  %mirror x or xy
        imfit(:,:,:,k)=imstack(end:-1:1,:,k:numberOfChannels:end)/normf(k);
        dT(1,k,:)=-dT(1,k,:);
    else
        imfit(:,:,:,k)=imstack(:,:,k:numberOfChannels:end)/normf(k);
    end
    
end
sharedA = int32(fitpar.link);
sharedA(end+1)=1; %link phase. Not part of GUI

out.indused=1:numberOfChannels:numframes;

% [P,CRLB1 LL] = CPUmleFit_LM_MultiChannel(single(imstack(:, :, :, :)),uint32(shared),iterations,single(PSF.Ispline), single(PSF.Aspline),single(PSF.Bspline),single(dTAll),single(phi0),z0);
% imageslicer(residual)

arguments{1}=single(imfit/EMexcess);
arguments{2}=uint32(sharedA);
arguments{3}=fitpar.iterations;
arguments{4}=single(fitpar.splinefithere.coeff.Ispline);
arguments{5}=single(fitpar.splinefithere.coeff.Aspline);
arguments{6}=single(fitpar.splinefithere.coeff.Bspline);
arguments{7}=single(dT); %XXXXXXXXXX
%imstack, sharedflag, iterations, spline coefficients, channelshift,
%fitmode, varmap
arguments{8}=single(fitpar.splinefithere.coeff.phaseshifts);
arguments{9}=fitpar.zstart/fitpar.splinefithere.dz+fitpar.splinefithere.z0;
[P CRLB LogL]=fitpar.fitfunction(arguments{:});

out.P=P;
out.CRLB=CRLB;
 out.LogL=LogL;
end
        

function out=fitwrapper(imstack,fitpar,varstack)
s=size(imstack);
if length(s)==2 
 s(3)=1;
end
if s(3)==0
    out=[];
 return
end
% fitpar=obj.fitpar;
EMexcess=fitpar.EMexcessNoise;
if isempty(EMexcess)
    EMexcess=1;
end

arguments{2}=fitpar.fitmode;
arguments{3}=fitpar.iterations;

arguments{5}=varstack;
arguments{6}=1;

    switch fitpar.fitmode
        case {1,2,4} %fix
            arguments{4}=fitpar.PSFx0;
            arguments{1}=single(imstack/EMexcess);
%         case 2 %free
        case 3 %z
            arguments{1}=single(imstack/EMexcess);
            arguments{4}=single(fitpar.zparhere);
%         case 4 %sx sy
        case {5,6} %spline   
            if fitpar.mirrorstack
                arguments{1}=single(imstack(:,end:-1:1,:)/EMexcess);
            else
                arguments{1}=single(imstack/EMexcess);
            end
            arguments{4}=single(fitpar.splinefithere.cspline.coeff);
    end
   
    if fitpar.fitmode==6
          
        [P1 CRLB1 LL1 P2 CRLB2 LL2 ]=fitpar.fitfunction(arguments{:});
        P = repmat(single(LL1>=LL2),1,6).*P1+repmat(single(LL1<LL2),1,6).*P2;
        CRLB = repmat(single(LL1>=LL2),1,5).*CRLB1+repmat(single(LL1<LL2),1,5).*CRLB2;
        LogL = repmat(single(LL1>=LL2),1,1).*LL1+repmat(single(LL1<LL2),1,1).*LL2;
    else
        [P CRLB LogL]=fitpar.fitfunction(arguments{:});
    end

out.P=P;
out.CRLB=CRLB;
out.LogL=LogL;
end
        
        
function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
ph=fileparts(p.cal_3Dfile);
if ~exist(ph,'file')
    p.cal_3Dfile= [fileparts(obj.getPar('loc_outputfilename')) filesep '*.mat'];
end
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY') && ~isfield(l,'cspline')
        msgbox('no 3D data recognized. Select other file.');
    end
    obj.setGuiParameters(struct('cal_3Dfile',[p f]));
     obj.setPar('cal_3Dfile',[p f]);
    if isfield(l,'transformation')
        obj.setPar('transformationfile',[p f]);
    end
       
end
end

function fitpar=getfitpar(obj)
p=obj.getAllParameters;
fitpar.iterations=p.iterations;
fitpar.fitmode=p.fitmode.Value;
fitpar.roisperfit=p.roisperfit;
fitpar.issCMOS=false;
fitpar.mainchannel=p.mainchannel.Value;
fitpar.weightsch=p.weightsch;
if length(fitpar.weightsch)==1
    fitpar.weightsch(2)=1;
end
for k=1:size(p.globaltable.Data,1)
    fitpar.factor{k}=str2num(p.globaltable.Data{k,2});
end
fitpar.link=[p.globaltable.Data{:,1}];
fitpar.zstart=p.zstart;
if fitpar.fitmode==3||fitpar.fitmode==5
    fitpar.issCMOS=p.isscmos;
%     fitpar.PSF2D=p.fit2D;
%     if p.fit2D
%         fitpar.fitmode=6;
%     end
    
    calfile=p.cal_3Dfile;
    cal=load(calfile);

    fitpar.objPos=0;
    if isfield(cal,'outforfit')
        fitpar.zpar{1,1}=cal.outforfit;
        fitpar.mode='zcal';
    elseif isfield(cal,'SXY')
        if p.isglobal 
            SS=cal.SXY_g;
        else
            SS=cal.SXY;
        end
        s=size(SS);
        Z=1;
%         if p.useObjPos
%             zr=cal.SXY(1).Zrangeall;
%             zr(1)=[];zr(end)=inf;
%             Z=find(p.objPos<=zr,1,'first');
%         end
        for X=s(1):-1:1
            for Y=s(2):-1:1
%                     if isfield(cal.SXY(X,Y,Z),'gauss_zfit')
                zpar{X,Y}=SS(X,Y,Z).gauss_zfit;
%                     else
%                         zpar{X,Y}=[];
%                     end
                %global:combine splines
%                 cs=cal.SXY(X,Y,Z).cspline_all;
                cs=SS(X,Y,Z).cspline;
                if iscell(cs.coeff)
                    coeff(:,:,:,:,1)=cs.coeff{1};
                    coeff(:,:,:,:,2)=cs.coeff{2};
                    cs.coeff=single(coeff);
                end
                splinefit{X,Y}=cs;
                obj.spatialXrange{X,Y}=SS(X,Y).Xrange;
                obj.spatialYrange{X,Y}=SS(X,Y).Yrange;
                

            end
        end
        if ~isempty(splinefit{1})
            fitpar.dz=splinefit{1}.dz;
            fitpar.z0=splinefit{1}.z0;
            fitpar.splinefit=splinefit;
        end
        fitpar.zpar=zpar;

        if numel(SS)>1
            obj.spatial3Dcal=true;
        else
            obj.spatial3Dcal=false;
        end
%         xr=SS(1,1).Xrangeall;
%         xr(1)=-inf;xr(end)=inf;
%         yr=SS(1,1).Yrangeall;
%         yr(1)=-inf;yr(end)=inf;
%         obj.spatialXrange=xr;
%         obj.spatialYrange=yr;
        fitpar.EMon=SS(1).EMon;
        fitpar.mode='cspline';
    elseif isfield(cal,'cspline')
        fitpar.zpar{1}=cal.gauss_zfit;
        fitpar.dz=cal.cspline.dz;
        fitpar.z0=cal.cspline.z0;
        fitpar.splinefit{1}=cal.cspline_all;
        if ~isfield(fitpar.splinefit{1}.cspline,'isEM')
            fitpar.splinefit{1}.cspline.isEM=false;
        end
        fitpar.EMon=false;
        fitpar.mode='csplineold';
    elseif isfield(cal,'cal4pi')
        fitpar.mode='4pi';
        fitpar.splinefit{1}=cal.cal4pi;
        fitpar.EMon=cal.EMon;
        fitpar.dz=cal.cal4pi.dz;
        fitpar.z0=cal.cal4pi.z0;
        
        savefit=copyfields([],cal.cal4pi,{'dz','z0','x0','transformation'});
        savefit=copyfields(savefit,cal.cal4pi.coeff,{'frequency','phaseshifts'});
        savefit.cal_3Dfile=p.cal_3Dfile;
        obj.setPar('savefit',struct('cal3D',savefit));
    else
        disp('no calibration found')

    end

    %load sCMOS
    if p.isscmos
        [~,~,ext]=fileparts(p.scmosfile);
        switch ext
            case '.tif'
                varmaph=imread(p.scmosfile);
            case '.mat'
                varmaph=load(p.scmosfile);
                if isstruct(varmaph)
                    fn=fieldnames(varmaph);
                    varmaph=varmaph.(fn{1});
                end
            otherwise
                disp('could not load variance map. No sCMOS noise model used.')
                p.isscmos=false;
                fitpar.issCMOS=false;
                varstack=0;
                varmaph=[];
        end
        if ~isempty(varmaph)
            roi=p.loc_cameraSettings.roi;
            varmap=varmaph(max(1,roi(1)):roi(3),max(1,roi(2)):roi(4));
        end
    else 
        varmap=[];
    end
    fitpar.varmap=varmap*p.loc_cameraSettings.pix2phot;
    if p.userefractive_index_mismatch
        fitpar.refractive_index_mismatch=p.refractive_index_mismatch;
    else
        fitpar.refractive_index_mismatch=1;
    end


% elseif fitpar.fitmode==5
%     calfile=p.cal_3Dfile;
%     cal=load(calfile);
%     fitpar.splinecoefficients=single(cal.cspline.coeff);
%     fitpar.z0=cal.z0;
%     fitpar.dz=cal.dz; 
%     fitpar.refractive_index_mismatch=p.refractive_index_mismatch;
%     fitpar.objPos=p.objPos;
    
else
    fitpar.mode='Gauss';
    fitpar.PSFx0=p.PSFx0;
end
if p.loc_cameraSettings.EMon
    fitpar.EMexcessNoise=2;
else
fitpar.EMexcessNoise=1;
end

end

function varstack=getvarmap(varmap,stackinfo,roisize)
numim=length(stackinfo.xpix);
varstack=zeros(roisize,roisize,numim,'single');
dn=floor(roisize/2);
for k=1:numim
%     stackinfo.x(k)
%     stackinfo.y(k)
    varstack(:,:,k)=varmap(stackinfo.xpix(k)-dn:stackinfo.xpix(k)+dn,stackinfo.ypix(k)-dn:stackinfo.ypix(k)+dn);
end
end

function fitmode_callback(a,b,obj)
p=obj.getGuiParameters;
fitmode=p.fitmode.Value;
% fitz={'loadcal','cal_3Dfile','trefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize','pixelsizex','pixelsizey','automirror','fit2D','isscmos','selectscmos','scmosfile'};
% fitxy={'PSFx0','tPSFx0'};
% switch fitmode
%     case {3,5}
%         ton=fitz;
%         toff=fitxy;
%     otherwise
%         toff=fitz;
%         ton=fitxy;
% end

switch fitmode
    case {1,2}
        roisize=7;
        iterations=30;
        RowName={'x','y','N','Bg','PSF',' '};
      
    otherwise
        roisize=13;
        iterations=50;
        RowName={'x','y','z','N','Bg',' '};
end

obj.setPar('loc_ROIsize',roisize);

% obj.fieldvisibility('on',ton,'off',toff);
obj.setGuiParameters(struct('iterations',iterations));
obj.guihandles.globaltable.RowName=RowName;
end

function pard=guidef(obj)
p1(1).value=1; p1(1).on={'PSFx0','tPSFx0'}; 
p1(1).off={'loadcal','cal_3Dfile','userefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize',...
    'fit2D','isscmos','pixelsizex','pixelsizey','selectscmos','scmosfile'};
p1(2)=p1(1);p1(2).value=2;
p1(3).value=3;p1(3).off={'PSFx0','tPSFx0'};p1(3).on={'loadcal','cal_3Dfile','userefractive_index_mismatch','refractive_index_mismatch','overwritePixelsize','fit2D','isscmos'};
p1(4)=p1(1);p1(4).value=4;
p1(5)=p1(3);p1(5).value=5;
p1(6)=p1(5);p1(6).value=6;

pard.fitmode.object=struct('Style','popupmenu','String',{{'PSF fix','PSF free','3D z','ellipt: PSFx PSFy','Spline'}},'Value',2,'Callback',{{@obj.switchvisible,p1,{@fitmode_callback,0,0,obj}}});
pard.fitmode.position=[1,1];
pard.fitmode.Width=1.5;
pard.fitmode.TooltipString=sprintf('Fit mode. Fit with constant PSF, free PSF, 3D with astigmatism, asymmetric PSF (for calibrating astigmatic 3D)');

pard.text.object=struct('Style','text','String','Iterations:');
pard.text.position=[1,2.5];
pard.text.Width=0.7;
pard.text.Optional=true;
pard.iterations.object=struct('Style','edit','String','50');
pard.iterations.position=[1,3.2];
pard.iterations.TooltipString=sprintf('number of iterations for the GPU fitter (typical 50, use 100-250 for ellipt: PSFx PSFy or 3Dz).');
pard.iterations.Optional=true;
pard.iterations.Width=0.5;

pard.roisperfitt.object=struct('Style','text','String','ROIs/fit:');
pard.roisperfitt.position=[1,3.9];
pard.roisperfitt.Width=0.6;
pard.roisperfitt.Optional=true;
pard.roisperfit.object=struct('Style','edit','String','15000');
pard.roisperfit.position=[1,4.5];
pard.roisperfit.TooltipString=sprintf('Number of 10 x 10 pixel ROIs passed to GPU for fitting. For other ROI sizes, the number is adjusted accordingly.');
pard.roisperfit.Optional=true;
pard.roisperfitt.TooltipString=pard.roisperfit.TooltipString;
pard.roisperfit.Width=0.5;

pard.tPSFx0.object=struct('Style','text','String','PSFx start (pix)');
pard.tPSFx0.position=[2,1];
pard.tPSFx0.Width=1.25;
pard.tPSFx0.Optional=true;

pard.PSFx0.object=struct('Style','edit','String','1');
pard.PSFx0.position=[2,2.25];
pard.PSFx0.Width=0.75;
pard.PSFx0.TooltipString=sprintf('start value for PSF, or size of PSF when PSF fixed (in camera pixels)');
pard.PSFx0.Optional=true;

% pard.fit2D.object=struct('Style','checkbox','String','2D PSF','Value',0);
% pard.fit2D.position=[2,3.5];
% pard.fit2D.Width=.75;
% pard.fit2D.TooltipString=sprintf('Check if PSF model is 2D (no specific PSF engineering), or displays a high degree of similarity above and below the focal plane');
% pard.fit2D.Optional=true;
pard.zstartt.object=struct('Style','text','String','z start (nm)');
pard.zstartt.position=[2,3.5];
pard.zstartt.Width=.75;
pard.zstartt.TooltipString=sprintf('z start parameter. Use vector with several values for 2D PSF');
pard.zstartt.Optional=true;

pard.zstart.object=struct('Style','edit','String','0');
pard.zstart.position=[2,4.25];
pard.zstart.Width=.75;
pard.zstart.TooltipString=pard.zstartt.TooltipString;
pard.zstart.Optional=true;

% pard.automirror.object=struct('Style','popupmenu','String',{{'auto','no mirror','mirror'}});
% pard.automirror.position=[2,4.25];
% pard.automirror.Width=.75;
% pard.automirror.Optional=true;

pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal');
pard.loadcal.position=[2,1];
pard.loadcal.Width=.75;
pard.cal_3Dfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.cal_3Dfile.position=[2,1.75];
pard.cal_3Dfile.Width=1.75;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


p(1).value=0; p(1).on={}; p(1).off={'globaltable','linkt','link','mainchannelt','mainchannel','weightch2t','weightch2'};
p(2).value=1; p(2).on={'globaltable','linkt','link','mainchannelt','mainchannel','weightch2t','weightch2'}; p(2).off={};

pard.isglobal.object=struct('Style','checkbox','String','Global fit','Callback',{{@obj.switchvisible,p}});
pard.isglobal.position=[3,3.];
pard.isglobal.Width=.75;
pard.isglobal.Optional=true;


pard.mainchannelt.object=struct('Style','text','String','main x,y:');
pard.mainchannelt.position=[4,3.];
pard.mainchannelt.Width=.75;
pard.mainchannelt.Optional=true;

pard.mainchannel.object=struct('Style','popupmenu','String',{{'mean','ch1','ch2'}});
pard.mainchannel.position=[4,3.5];
pard.mainchannel.Width=.75;
pard.mainchannel.Optional=true;

pard.weightscht.object=struct('Style','text','String','Weights r t');
pard.weightscht.position=[5,3.];
pard.weightscht.Width=.75;
pard.weightscht.Optional=true;
pard.weightsch.object=struct('Style','edit','String','1 1');
pard.weightsch.position=[5,3.75];
pard.weightsch.Width=.5;
pard.weightsch.Optional=true;


pard.globaltable.object=struct('Style','listbox','String','x');
pard.globaltable.position=[7,4.25];
pard.globaltable.Width=.75;
pard.globaltable.Optional=true;
pard.globaltable.Height=5;
% pard.linkt.object=struct('Style','text','String','link: x y z N bg');
% pard.linkt.position=[3,2];
% pard.linkt.Width=1.5;
% pard.linkt.Optional=true;

% pard.link.object=struct('Style','edit','String','1 1 1 0 0');
% pard.link.position=[3,3.5];
% pard.link.Width=1.5;
% pard.link.Optional=true;
p(1).value=0; p(1).on={}; p(1).off={'refractive_index_mismatch'};
p(2).value=1; p(2).on={'refractive_index_mismatch'}; p(2).off={};
pard.userefractive_index_mismatch.object=struct('Style','checkbox','String','RI mismatch:','Callback',{{@obj.switchvisible,p}});
pard.userefractive_index_mismatch.position=[3,1];
pard.userefractive_index_mismatch.Width=1.25;
pard.userefractive_index_mismatch.Optional=true;


pard.refractive_index_mismatch.object=struct('Style','edit','String','.8');
pard.refractive_index_mismatch.position=[3,2.25];
pard.refractive_index_mismatch.TooltipString=sprintf('Correction factor to take into account the different refracrive indices of immersion oil and buffer. \n This leads to smaller distances inside the sample compared to bead calibration. \n Bead calibration: in piezo positions (nm). \n This factor transforms z positions to real-space z positions. \n For high-NA oil objectives: typical 0.72 (range 0.7-1).');
pard.refractive_index_mismatch.Optional=true;
pard.refractive_index_mismatch.Width=0.35;


p(1).value=0; p(1).on={}; p(1).off={'pixelsizex','pixelsizey'};
p(2).value=1; p(2).on={'pixelsizex','pixelsizey'}; p(2).off={};
pard.overwritePixelsize.object=struct('Style','checkbox','String','pixelsize X,Y (um):','Callback',{{@obj.switchvisible,p}});
pard.overwritePixelsize.position=[4,1];
pard.overwritePixelsize.Width=1.25;
pard.overwritePixelsize.Optional=true;

pard.pixelsizex.object=struct('Style','edit','String','.1');
pard.pixelsizex.position=[4,2.25];
pard.pixelsizex.Width=0.35;
pard.pixelsizex.Optional=true;

pard.pixelsizey.object=struct('Style','edit','String','.1');
pard.pixelsizey.position=[4,2.6];
pard.pixelsizey.Width=0.35;
pard.pixelsizey.Optional=true;

p(1).value=0; p(1).on={}; p(1).off={'selectscmos','scmosfile'};
p(2).value=1; p(2).on={'selectscmos','scmosfile'}; p(2).off={};
pard.isscmos.object=struct('Style','checkbox','String','sCMOS','Callback',{{@obj.switchvisible,p}});   
pard.isscmos.position=[7,1];
pard.isscmos.Optional=true;
pard.selectscmos.object=struct('Style','pushbutton','String','Load var map','Callback',{{@loadscmos_callback,obj}});   
pard.selectscmos.TooltipString='Select sCMOS variance map (in ADU^2) of same size ROI on chip as image stack';
pard.selectscmos.position=[7,1.6];
pard.selectscmos.Optional=true;
pard.selectscmos.Width=.9;
pard.scmosfile.object=struct('Style','edit','String','');
pard.scmosfile.TooltipString='Tiff/.mat image containing sCMOS variance map (same ROI on camera as tiff).';
pard.scmosfile.position=[7,2.5];
pard.scmosfile.Optional=true;
    pard.scmosfile.Width=.5;
    
pard.asymmetry.object=struct('Style','checkbox','String','get asymmetry');   
pard.asymmetry.position=[5,1];
pard.asymmetry.Optional=true;
    


pard.syncParameters={{'cal_3Dfile','',{'String'}}};

pard.plugininfo.type='WorkflowFitter';
pard.plugininfo.description='Global fitter for multiple channels.';
end