classdef FRCresolution<interfaces.DialogProcessor
    %FRCresolution calculates FRC resolution
    %     FRC resolution implementation based on the matlab code provided with: 
    %[1]	R. P. J. Nieuwenhuizen, K. A. Lidke, M. Bates, D. L. Puig, D. Gr?nwald,
    %S. Stallinga, and B. Rieger, ?Measuring image resolution in optical nanoscopy.,? 
    %Nat Methods, vol. 10, no. 6, pp. 557?562, Jun. 2013.
    methods
        function obj=FRCresolution(varargin)           
            obj@interfaces.DialogProcessor(varargin{:}) 
            obj.inputParameters=anyRender;
            obj.showresults=true;
            obj.history=false;
        end
        
        function out=run(obj,p)
            %binning: works only with layer1.
            out=[];
             if p.takeimage
                 p.pixrec_frc=p.sr_pixrec;
             end
            lochere=obj.locData.copy;
            %function for following
%             bfrc=round(rand(length(lochere.loc.xnm),1));
            bfrc=makeblocks(p,length(lochere.loc.xnm));
%             bfrc=zeros(length(lochere.loc.xnm),1);
%             bfrc(1:end/2)=1;
            lochere.addloc('FRCblocks',bfrc);
            lochere.regroup;
            lochere.filter('FRCblocks',[],'minmax',[-.1 .5])
            image1=getimage(lochere,p);
            lochere.filter('FRCblocks',[],'minmax',[.5 1.1])
            image2=getimage(lochere,p);
            
            image1f=filterimage(image1);
            image2f=filterimage(image2);
            frc_curve=getFRC(image1f,image2f);
            [FRC_resolutionp,rl,rh]=frctoresolution(frc_curve,size(image1));
            FRC_resolution=FRC_resolutionp*p.pixrec_frc;
            errres=(rh-rl)/2*p.pixrec_frc;
            ax1=obj.initaxis('frc');
             qmax = 0.5/(p.pixrec_frc);
            
            plot([0 qmax],[0 0],'k')
            hold on
            plot(linspace(0,qmax, length(frc_curve)), frc_curve,'-')
%             plot(linspace(0,qmax*sqrt(2), length(frc_curve)), frc_curve,'-')
            plot([0 qmax],[1/7 1/7],'m')
            plot(1/(FRC_resolution),1/7,'rx')
            plot(1/(FRC_resolution)*[1 1],[-0.2 1/7],'r')
            hold off
            xlim([0,qmax]);
            ylim([-0.2 1.2])
            xlabel('Spatial frequency (nm^{-1})');
            ylabel('FRC')
            title(['FRC reoslution (nm): ' num2str(FRC_resolution,'%4.1f') ' +/- ' num2str(errres,'%4.1f')])
            
            out.clipboard=['FRC reoslution (nm): ' 9 num2str(FRC_resolution,'%4.1f') 9 ' +/- ' 9 num2str(errres,'%4.1f')];
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function  bfrc=makeblocks(p,l)
bfrc=false(l,1);
blocksize=floor(l/p.frc_blocks);
for k=1:p.frc_blocks
    indrange=(k-1)*blocksize+1:k*blocksize;
    if strcmp(p.blockassignment.selection,'random')
        side=rand>0.5;
    else
        side=mod(k,2);
    end
    bfrc(indrange)=side;
end
end

function frc_out=getFRC(image1,image2)
in1=fftshift(fft2(image1));
in2=fftshift(fft2(image2));
inc=in1.*conj(in2);
frc_num=real(radialsum(inc));
in1 = abs(in1).^2;
in2 = abs(in2).^2;
frc_denom = sqrt(abs(radialsum(in1).*radialsum(in2)));                      % Denominator
frc_out = double(frc_num)./double(frc_denom);                               % FRC
frc_out(isnan(frc_out)) = 0;  
end

function rs=radialsum(img)
s=size(img);
center=floor((s+1)/2);
rs=zeros(ceil(s(1)/2)+1,1);
for k=1:s(1)
    for l=1:s(2)
        d=sqrt((k-center(1)).^2+(l-center(2)).^2);
        ind=round(d)+1;
        if ind<=length(rs)
        rs(ind)=rs(ind)+img(k,l);
        end
    end
end
end
function imageo=filterimage(imagei)
sz=size(imagei);
nfac = 8;    
% Image width / Width of edge region
x=((1:sz(1))-sz(1)/2)/sz(1);
x_im = meshgrid(x);
mask = 0.5-0.5*cos(pi*nfac*x_im);          
mask(abs(x_im)<((nfac-2)/(nfac*2))) = 1;
maskt=mask.*mask';
imageo=maskt.*imagei;

end
function image=getimage(lochere,p)
if p.takeimage
    p.sr_axes=[];
    imt=anyRender(lochere,p);
    image=sum(imt.composite,3);
else
    
    locs=lochere.getloc({'xnm','ynm','FRCblocks'},'position','roi','layer',1);
    
    rx=p.sr_pos(1)+p.sr_size(1)*[-1 1];
    ry=p.sr_pos(2)+p.sr_size(2)*[-1 1];
    image=myhist2(locs.xnm,locs.ynm,p.pixrec_frc,p.pixrec_frc,rx,ry);
    
end
s=size(image);
if s(2)~=s(1)
    image(max(s(1:2)),max(s(1:2)))=0;
end
end

function varargout=frctoresolution(frc_in,sz)
sz=sz(1);
%% Check for undesirable values in FRC curve
frc_in = real(frc_in);
% frc_in(isnan(frc_in)) = 0;
frc_in(frc_in>1) = 1;
frc_in(frc_in<-1) = -1;

%% Calculate the number of pixels in a Fourier ring
% Number of pixels in a Fourier ring
nr = double(radialsum(ones(sz)));       

if length(nr) ~= length(frc_in)
   error('The specified image size could not have resulted in this FRC curve.'); 
end

% Calculate the effective width of pixels in a Fourier ring
if nargout >1
    n1 = find(frc_in<0, 1,'first' );
    frc_tmp = frc_in;
    frc_tmp(isnan(frc_tmp)) = 0;
    frc_tmp(isinf(frc_tmp)) = 0;
    ne_inv = mean(frc_tmp(n1:end).^2.*nr(n1:end));
end

%% Smoot the FRC curve 
% Least squares interpolation for curve smoothing
sspan = ceil(sz/20);      % Smoothing span
if (sz/20)<5
    sspan = 5;
end
sspan = sspan  + (1-mod(sspan,2));

% if TB_curve
%    p = pwd; % hack to avoid the function shadwoing by smooth from dip_image
%    cd([TB_d filesep 'curvefit'])
   frc_in = double(smooth(frc_in,sspan,'loess'));        
%    cd(p)
% else
%    frc_in = [double(gaussf(frc_in,.9))]';
% end
q = (0:(length(frc_in)-1))'/sz;                   % Spatial frequencies

%% Calculate intersections between the FRC curve and the threshold curve
% isects = polyxpoly(q,frc_in,q,thresholdcurve);
thresholdcurve = 1/7*ones(size(frc_in));
isects = isect(q,frc_in,thresholdcurve);

%% Find first intersection to obtain the resolution

% Throw away intersections at frequencies beyond the Nyquist frequency
isects = isects(isects<0.5);

if ~isempty(isects)
    % Find the first intersection where the FRC curve is decreasing
    isect_inds = 1+floor(sz*isects);    % Indices of the intersections
    for ii = 1:length(isect_inds)
        isect_ind = isect_inds(ii);
        if frc_in(isect_ind+1) < frc_in(isect_ind)
            resolution = 1/isects(ii);
            break
        end
    end
end

if ~exist('resolution','var')
    resolution = -1;
end

if resolution == -1
    varargout{1} = NaN;
    varargout{2} = 0;
    varargout{3} = 0.5;
    fprintf(' -- Could not find the resolution --\n')
    return
end

if nargout < 2
    varargout{1} = resolution;
    return
end

%% Calculate the resolution uncertainty

% Variance of the unsmoothed FRC values
frc_var = (1+2*frc_in-frc_in.^2).*(1-frc_in).^2*ne_inv./nr;

% Indices of FRC values inside the smoothing window
indexvec = (ceil(isect_ind-sspan/2)):(floor(isect_ind+sspan/2));    
indexvec = indexvec(indexvec>0);
indexvec = indexvec(indexvec<=length(frc_in));

% Calculate the covariance matrix of the LMS fit parameters.
S = diag(frc_var(indexvec));    % Approximate covariance matrix of the FRC values
X = cat(2,ones(size(indexvec')),q(indexvec),q(indexvec).^2);
C = inv(X'*X);
C = C*X'*S*X*C;

% Estimate the variance of FRC(fres) resulting from the least mean
% squares fit. The factor ne_inv arises to correct for the
% correlations between neighbouring FRC values.
frc_newvar = ne_inv*(1./[1 resolution resolution^2])*C*(1./[1; resolution; resolution^2]);

if isnan(frc_newvar) || isinf(frc_newvar)
    frc_newvar = frc_var(isect_ind);
end

% Calculate the crossing between FRC + 1*stdev(FRC) and threshold
frc_high = frc_in+sqrt(frc_newvar);
isects_high = isect(q,frc_high,thresholdcurve);
isects_high = isects_high(isects_high<0.5);
% Calculate the crossing between FRC - 1*stdev(FRC) and threshold 
frc_low = frc_in-sqrt(frc_newvar);
isects_low = isect(q,frc_low,thresholdcurve);
isects_low = isects_low(isects_low<0.5);

% Find resolution + 1 standard deviation
if ~isempty(isects_high)
    isect_inds = 1+floor(sz*isects_high);
    for ii = 1:length(isect_inds)
        isect_ind = isect_inds(ii);
        if frc_high(isect_ind+1) < frc_high(isect_ind)
            resolution_high = 1/isects_high(ii);
            break
        end
    end
end
if ~exist('resolution_high','var')
    resolution_high = 0;
end

% Find resolution + - standard deviation
if ~isempty(isects_low)
    isect_inds = 1+floor(sz*isects_low);
    for ii = 1:length(isect_inds)
        isect_ind = isect_inds(ii);
        if frc_low(isect_ind+1) < frc_low(isect_ind)
            resolution_low = 1/isects_low(ii);
            break
        end
    end
end
if ~exist('resolution_low','var')
    resolution_low = 0.5;
end

%% Write outputs
varargout{1} = resolution;
varargout{2} = resolution_high;
varargout{3} = resolution_low;


end


function pard=guidef
pard.t0.object=struct('String','FRC resolution. Uses ROI or FoV (if rendered image)','Style','text');
pard.t0.position=[1,1];
pard.t0.Width=4;

pard.takeimage.object=struct('String','use settings from rendered image','Style','checkbox','Value',0);
pard.takeimage.position=[2,1];
pard.takeimage.Width=2;
pard.t1.object=struct('String','otherwise:','Style','text');
pard.t1.position=[3,1];
pard.t2.object=struct('String','pixelsize (nm)','Style','text');
pard.t2.position=[4,1];

pard.pixrec_frc.object=struct('String','5','Style','edit');
pard.pixrec_frc.position=[4,2];

pard.t3.object=struct('String','binning of the data in blocks','Style','text');
pard.t3.position=[2,3];
pard.t3.Width=2;
pard.t4.object=struct('String','number of blocks','Style','text');
pard.t4.position=[3,3];

pard.frc_blocks.object=struct('String','10','Style','edit');
pard.frc_blocks.position=[3,4];

pard.blockassignment.object=struct('String',{{'alternating','random'}},'Style','popupmenu');
pard.blockassignment.position=[4,4];

pard.plugininfo.name='FRC resolution';
pard.plugininfo.description='FRC resolution. Uses ROI or FoV (if rendered image is ticked). implementation based on the matlab code provided with: [1]	R. P. J. Nieuwenhuizen, K. A. Lidke, M. Bates, D. L. Puig, D. Gr?nwald, S. Stallinga, and B. Rieger, ?Measuring image resolution in optical nanoscopy.,? Nat Methods, vol. 10, no. 6, pp. 557?562, Jun. 2013.';
pard.plugininfo.type='ProcessorPlugin';
end