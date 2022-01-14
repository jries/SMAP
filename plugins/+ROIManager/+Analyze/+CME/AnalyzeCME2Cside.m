classdef AnalyzeCME2Cside<interfaces.DialogProcessor&interfaces.SEProcessor
    methods
        function obj=AnalyzeCME2Cside(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            if nargin>0
                obj.handle=varargin{1};
            end
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
          out=analyze2Cside(obj.SE,p);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.export_selected.object=struct('String','some parameters','Style','checkbox');
pard.export_selected.position=[1,1];
pard.export_selected.Width=2;

pard.plugininfo.type='ROI_Analyze';
end

function out=analyze2Cside(SE,p)
global ccamplitude dx dy
out=[]
%XXXX introduce symmetry by using mirrored images
%in between: at some point perfrom rotational alignment with class
%averages, then mirror again.

% calculate feature vector from evaluation results
featureVectorEv=buildFeatureVector(SE.sites);
% k-means to determine classes
timewindows=5;
[idx,C]=kmeans(featureVectorEv,timewindows,'Replicates',100);
idxnew=sortFV(featureVectorEv,timewindows,idx);
images=getFieldAsVector(SE.sites,'evaluation','CME2CSide','imageshift');

% for CC: reweigh channels
weights=[1,1];
imagesw=reweighImages(images,weights);

% calculate class averages
hfig=figure(79);
classaverages=getShiftClassAverages(imagesw,idxnew, [],[],hfig);


% CC each image with each class -> similarity vectors
[ccamplitude, dx, dy]=correlateImages(imagesw,classaverages);

[n,indmax]=max(ccamplitude,[],2);
ds=randn(length(idxnew),1)*.1;
ds2=randn(length(idxnew),1)*.1;
figure(88);plot(idxnew+ds,indmax+ds2,'*')
%apply dx, dy, recalculate classaverages
hfig=figure(80);
[classaverages,shiftimages]=getShiftClassAverages(imagesw,idxnew, dx,dy, hfig);
%correlate again -> for FV
%faster: just take shifted image and calculate product with class (no extra
%determination of shift needed?) or onnly next iteration
% add similarities to feature vectors
ccscale=.001;
sfv=size(featureVectorEv);
featureVectorCC=horzcat(featureVectorEv,(ccamplitude-min(ccamplitude(:)))/ccscale);
%k -means -> new classes
[idxcc,C]=kmeans(featureVectorCC,timewindows,'Replicates',100);
idxnewcc=sortFV(featureVectorCC(:,1:sfv(2)),timewindows,idxcc);
% calculate shifted images from CC

% calculate new class averages
hfig=figure(83);
classaverages=getShiftClassAverages(shiftimages,idxnewcc,[],[], hfig);




return

% global abga dxa dya
% im=getImageAligned(SE.sites(1),p);
ch=2;

if isempty(abga)
[dxa,dya,abga]=getCorrelationMatrix(SE,ch);
end

featureV=buildFeatureVector(abga,SE.sites);
timewindows=15;
[idx,C]=kmeans(featureV,timewindows,'Replicates',100);
sfv=size(featureV);

% y=mdscale(abga,1);
% [ys,ind]=sort(y);

%try sorting by PCA
[coeff,score]=pca(featureV);
[~,ind]=sort(score(:,1));

figure(99)
for k=1:length(ind)
    subplot(6,5,k)
    imh=SE.sites(ind(k)).evaluation.CME2CSide.imageshift;
    imh(:,:,2)=imh(:,:,2)/5;
    imagesc(imh);
    
end

out.a=abga;out.dx=dxa;out.dy=dya;
end

function imout=shiftimage(imin,dx,dy)
if dx==0&&dy==0
    imout=imin;
    return
end
sim=size(imin);
imout=zeros(sim);
dxr=round(dx);dyr=round(dy);

rx=1:sim(1)-abs(dxr);
rx2=abs(dxr)+1:sim(1);
if dxr<0
    rxa=rx;rxb=rx2;
else
    rxa=rx2;rxb=rx;
end

ry=1:sim(2)-abs(dyr);
ry2=abs(dyr)+1:sim(2);
if dyr<0
    rya=ry;ryb=ry2;
else
    rya=ry2;ryb=ry;
end

imout(rxa,rya,:)=imin(rxb,ryb,:);
end

function idxnew=sortFV(fv,timewindows,idx)
sfv=size(fv);
meanfv=zeros(timewindows,sfv(2));
for k=1:timewindows
    incluster=(idx==k);
    meanfv(k,:)=mean(fv(incluster,:),1);
end
sortfv=mean(meanfv,2);
[~,indsort]=sort(sortfv);
[~,indsort2]=sort(indsort);
 idxnew=indsort2(idx);
end


% function cav=getClassAverages(images,idx,hfig)
% timewindows=max(idx);
% imlines={};
% mh=max(hist(idx));
% for k=1:timewindows
%     incluster=find(idx==k);
%     imline=[];
%     avim=0;
%     for l=1:length(incluster)
%      imh=images{incluster(l)};
%      imh(:,:,2)=imh(:,:,2)/3;
%      imline=horzcat(imline,imh);
%      avim=imh+avim;
%     end
%     cav{k}=getSymmetricImage(avim/length(incluster));
%     imblank=0*imh;imblank(:,:,3)=.2;
%     imline=horzcat(imline,repmat(imblank,1,mh-length(incluster)));
%     imlines{k}=horzcat(imline,cav{k});
% end
% 
% if ~isempty(hfig)
% imall=[];
% for k=1:timewindows
%      imall=vertcat(imall,imlines{k});
% end
% figure(hfig)
% imagesc(imall)
% axis equal
% axis tight
% end
% end

function [cav,shiftimages]=getShiftClassAverages(images,idx,dx,dy,hfig)
if isempty(dx)
    dx=zeros(length(idx),max(idx));
end
if isempty(dy)
    dy=zeros(length(idx),max(idx));
end
%later only one get Class Averages
timewindows=max(idx);
imlines={};
mh=max(hist(idx));
shiftimages{length(images)}=0*images{1};
for k=1:timewindows
    incluster=find(idx==k);
    imline=[];
    avim=0;
    for l=1:length(incluster)
     imh=images{incluster(l)};
     shiftimages{incluster(l)}=imh;
     imh(:,:,2)=imh(:,:,2)/3;
     imline=horzcat(imline,imh);
     dxh=dx(incluster(l),k);
     dyh=dy(incluster(l),k);
     imh=shiftimage(imh,dxh,dyh);
     
     avim=imh+avim;
    end
    cav{k}=getSymmetricImage(avim/length(incluster));
%     cav{k}=avim/length(incluster);
    imblank=0*imh;imblank(:,:,3)=.2;
    imline=horzcat(imline,repmat(imblank,1,mh-length(incluster)));
    imlines{k}=horzcat(imline,cav{k});
end

if ~isempty(hfig)
imall=[];
for k=1:timewindows
     imall=vertcat(imall,imlines{k});
end
figure(hfig)
imagesc(imall)
axis equal
axis tight
end
end

function imo=getSymmetricImage (in)
% in=shiftimage(in,0,30);
inm=in(:,end:-1:1,:);
[dx,dy,abg]=getShiftCorr(sum(in,3),sum(inm,3),0,100);
dy1=round(dy/2);dy2=-round(dy-dy1);
imo=(shiftimage(inm,0,dy1)+shiftimage(in,0,dy2))/2;
end

function fv=buildFeatureVector(sites)
lsN=10000;lsLengthz=30;lsLengthx=200; %length scales
% facCC=1;
Nmds=0;


fields=struct('N1',lsN,'N2',lsN,'zmedd',lsLengthz,'zmd',lsLengthz,'qdx1',lsLengthx,'qdx2',lsLengthx,'qdz1',lsLengthz,'qdz2',lsLengthz);
fn=fieldnames(fields);

lsites=length(sites);
fv=zeros(lsites,Nmds+length(fn));

 for k=1:lsites
     sites(k).evaluation.CME2CSide.qdx1=sites(k).evaluation.CME2CSide.qx1(2)-sites(k).evaluation.CME2CSide.qx1(1);
     sites(k).evaluation.CME2CSide.qdx2=sites(k).evaluation.CME2CSide.qx2(2)-sites(k).evaluation.CME2CSide.qx2(1);
     sites(k).evaluation.CME2CSide.qdz1=sites(k).evaluation.CME2CSide.qz1(2)-sites(k).evaluation.CME2CSide.qz1(1);
     sites(k).evaluation.CME2CSide.qdz2=sites(k).evaluation.CME2CSide.qz2(2)-sites(k).evaluation.CME2CSide.qz2(1);
 end
for k=1:lsites
    evc=sites(k).evaluation.CME2CSide;
    for l=1:length(fn)
        fv(k,l)=evc.(fn{l})/fields.(fn{l});
    end   
end
end

function imout=reweighImages(images,weights)
if length(weights)==2
    weights(3)=0;
end
for k=length(images):-1:1
    imtemp=images{k};
    imtemp(:,:,1)=imtemp(:,:,1)*weights(1)/sum(weights);
    imtemp(:,:,2)=imtemp(:,:,2)*weights(2)/sum(weights);
    imtemp(:,:,3)=imtemp(:,:,3)*weights(3)/sum(weights);
    imout{k}=(imtemp);
end
end

function [ccamplitude, dxa, dya]=correlateImages(images,classaverages)
li=length(images);
lc=length(classaverages);
dxa=ones(li,lc);dya=ones(li,lc);ccamplitude=ones(li,lc);
for k2=1:lc
    for k1=1:li    
         figure(77)
        [dx,dy,abg]=getShiftCorr(sum(images{k1},3),sum(classaverages{k2},3),false,100,false);
        title(k2)
        drawnow
        ccamplitude(k1,k2)=abg;
        dxa(k1,k2)=dx;
        dya(k1,k2)=dy;
    end
end
end


% function [dxa,dya,abga]=getCorrelationMatrix(SE,ch)
% l=length(SE.sites);
% dxa=ones(l);dya=ones(l);abga=ones(l);
% 
% for k1=1:length(SE.sites);
% %     k1
%     for k2=k1+1:length(SE.sites);
% %         k2
%         im1=SE.sites(k1).evaluation.CME2CSide.imageshift;
%         im2=SE.sites(k2).evaluation.CME2CSide.imageshift;
%         figure(81)
%         imagesc([im1(:,:,ch),im2(:,:,ch)])
%         
%         figure(80)
%         [dx,dy,abg]=getShiftCorr(1*im1(:,:,ch),1*im2(:,:,ch),1,100);
%         abga(k1,k2)=abg;
%         dxa(k1,k2)=dx;
%         dya(k1,k2)=dy;
%         
%         abga(k2,k1)=abg;
%         dxa(k2,k1)=dx;
%         dya(k2,k1)=dy;
%         drawnow
%     end
% end
% end
% 
% function im=getImageAligned(site,p)
% 
% end