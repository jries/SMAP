%Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7

function simplefitter_cspline(p)
%parameters:
% p.imagefile: fiename of data (char);
% p.calfile: filename of calibration data (char);
% p.offset=ADU offset of data;
% p.conversion=conversion e-/ADU;
% p.preview: true if preview mode (fit only current image and display
% results).
% p.previewframe=frame to preview;
% p.peakfilter=filtersize (sigma, Gaussian filter) for peak finding;
% p.peakcutoff=cutoff for peak finding
% p.roifit=size of the ROI in pixels
% p.bidirectional= use bi-directional fitting for 2D data
% p.mirror=mirror images if bead calibration was taken without EM gain
% p.status=handle to a GUI object to display the status;
% p.outputfile=file to write the localization table to;
% p.outputformat=Format of file;
% p.pixelsize=pixel size in nm;

% p.loader which loader to use
% p.mij if loader is fiji: this is the fiji handle
% p.isscmos scmos camera used
% p.scmosfile file containgn scmos varmap


global simplefitter_stop

fittime=0;
fitsperblock=50000;
imstack=zeros(p.roifit,p.roifit,fitsperblock,'single');
peakcoordinates=zeros(fitsperblock,3);
indstack=0;
resultsind=1;
% bgmode='wavelet'; 
if contains(p.backgroundmode,'avelet')
    bgmode=2;
elseif contains(p.backgroundmode,'aussian')
     bgmode=1;
else
    bgmode=0;
end

%scmos
varmap=[];
if p.isscmos
    varstack=ones(p.roifit,p.roifit,fitsperblock,'single');
    [~,~,ext]=fileparts(p.scmosfile);
    switch ext
        case '.tif'
            varmap=imread(p.scmosfile);
        case '.mat'
            varmap=load(p.scmosfile);
            if isstruct(varmap)
                fn=fieldnames(varmap);
                varmap=varmap.(fn{1});
            end
        otherwise
            errordlg('could not load variance map. No sCMOS noise model used.')
            p.isscmos=false;
            varstack=0;
    end
else
    varstack=0;
end
varmap=varmap*p.conversion^2;
%results
% frame, x,y,z,phot,bg, errx,erry, errz,errphot, errbg,logLikelihood
%load calibration
if exist(p.calfile,'file')
    cal=load(p.calfile);
    p.dz=cal.cspline.dz;  %coordinate system of spline PSF is corner based and in units pixels / planes
    p.z0=cal.cspline.z0;
    p.coeff=cal.cspline.coeff;
    if iscell(p.coeff)
        p.coeff=p.coeff{1};
    end
    p.isspline=true;
else
%     errordlg('please select 3D calibration file')
    warndlg('3D calibration file could not be loaded. Using Gaussian fitter instead.','Using Gaussian fit','replace');
    p.isspline=false;
end

p.dx=floor(p.roifit/2);
% readerome=bfGetReader(p.imagefile);
p.status.String=['Open tiff file' ]; drawnow

 switch p.loader
     case 1
        reader=mytiffreader(p.imagefile);
        numframes=reader.info.numberOfFrames;
     case 2
        reader=bfGetReader(p.imagefile);
        numframes=reader.getImageCount;
     case 3 %fiji
          ij=p.mij.imagej;
          ijframes=ij.getFrames;
          for k=1:length(ijframes)
            if strcmp(ijframes(k).class,'ij.gui.StackWindow')&&~isempty(ijframes(k).getImagePlus)    
                reader=ijframes(k).getImagePlus.getStack;
                break
            end
          end
          if ~exist('reader','var')
              p.status.String='Error... Check if image is loaded in ImageJ'; drawnow
              return
          end
%           numframes=reader.size;
          numframes=reader.getSize;
 end


if p.preview
    frames=min(p.previewframe,numframes);
else
    frames=1:numframes;
end

%loop over frames, do filtering/peakfinding
hgauss=fspecial('gaussian',max(3,ceil(3*p.peakfilter+1)),p.peakfilter);
rsize=max(ceil(6*p.peakfilter+1),3);
hdog=fspecial('gaussian',rsize,p.peakfilter)-fspecial('gaussian',rsize,max(1,2.5*p.peakfilter));
tshow=tic;
for F=frames
    image=getimage(F,reader,p);

    sim=size(image);
    imphot=(single(image)-p.offset)*p.conversion;
    
    %background determination
    if 0 %bgmode==3% wavelet
%         bg=mywaveletfilter(imphot,3,false,true);
%         
%         impf=filter2(hgauss,(imphot)-(bg));
    elseif bgmode==1
%         impf=filter2(hdog,sqrt(imphot));
%          impf=filter2(hdog,(imphot));
         impf=filter2(hdog,(imphot- min(imphot(:,1))));
    elseif bgmode==0
        impf=filter2(hgauss,(imphot));
    end
    maxima=maximumfindcall(impf);
    indmgood=maxima(:,3)>(p.peakcutoff);
    indmgood=indmgood&maxima(:,1)>p.dx &maxima(:,1)<=sim(2)-p.dx;
    indmgood=indmgood&maxima(:,2)>p.dx &maxima(:,2)<=sim(1)-p.dx;
    maxgood=maxima(indmgood,:);

    
    if p.preview && size(maxgood,1)>2000
        p.status.String=('increase cutoff');
        return
    elseif p.preview && size(maxgood,1)==0
                p.status.String=('No localizations found, decrease cutoff');
        return
    end
    
    %cut out images
    for k=1:size(maxgood,1)
        if maxgood(k,1)>p.dx && maxgood(k,2)>p.dx && maxgood(k,1)<= sim(2)-p.dx && maxgood(k,2)<=sim(1)-p.dx 
            indstack=indstack+1;
            if p.mirror
                imstack(:,:,indstack)=imphot(maxgood(k,2)-p.dx:maxgood(k,2)+p.dx,maxgood(k,1)+p.dx:-1:maxgood(k,1)-p.dx);
            else
                imstack(:,:,indstack)=imphot(maxgood(k,2)-p.dx:maxgood(k,2)+p.dx,maxgood(k,1)-p.dx:maxgood(k,1)+p.dx);
            end
            if p.isscmos
                varstack(:,:,indstack)=varmap(maxgood(k,2)-p.dx:maxgood(k,2)+p.dx,maxgood(k,1)-p.dx:maxgood(k,1)+p.dx);
            end
            peakcoordinates(indstack,1:2)=maxgood(k,1:2);
            peakcoordinates(indstack,3)=F;
   
            if indstack==fitsperblock
                p.status.String=['Fitting...' ]; drawnow
                t=tic;
                resultsh=fitspline(imstack,peakcoordinates,p,varstack);
                fittime=fittime+toc(t);
                
                results(resultsind:resultsind+fitsperblock-1,:)=resultsh;
                resultsind=resultsind+fitsperblock;
                
                indstack=0;
            end
        end
       
    end
    if toc(tshow)>1
        tshow=tic;
        p.status.String=['Loading frame ' num2str(F) ' of ' num2str(numframes)]; drawnow
    end
    if  simplefitter_stop
        break
    end
end
closereader(reader,p);
p.status.String=['Fitting last stack...' ]; drawnow
if indstack<1
    p.status.String=['No localizations found. Increase cutoff?' ]; drawnow
else

t=tic;
if p.isscmos
    varh=varstack(:,:,1:indstack);
else
    varh=0;
end
resultsh=fitspline(imstack(:,:,1:indstack),peakcoordinates(1:indstack,:),p,varh); %fit all the rest
fittime=fittime+toc(t);

results(resultsind:resultsind+indstack-1,:)=resultsh;
end
if p.preview
    figure(201)
%     imagesc(impf.^2);
    imagesc(impf);
     colorbar
    hold on
    plot(maxgood(:,1),maxgood(:,2),'wo')
    plot(results(:,2),results(:,3),'k+')
    hold off
   
    p.status.String=['Preview done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations.']; drawnow
else
    p.status.String=['Fitting done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations. Saving now.']; drawnow
    results(:,[13,15])=results(:,[2, 7])*p.pixelsize(1);
    results(:,[14,16])=results(:,[3, 8])*p.pixelsize(end);
    
    if p.isspline
        resultstable=array2table(results,'VariableNames',{'frame','x_pix','y_pix','z_nm','photons','background',' crlb_x','crlb_y','crlb_z','crlb_photons','crlb_background','logLikelyhood','x_nm','y_nm','crlb_xnm','crlb_ynm'});
    else
        resultstable=array2table(results,'VariableNames',{'frame','x_pix','y_pix','sx_pix','sy_pix','photons','background',' crlb_x','crlb_y','crlb_photons','crlb_background','logLikelyhood','x_nm','y_nm','crlb_xnm','crlb_ynm'});
    end
    % 
    writenames=true;
    if contains(p.outputformat,'pointcloud')
        resultstable=resultstable(:,[1 13 14 4 15 9]); %for pointcloud-loader
        del='\t';
        disp('Load in http://www.cake23.de/pointcloud-loader/')
    elseif contains(p.outputformat,'ViSP')
        writenames=false;
        del='\t';
         resultstable=resultstable(:,[13 14 4 15 16 9 5 1]);
         [path,file]=fileparts(p.outputfile);
         p.outputfile=fullfile(path, [file '.3dlp']);
         disp('Load in Visp: https://science.institut-curie.org/research/multiscale-physics-biology-chemistry/umr168-physical-chemistry/team-dahan/softwares/visp-software-2/')
    else
        del=',';
         disp('Generic output. Can be imported e.g. in PALMsiever: https://github.com/PALMsiever/palm-siever')
    end

    writetable(resultstable,p.outputfile,'Delimiter',del,'FileType','text','WriteVariableNames',writenames);
    p.status.String=['Fitting done. ' num2str(size(results,1)/fittime,'%3.0f') ' fits/s. ' num2str(size(results,1),'%3.0f') ' localizations. Saved.']; drawnow
end
end

function results=fitspline(imstack,peakcoordinates,p,varstack)
if p.isspline
    if p.bidirectional
        fitmode=5;
        zstart=[-300 300]/p.dz;
    else
        fitmode=5;
        zstart=0;
    end
    fitpar=single(p.coeff);
else
    if p.bidirectional
        fitmode=2;
    else
        fitmode=4;
    end
    fitpar=single(1);
end

[Pcspline,CRLB,LL]=mleFit_LM(imstack,fitmode,50,fitpar,varstack,1,zstart);
results=zeros(size(imstack,3),12);
results(:,1)=peakcoordinates(:,3);
if  p.mirror
    results(:,2)=p.dx-Pcspline(:,2)+peakcoordinates(:,1);
else
    results(:,2)=Pcspline(:,2)-p.dx+peakcoordinates(:,1);
    
end


if p.isspline
    % frame, x,y,z,phot,bg, errx,erry, errz,errphot, errbg,logLikelihood
results(:,3)=Pcspline(:,1)-p.dx+peakcoordinates(:,2); %x,y in pixels 
results(:,4)=(Pcspline(:,5)-p.z0)*p.dz;
results(:,5:6)=Pcspline(:,3:4);
results(:,7:8)=real(sqrt(CRLB(:,[2 1])));
results(:,9)=real(sqrt(CRLB(:,5)*p.dz));
results(:,10:11)=real(sqrt(CRLB(:,3:4)));
results(:,12)=LL;
else
        % frame, x,y,sx,sy,phot,bg, errx,erry,errphot, errbg,logLikelihood
results(:,3)=Pcspline(:,1)-p.dx+peakcoordinates(:,2); %x,y in pixels 
results(:,4)=Pcspline(:,5);
if p.bidirectional
    results(:,5)=results(:,4);
else
    results(:,5)=Pcspline(:,6);
end
results(:,6:7)=Pcspline(:,3:4);

results(:,8:9)=real(sqrt(CRLB(:,[2 1])));

results(:,10:11)=real(sqrt(CRLB(:,3:4)));
results(:,12)=LL;
end

end

function img=getimage(F,reader,p)
switch p.loader
     case 1
         img=reader.read(F);
    case 2
        img=bfGetPlane(reader,F);
    case 3        
        ss=[reader.getWidth reader.getHeight reader.getSize];
        if F>0&&F<=ss(3)
            pixel=reader.getPixels(F);
            img=reshape(pixel,ss(1),ss(2))';
        else
            img=[];
        end
end

end

function closereader(reader,p)
switch p.loader
     case 1
         reader.close;
    case 2
    case 3        

end

end