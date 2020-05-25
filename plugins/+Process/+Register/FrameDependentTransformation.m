classdef FrameDependentTransformation<interfaces.DialogProcessor
%     calculates transformation (global or local) based on localizations
%     (e.g. multi-color beads or fluorophores in ratiometric imaging)
    properties
        isz=0;
        transformation=[];
    end
    methods
        function obj=FrameDependentTransformation(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.showresults=true;
            obj.history=true;
        end
        
        function out=run(obj,p)          
            %implement: if not globally fitted: match locs
            % target coordinates.        
            % get coordinates: coordref, coordtarget in units of pixeles
            
            % global fit
            if isfield(obj.locData.loc,'xpix1')
                locs=obj.locData.getloc({'xpix1','xpix2','ypix1','ypix2','xpix1err','xpix2err','ypix1err','ypix2err','frame','filenumber'},'layer',find(obj.getPar('sr_layerson')),'Position','all');
                filenumber=locs.filenumber(1);
                Tinitial=obj.locData.files.file(filenumber).savefit.fitparameters.loc_globaltransform;
                roi=obj.locData.files.file(filenumber).info.roi;
                indbad=isnan(locs.xpix1)|isnan(locs.xpix2)|isnan(locs.ypix1)|isnan(locs.ypix2);
                coordref=horzcat(locs.xpix1(~indbad)+roi(1),locs.ypix1(~indbad)+roi(2));
                coordt2ref=horzcat(locs.xpix2(~indbad)+roi(1),locs.ypix2(~indbad)+roi(2));
                coordtarget=Tinitial.transformToTarget(2,coordt2ref,'pixel');
                frames=locs.frame(~indbad);
                
                %weights for fitting from localization precision. w~1/sigma^2            
                wx=1./(locs.xpix1err(~indbad).^2+locs.xpix2err(~indbad).^2);  %weights: It is variance! see Wikipedia
                wy=1./(locs.ypix1err(~indbad).^2+locs.ypix2err(~indbad).^2);
                allerror=((locs.xpix1err(~indbad)+(locs.ypix1err(~indbad))+(locs.xpix2err(~indbad))+(locs.ypix2err(~indbad)))); %average 1D loc prec
            else %single fit
                 
            
            %use initial transformation to put back coordinates 2 to
            %original position in chip
                locs=obj.locData.getloc({'xpix','ypix','xpixerr','ypixerr','frame','filenumber'},'layer',find(obj.getPar('sr_layerson')),'Position','all');
                filenumber=locs.filenumber(1);
                roi=obj.locData.files.file(filenumber).info.roi;
                %get Tinitial
                if isfield(obj.locData.files.file(filenumber),'transformationfit')
                    Tinitial=obj.locData.files.file(filenumber).transformationfit; %needed to transform global fitted data to target, this uses explicitely the transformation
                elseif isfield(obj.locData.files.file(filenumber),'transformation')
                    Tinitial=obj.locData.files.file(filenumber).transformation;
                else
                    errordlg('attach transformation')
    %                 Tinitial=obj.locData.files.file(1).savefit.fitparameters.
                end
                cpix=horzcat(locs.xpix+roi(1),locs.ypix+roi(2));
                indref=Tinitial.getPart(1,cpix); %get ref
                cref=cpix(indref,:);
                locsr.x=cref(:,1);locsr.y=cref(:,2);locsr.frame=locs.frame(indref);
                indtarget=Tinitial.getPart(2,cpix);%get target
                ctar=cpix(indtarget,:);
                ct2r=Tinitial.transformToReference(2,ctar);% transform target to ref
                locst.x=ct2r(:,1);locst.y=ct2r(:,2);locst.frame=locs.frame(indtarget);
                 %match
                [iAa,iBa,nA,nB,nseen]=matchlocsall(locsr,locst,0,0,1,inf);
                coordref=cref(iAa,:);
                coordtarget=ctar(iBa,:);
                
                %weights
                indrefi=find(indref);indtari=find(indtarget);
                wx=1./(locs.xpixerr(indrefi(iAa)).^2+locs.xpixerr(indtari(iBa)).^2);  %weights: It is variance! see Wikipedia
                wy=1./(locs.ypixerr(indrefi(iAa)).^2+locs.ypixerr(indtari(iBa)).^2);
                frames=locs.frame(indrefi(iAa));            
                allerror=(locs.xpixerr(indrefi(iAa))+locs.xpixerr(indtari(iBa))+locs.ypixerr(indrefi(iAa))+locs.ypixerr(indtari(iBa))); %average 1D loc prec
                
                %test
%                 fff=locs.frame(indrefi(iAa))-locs.frame(indtari(iBa)); is
%                 zero
%                 coordt2ref=ct2r(iBa,:);                
%                 r=1:100;
%                 figure(88);plot(coordref(r,1),coordref(r,2),'.',coordt2ref(r,1),coordt2ref(r,2),'.')
%                 %ref, target from match
            end


            ff=min(frames):max(frames);
            % initial shift
            coordt2refi=Tinitial.transformToReference(2,coordtarget,'pixel'); %put back to reference part
            dxdyi=coordt2refi-coordref; %differnece between reference and target localizations
            fitxi=getsmoothcurve(frames,dxdyi(:,1),wx); %get smoothed approximation of this difference
            fityi=getsmoothcurve(frames,dxdyi(:,2),wy); 
            axi=obj.initaxis('initial T');
            plot(axi,ff,fitxi(ff),ff,fityi(ff))
            
            
            %recalculate transformation based on localizations (might be
            %different from bead calibration)
            Trefine=Tinitial.copy;
            Trefine.findTransform(2,coordref,coordtarget,p.transform.selection,p.transformparam);
            
            coordt2ref2=Trefine.transformToReference(2,coordtarget,'pixel'); %put back to reference part

            dxdy=coordt2ref2-coordref; %differnece between reference and target localizations
            fitx=getsmoothcurve(frames,dxdy(:,1),wx); %get smoothed approximation of this difference
            fity=getsmoothcurve(frames,dxdy(:,2),wy);
            % plot shift vs frame:
            ax=obj.initaxis('shift');
            plot(ax,ff,fitx(ff),ff,fity(ff))
            
            % write shifts into transformation for testing
            Trefine.frameshift.fitx=fitx; %Trefine2 is calculated with these shifts.
            Trefine.frameshift.fity=fity;
            
            %  test
            coordrefcorr=Trefine.transformToReferenceFramecorrection(2,coordtarget,frames);
            dxdytest=coordrefcorr-coordref; %differnece between reference and target localizations
            fitx2=getsmoothcurve(frames,dxdytest(:,1),wx); %get smoothed approximation of this difference
            fity2=getsmoothcurve(frames,dxdytest(:,2),wy);
            ax2=obj.initaxis('test shift');
            plot(ax2,ff,fitx2(ff),ff,fity2(ff))            
            
            % SECOND ITERATION
            %manually correct shift to re-calculate Transformation based on
            %shifted coordinates
            mirrorfac=1-2*Tinitial.mirrorchannel(2);   %mirror of target. Assumption: reference is not mirrored.
            coordtargetcorr=coordtarget-horzcat(mirrorfac(1)*fitx(frames),mirrorfac(2)*fity(frames)); %subtract shift from target coordinates
            
            %filter out bad coordinates
            % Gaussian weighting
            sdist=mean(allerror)/4;
%             sdist=((mean(locs.xpix1err(~indbad))+mean(locs.ypix1err(~indbad))+mean(locs.xpix2err(~indbad))+mean(locs.ypix2err(~indbad))))/4; %average 1D loc prec
            %keep localizations that are closer together than 2*sigma_av
            %and have localization precision < 2*sigma_av
            indgood=sum(dxdytest.^2,2)<(sdist*2)^2 & allerror < sdist*2 ;

%             indgood=sum(dxdytest.^2,2)<(sdist*2)^2 & ((locs.xpix1err(~indbad))+(locs.ypix1err(~indbad))+(locs.xpix2err(~indbad))+(locs.ypix2err(~indbad)))< sdist*2 ;
            %re-calculate transformation
            Trefine2=Trefine.copy; %redo transformation based on origingal shift
            Trefine2.findTransform(2,coordref(indgood,:),coordtargetcorr(indgood,:),p.transform.selection,p.transformparam);
            
            %test forward and backward transform with correction
            coordrefcorr=Trefine2.transformToReferenceFramecorrection(2,coordtarget,frames);
            dxdytest=coordrefcorr-coordref; %differnece between reference and target localizations
            fitxr=getsmoothcurve(frames,dxdytest(:,1),wx); %get smoothed approximation of this difference
            fityr=getsmoothcurve(frames,dxdytest(:,2),wy);
            axr=obj.initaxis('2nd iteration to reference');
            plot(axr,ff,fitxr(ff),ff,fityr(ff))  
            
            coordref2tcorr=Trefine2.transformToTargetFramecorrection(2,coordref,frames);
            dxdytest4=coordref2tcorr-coordtarget; %differnece between reference and target localizations
            fitx4=getsmoothcurve(frames,dxdytest4(:,1),wx); %get smoothed approximation of this difference
            fity4=getsmoothcurve(frames,dxdytest4(:,2),wy);
            ax4=obj.initaxis('2nd iteration to target');
            plot(ax4,ff,fitx4(ff),ff,fity4(ff))
            
            % save and update
            transformation=Trefine2;
            fout=strrep(obj.locData.files.file(1).name,'_sml.mat','_T.mat');
            save(fout,'transformation');
            obj.setPar('transformationfile',fout);
            if ~isfield(obj.locData.files.file(1),'transformationfit')
                obj.locData.files.file(1).transformationfit=obj.locData.files.file(1).transformation;
            end
            obj.locData.files.file(1).transformation=transformation;
            if p.updatesmlfile  
                obj.locData.savelocs(obj.locData.files.file(1).name);
            end

            %for testing: write dx, dy for all localizations to
            %localization data
            if p.writetoloc
                coordallr=horzcat(obj.locData.loc.xpix1+roi(1),obj.locData.loc.ypix1+roi(2));

                coordalltar=horzcat(obj.locData.loc.xpix2+roi(1),obj.locData.loc.ypix2+roi(2));
                coordallt=Tinitial.transformToTarget(2,coordalltar,'pixel');
                coordalltc=coordallt-horzcat(mirrorfac(1)*fitx(obj.locData.loc.frame),mirrorfac(2)*fity(obj.locData.loc.frame));
                coordallt2ref=Trefine2.transformToReference(2,coordalltc,'pixel');
                dxdyall=coordallt2ref-coordallr;
                obj.locData.loc.dxfcorr=dxdyall(:,1);
                obj.locData.loc.dyfcorr=dxdyall(:,2);
                obj.locData.regroup;
            end
            out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function initGui(obj)
        end
    end
end

function fitx=getsmoothcurve(frames,dx,wx)           
dframe=100;
ff=(min(frames): dframe:max(frames))';
if nargin<3
    wx=[];
end
[dxb,sb]=bindatamean(frames,dx,ff,wx);
[fitx,gof,p]=fit(ff,dxb,'smoothingspline','Weights',sb,'SmoothingParam',1e-11*dframe);
end
            
            

function pard=guidef(obj)
% pard.Tfile.object=struct('Style','edit','String','');
% pard.Tfile.position=[8,1];
% pard.Tfile.Width=3;
% pard.Tfile.object.TooltipString='default file for transformation matrix. You can select new file after transformation has been calculated.';
% 
% pard.browse.object=struct('Style','pushbutton','String','load T','Callback',@obj.browse_callback);
% pard.browse.position=[8,4];
% pard.browse.object.TooltipString='Save the newly calculated transformation matrix.';


pard.texttt.object=struct('String','Transformation:','Style','text');
pard.texttt.position=[1,1];
pard.transform.object=struct('Style','popupmenu','String','projective|affine|similarity|polynomial|lwm|pwl');
pard.transform.position=[1,2];
pard.transform.object.TooltipString='select one of Matlabs transformations. Not all might work.';

pard.transformparam.object=struct('Style','edit','String','3');
pard.transformparam.position=[1,3];
pard.transformparam.object.TooltipString='Parameter for lwm and polynomial';
pard.transformparam.Width=0.5;

pard.updatesmlfile.object=struct('Style','checkbox','String','write T to .sml','Value',0);
pard.updatesmlfile.position=[2,1];
pard.updatesmlfile.object.TooltipString='If checked, the transformation file is appended to the .sml file and saved there as well when you click save T';

pard.writetoloc.object=struct('Style','checkbox','String','write dx,dy to locData.loc','Value',0);
pard.writetoloc.position=[3,1];
pard.writetoloc.object.TooltipString='If checked, the transformation file is appended to the .sml file and saved there as well when you click save T';
pard.writetoloc.Width=3;

pard.inputParameters={'currentfileinfo'};
pard.plugininfo.description='calculates transformation (global or local) based on localizations (e.g. multi-color beads or fluorophores in ratiometric imaging)';
pard.plugininfo.type='ProcessorPlugin';
end