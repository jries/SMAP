%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
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
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.
%%
function [SXY,beadpos,parameters]=calibrate3D_gs(p)
% p.filelist
% p.outputfile
% p.dz
% p.modality
% p.zcorr
% p.ROIxy
% p.ROIz
% p.smoothxy
% p.smoothz
% p.gaussrange
% p.filter;
% p.zcorrframes
% p.gaussroi
% p.fov =[x1 y1 x2 y2] 
% p.mindistance
% p.emgain
% p.xrange
% p.yrange
% p.positions (get beads only here)
% p.smap (called from SMAP: extended functionality)
% isglobalfit
% Tfile

%initialize paramters that were not transferred:
if ~isfield(p,'smap')&& (~p.isglobalfit)
    p.smap=false;
end

if ~isfield(p,'xrange')
    p.xrange=[-inf inf]; p.yrange=[-inf inf]; 
end
if ~isfield(p,'emmirror')
    p.emmirror=0;
end

if ~isfield(p,'smoothxy')
    p.smoothxy=0;
end

if ~isfield(p,'isglobalfit')
    p.isglobalfit=0;
end
if ~isfield(p,'transformation')
    p.transformation=[];
end
if ~isfield(p,'filechannel')
    p.filechannel=1;
end

%get bead positions
p.status.String='Load files and segment beads';drawnow

%determine out put figures
if ~isfield(p,'tabgroup')
    f=figure('Name','Bead calibration');
    p.tabgroup=uitabgroup(f);
    calibrationfigure=f;
else
    calibrationfigure=p.tabgroup.Parent;
end

%get beads from images
[beads,p]=images2beads_globalfit(p);

%if only to take beads in a certain range, remove others
imageRoi=p.roi{1};

if isempty(beads)
    warndlg('Could not find and segment any bead. ROI size too large?')
    p.status.String='error: could not find and segment any bead...';
    return
end

p.midpoint=round(size(beads(1).stack.image,3)/2); %reference for beads
p.ploton=false;

%write reference frame to beads:
% for k=1:length(beads)
%     beads(k).f0=p.midpoint;
% end

%get positions of beads
for k=length(beads):-1:1
    beadposxs(k)=beads(k).pos(1);
    beadposys(k)=beads(k).pos(2);
    beadfilenumber(k)=beads(k).filenumber;
end

% spatially dependent calibration. Used e.g. for two separae channels, but
% also for field-depndent PSF calibration.
tgmain=p.tabgroup;
for X=1:length(p.xrange)-1
    for Y=1:length(p.yrange)-1
        if length(p.xrange)>2||length(p.yrange)>2
            ht=uitab(tgmain,'Title',['X' num2str(X) 'Y' num2str(Y)]);
            p.tabgroup=uitabgroup(ht);
        end
        
        indgood=beadposxs< p.xrange(X+1) & beadposxs>p.xrange(X) & beadposys<p.yrange(Y+1) & beadposys>p.yrange(Y);
        beadsh=beads(indgood);
        
        for k=1:max(beadfilenumber)
            indfile=(beadfilenumber==k)&indgood;
            p.fileax(k).NextPlot='add';
            scatter(p.fileax(k),beadposxs(indfile),beadposys(indfile),60,[1 1 1])
            scatter(p.fileax(k),beadposxs(indfile),beadposys(indfile),50)
        end
        if isempty(beadsh)
            disp(['no beads found in part' num2str(p.xrange(X:X+1)) ', ' num2str(p.yrange(Y:Y+1))])
            continue
        end

        % get cspline calibration
        p.status.String='get cspline calibration';drawnow
        [csplinecal,indgoods,beadpos{X,Y},~,testallrois]=getstackcal_gs(beadsh,p); %calculate average PSF
        
        if ~isempty(beadpos{X,Y})
            for f=1:max(beadpos{X,Y}.filenumber(:))
                indfile=(beadpos{X,Y}.filenumber==f);
                p.fileax(f).NextPlot='add';
                plot(p.fileax(f),beadpos{X,Y}.xim(indfile),beadpos{X,Y}.yim(indfile),'m+')
            end
        end
        
        cspline=csplinecal.cspline;
        cspline.global.isglobal=p.isglobalfit;
        cspline.global.transformation=p.transformation;
        if p.isglobalfit
            cspline.global.coeffrawref=single(csplinecal.cspline.coeffrawref);
            cspline.global.coeffrawtar=single(csplinecal.cspline.coeffrawtar);
            cspline.normf=csplinecal.cspline.normf;
        end
        
        PSF=csplinecal.PSF;
        SXY(X,Y)=struct('cspline',cspline,'Xrangeall',p.xrange+imageRoi(1),'Yrangeall',p.yrange+imageRoi(2),'Xrange',p.xrange([X X+1])+imageRoi(1),...
            'Yrange',p.yrange([Y Y+1])+imageRoi(2),'posind',[X,Y],'EMon',p.emgain,'PSF',{PSF});
        % ZERNIKE fitting taken out, now in calibrate3D_g.m
    end
end
axcrlb=axes(uitab(p.tabgroup,'Title','CRLB'));
plotCRLBcsplinePSF(csplinecal.cspline,axcrlb)

parameters=myrmfield(p,{'tabgroup','status','ax_z','ax_sxsy','fileax'});
    
p.status.String='save calibration';drawnow

if ~isempty(p.outputfile)
    if p.smap
        parameters.smappos.P=[];
        save(p.outputfile,'SXY','parameters');
    else
        save(p.outputfile,'cspline','parameters');
    end
    filefig=strrep(p.outputfile,'.mat','.fig');
    savefig(calibrationfigure,filefig,'compact');
end
p.status.String='Calibration done';drawnow
end




