classdef Phase2z4Pi<interfaces.DialogProcessor
    % calculates z from size of PSF in x or y / channel 1/2. For biplane
    % and astigmatic 3D.
%     Related to: Huang, Fang, George Sirinakis, Edward S. Allgeyer, Lena
%     K. Schroeder, Whitney C. Duim, Emil B. Kromann, Thomy Phan, et al.
%     “Ultra-High Resolution 3D Imaging of Whole Cells.” Cell 166, no. 4
%     (August 2016): 1028–40. https://doi.org/10.1016/j.cell.2016.06.016.

    methods
        function obj=Phase2z4Pi(varargin)        
            obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.showresults=true;
             obj.history=true;

        end
        
        function out=run(obj,p)
            obj.setPar('undoModule','Phase2z4Pi');
            notify(obj.P,'backup4undo');
%             fitterpath=[fileparts(obj.getPar('maindirectory')) filesep 'ries-private' filesep 'PSF4Pi'];
%             if ~isdeployed
%                 addpath(fitterpath)
%             end
            locsall=obj.locData.getloc({'znm','zastig','zastigerr','phase','znmerr','phaseerr','frame'});
            locs=obj.locData.getloc({'znm','zastig','zastigerr','phase','znmerr','phaseerr','frame','filenumber'},'layer',find(obj.getPar('sr_layerson')),'position','fov');
            if isempty(locs.zastig)
                zastig=locs.znm;
                zastigerr=locs.znmerr;
                zastigall=locsall.znm;
                zastigerrall=locsall.znmerr;
            else
                zastig=locs.zastig;
                zastigerr=locs.zastigerr;
                zastigall=locsall.zastig;
                zastigerrall=locsall.zastigerr;
            end
            zastig=-zastig;zastigall=-zastigall; %in the fitter z inverted
            phase=mod(locs.phase,2*pi);
            phaseerr=locs.phaseerr;
            phaseerrall = locsall.phaseerr;
            cal3D=obj.locData.files.file(locs.filenumber(1)).savefit.cal3D;
%             frequency=1/cal3D.zT/cal3D.pixelsize_z;
            periodnm=cal3D.zT*1000;
            frequency=pi/periodnm;
%             numwindows=1; windowsize=ceil(max(locs.frame)/numwindows);
            framepos=[0:p.windowframe:max(locs.frame) max(locs.frame)+1];
            if framepos(end)-framepos(end-1)<(framepos(2)-framepos(1))*0.75
                framepos(end-1)=[];
            end
            
            z0=0;
            axp=obj.initaxis('phase vs z');
            axph=axp;
            %ax2=obj.initaxis('z0 corrected phase');

            %hold(ax2,'off')
            for k=1:length(framepos)-1
                inframe=locs.frame>framepos(k)&locs.frame<framepos(k+1);
                z0=getz0phase(zastig(inframe),phase(inframe),frequency,z0,axph);
                axph=[];
                z0all(k)=z0;
                [zph,phicor]=z_from_phi_JR(zastig(inframe),mod(phase(inframe),2*pi),frequency,z0);
%                 [zph,phicor]=z_from_phi_JR(zastig,mod(phase,2*pi),frequency,z0);
                %plot(ax2,zastig(inframe),mod(phicor,2*pi),'.','markersize',1)
                 %hold(ax2,'on')
                
            end
            
            
%             z0=getz0phase(zastig(inframe),phase(inframe),frequency,z0,axp);
            
            frameposc=framepos(1:end-1)+(framepos(2)-framepos(1))/2;
            if length(z0all)>3
                z0int=fit(double(frameposc)',double(z0all)','smoothingspline');
                ax=obj.initaxis('z0');
                plot(ax,frameposc,z0all,frameposc,z0int(frameposc))
                z0=z0int(locsall.frame);
            else
                z0=mean(z0all);
            end
            [zph,phcor]=z_from_phi_JR(zastigall,mod(locsall.phase,2*pi),frequency,z0); %minus from fitter inversion
            obj.locData.setloc('zphase',-zph);
            if isempty(locs.zastig)
                obj.locData.setloc('zastig',zastigall);
                obj.locData.setloc('zastigerr',zastigerrall);
            end
            ax3=obj.initaxis('phase corr all');
            dscatter(zastigall,phcor)
            hold on
            zplot=min(zastigall):max(zastigall);
        
            plot(zplot,mod(zplot*2*frequency,2*pi),'k.')
            hold off
            
            obj.locData.setloc('znm',-(zph+z0));
            obj.locData.setloc('znmerr',phaseerrall/2/frequency);
            obj.locData.setloc('locprecznm',phaseerrall/2/frequency);
            obj.locData.regroup;
            obj.locData.filter;
            out=0;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

function pard=guidef
pard.windowframet.object=struct('String','window (frames)','Style','text');
pard.windowframet.position=[1,1];

pard.windowframe.object=struct('String','2000','Style','edit');
pard.windowframe.position=[1,2];

% 
% pard.text2.object=struct('String','zmin','Style','text');
% pard.text2.position=[2,1];
% pard.text3.object=struct('String','zmax','Style','text');
% pard.text3.position=[3,1];
% 
% pard.zmin.object=struct('Style','edit','String',-400); 
% pard.zmin.position=[2,2];
% pard.zmax.object=struct('Style','edit','String',400); 
% pard.zmax.position=[3,2];
% 
% pard.pixauto.object=struct('Style','checkbox','String','set pixelsize (x,z)','Value',0);
% pard.pixauto.position=[4,1];
% pard.pixrecset.object=struct('Style','edit','String','5, 5'); 
% pard.pixrecset.position=[4,2];

% pard.plugininfo.description= 'Side view from ROI';
pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='calculates z from size of PSF in x or y / channel 1/2. For biplane and astigmatic 3D. Related to: Huang, Fang, George Sirinakis, Edward S. Allgeyer, Lena K. Schroeder, Whitney C. Duim, Emil B. Kromann, Thomy Phan, et al. “Ultra-High Resolution 3D Imaging of Whole Cells.” Cell 166, no. 4 (August 2016): 1028–40. https://doi.org/10.1016/j.cell.2016.06.016.';
end

function z0=getz0phase(zastig,phase,frequency,z0,ax)
% phasez=mod((zastig-z0)*2*frequency,2*pi);
% cyclicaverage(mod(phase-phasez,2*pi),2*pi)
% z0=0
% frequency=pi/periodnm;
zfp=phase/2/frequency;
dz=zfp-zastig+pi/frequency/2+z0;
dzm=mod(dz,pi/frequency);
z0=-cyclicaverage(dzm,pi/frequency)+pi/frequency/2+z0;



if nargin>4 &&~isempty(ax)
    phasez=mod((zastig-z0)*2*frequency,2*pi);
    indgood=zastig>min(zastig)+1 & zastig<max(zastig)-1 ;
    axes(ax)
    dscatter(zastig(indgood),phase(indgood))
    hold on
    plot(zastig,phasez,'r.')
%     plot(zastig(indgood),zph(indgood),'g.')
    xlabel('zastig')
    ylabel('phase')
    legend('phase','phase from z astig')
%     figure(88)
%     plot(-zastig(indgood),zph(indgood),'k.')
end

% figure(88);histogram(dzm)
% waitforbuttonpress
% phasez=mod((zastig-z0)*2*frequency,2*pi);
% dphase=phase-phasez;
% if sum(dphase>0)>sum(dphase<0)
%     dphasem=mean(dphase(dphase>0));
%     phx=-pi/frequency;
% else
%     dphasem=mean(dphase(dphase<0));
%     phx=-pi/frequency;
% end
% dz=dphasem/2/frequency;
% z0=z0+dz-phx;
% 
% phasez1b=mod((zastig-z0)*2*frequency,2*pi);

% iter=10;
% err=1;
% for k=1:iter
% phasez2=mod((zastig-z0)*2*frequency,2*pi);
% dphasez2=phase-phasez2;
% dz=mean(dphasez2)/2/frequency;
% z0=z0+dz;
% if abs(dz)<err
%     break
% end%roubst mean later?
% 
% end
% 
% phasez2=mod((zastig-z0)*2*frequency,2*pi);
% figure(88);plot(zastig,phase,'.',zastig,phasez2,'+')
% % % figure(88);plot(zastig,phase,'.',zastig,phasez,'+',zastig,phasez2,'*',zastig,phasez1b,'x')
% waitforbuttonpress
% f=@(z0,x) mod((x-z0)*2*frequency,2*pi);
% fp=fit(zastig,phase,f,'StartPoint',z0);
% phasez2=mod((zastig-fp.z0)*2*frequency,2*pi);
% figure(88);plot(zastig,phase,'.',zastig,phasez2,'+')
% waitforbuttonpress
% z0=fp.z0;
end
