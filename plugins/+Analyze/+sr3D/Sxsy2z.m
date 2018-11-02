classdef Sxsy2z<interfaces.DialogProcessor
    properties
%         cal3D
%         outsx2sy2
%         spline
    end
    methods
        function obj=Sxsy2z(varargin)        
           obj@interfaces.DialogProcessor(varargin{:}) ;
             obj.inputParameters={'cam_pixelsize_nm'};
           obj.showresults=false;
            obj.history=true;
        end
        function initGui(obj)
            set(obj.guihandles.savebutton,'Callback',{@loadcalfile,obj})
        end
        function out=run(obj,p)
            out=[];
            locs=obj.locData.getloc({'PSFxnm','PSFynm'});
            if isempty(locs.PSFynm)
                error('PSFynm required')
            end
            
            ld=load(p.calfile);
            
            if ~p.usespline
                sx2sy2=(locs.PSFxnm/p.cam_pixelsize_nm(1)).^2-(locs.PSFynm/p.cam_pixelsize_nm(end)).^2;
                z=feval(ld.outsx2sy2,sx2sy2);
                z(locs.PSFynm==0)=-2;
            else %spline method
                spline=ld.outspline.s2z;
                sx=locs.PSFxnm/p.cam_pixelsize_nm;
                sy=locs.PSFynm/p.cam_pixelsize_nm;
                
                sxr=round((sx-spline.smin)/spline.ds);
                syr=round((sy-spline.smin)/spline.ds);
                
                ns=length(spline.s);
                goodind=~(sxr<1|sxr>ns|syr<1|syr>ns);
                
                sind=sub2ind(size(spline.z),sxr(goodind),syr(goodind));
                z=zeros(size(locs.PSFxnm),'single');
                z(goodind)=-spline.z(sind);
                z(~goodind)=-2;
                d=zeros(size(locs.PSFxnm),'single');
                d(goodind)=spline.d(sind);
                d(~goodind)=inf;
                obj.locData.loc.zfromPSFdistance=d;
                %quite arbitratry locprec.
                obj.locData.loc.locprecnmz=obj.locData.loc.locprecnm*2.5.*(1+10*d);
            end
%                 smin=0;
%                 smax=3;
%                 ds=0.03;
%                 s=smin:ds:smax;
%                 zmat=zeros(length(s));
%                 dmat=zeros(length(s));
%                 for k=1:length(s)
%                     k
%                     for l=1:length(s)
%                         [zmat(l,k),dmat(l,k)]=zfromspline(s(k),s(l),obj.spline.x,obj.spline.y);
%                     end
%                 end
%             end
%             z=-1:0.03:1;
%             
%             
%             figure(90);hold off
%             imagesc(s,s,dmat)
%             colorbar
%            hold on;
%            plot(obj.spline.x(z),obj.spline.y(z))
%             
%             figure(91);
%              hold off
%            imagesc(s,s,zmat)
%            colorbar
%            hold on;
%            plot(obj.spline.x(z),obj.spline.y(z))
%            
            
            
%             refractiveIndexMismatch=1.25
            obj.locData.loc.znm=z*1000*p.refractiveIndexMismatch;
            obj.locData.regroup;
            obj.setPar('locFields',fieldnames(obj.locData.loc))
     
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end

% function [z,d]=zfromspline(sx,sy,splinex,spliney)
% % zt=[-.5 0 .5];
% % zt=0;
% % for k=length(zt):-1:1
% %     zh(k)=fminsearch(@distsplinerr,zt(k),[],sx,sy,splinex,spliney);
% %     dh(k)=distspline(zh(k),sx,sy,splinex,spliney);
% % end
% % [d,indm]=min(dh);
% % z=zh(indm);
% 
% z=fminbnd(@distsplinerr,-1,1,[],sx,sy,splinex,spliney);
%  d=distspline(z,sx,sy,splinex,spliney);
% end
% 
% function d=distspline(z,sx,sy,splinex,spliney)
% d=sqrt(((sx)-(splinex(z))).^2+((sy)-(spliney(z))).^2);
% % d=((sqrt(sx)-sqrt(splinex(z))).^2+(sqrt(sy)-sqrt(spliney(z))).^2);
% end
% function err= distsplinerr(z,sx,sy,splinex,spliney)
% err=sum(distspline(z,sx,sy,splinex,spliney));
% end
function loadcalfile(a,b,obj)
fn=obj.guihandles.calfile.String;
[f,p]=uigetfile(fn);
if f
%     load([p f])
%     if exist('cal3D','var')
%         obj.cal3D=cal3D;
%     end
%     if exist('outspline','var')
%         obj.spline=outspline;
%     end
%     obj.outsx2sy2=outsx2sy2;
    obj.guihandles.calfile.String=[p f];
end
end

function pard=guidef
% pard.d3_color.object=struct('Style','checkbox','String','render in color');
% pard.d3_color.position=[2,1];

pard.zposobjective_check.object=struct('String','z-position (um) of objective','Style','checkbox');
pard.zposobjective_check.position=[1,1];
pard.zposobjective_check.Width=2;


pard.zposobjective.object=struct('Style','edit','String','30'); 
pard.zposobjective.position=[1,3];

pard.pol_check.object=struct('String','give polynome here:','Style','checkbox');
pard.pol_check.position=[2,1];
pard.pol_check.Width=2;


pard.fitpol.object=struct('Style','edit','String','1 0 0'); 
pard.fitpol.position=[2,3];


pard.t1.object=struct('Style','text','String','refractive Index Mismatch factor (<=1)'); 
pard.t1.position=[3,1];
pard.t1.Width=2;
pard.refractiveIndexMismatch.object=struct('Style','edit','String','1'); 
pard.refractiveIndexMismatch.position=[3,3];

pard.calfile.object=struct('Style','edit','String','settings/cal_3DAcal.mat');
pard.calfile.position=[4,1];
pard.calfile.Width=3;


pard.savebutton.object=struct('String','load','Style','pushbutton');
pard.savebutton.position=[4,4];

pard.usespline.object=struct('String','Use spline approximation','Style','checkbox');
pard.usespline.position=[6,1];
pard.usespline.Width=2;

pard.syncParameters={{'cal3D_file','calfile',{'String'}}};

pard.plugininfo.type='ProcessorPlugin';

end
