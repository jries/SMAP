classdef Register3Dstacks<interfaces.DialogProcessor
    properties

    end
    methods
        function obj=Register3Dstacks(varargin)        
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.history=true;
                obj.showresults=true;
                obj.guiselector.show=true;

        end
        
        function out=run(obj,p)
            obj.setPar('undoModule','Register3Dstacks');
            notify(obj.P,'backup4undo');
            out=[];
            locs=obj.locData.getloc({'xnm','ynm','znm','filenumber'},'layer',find(obj.getPar('sr_layerson')),'position','roi','removeFilter',{'filenumber'});
            numfiles=max(locs.filenumber);
            p.axxy=obj.initaxis('shift xy');
            p.axz=obj.initaxis('shift z');
            
%             zslice=0*locs.znm;
            for k=1:numfiles
                [~,filename]=fileparts(obj.locData.files.file(k).name);
                filename=strrep(filename,'.','_');
                ind1=strfind(filename,'L');
                ind2=strfind(filename,'S');
                ind3=strfind(filename,'_');
                L(k)=str2double(filename(ind1(1)+1:ind2(1)-1));
                S(k)=str2double(filename(ind2(1)+1:ind3(1)-1));
%                 ind=strfind(filename,'_Pos');
%                 ind2=strfind(filename(ind+4:end),'_');
%                 zpos(k)=str2double(filename(ind+4:ind+2+ind2(1)));
%                 zslice(=locs.znm+zpos(k)*p.dz;
            end
            filenumbers=1:numfiles;
            slices=min(S):max(S);
            
            if ~p.onlynominal
                for ks=1:length(slices)
                    sh=slices(ks);
                    indslice=find(S==sh);

                    coordall=[locs.xnm,locs.ynm,locs.znm];
                    in1=locs.filenumber==filenumbers(indslice(1));
                    numlocsall(ks,1)=sum(in1);
    %                 figure(100);hold off
    %                 plot(coordall(in1,1),coordall(in1,2),'.'); hold on
                    for kl=2:length(indslice)                
                        in2=locs.filenumber==filenumbers(indslice(kl));
                        numlocsall(ks,kl)=sum(in2);
                        coordref=coordall(in1,:);
                        coordtar=coordall(in2,:);
                        shift(kl,:,ks)=findshift(coordref,coordtar,p); %shift between two consecutive frames
                        coordall(in2,:)=coordall(in2,:)+shift(kl,:,ks); %%XXXXXXXXX
    %                     figure(100); plot(coordall(in2,1),coordall(in2,2),'.'); hold on
                        in1=in1|in2; %correlate with all previous aligned localizations;
                        endgit 
                    inslice{ks}=in1;
                end
                shift
                
                axn=obj.initaxis('localizations');plot(axn,numlocsall);
    %             [zpossort,filesortind]=sort(zpos);
                for k=1:length(slices)-1

                    in1f=inslice{k};
                    in2f=inslice{k+1};
                    inabove=(locs.znm>p.dz/2-p.overlap/2)& locs.znm<p.dz/2+p.overlap/2;
                    inbelow=(locs.znm<-p.dz/2+p.overlap/2) & (locs.znm>-p.dz/2-p.overlap/2);
                    in1=in1f&inabove;
                    in2=in2f&inbelow;
                    coordref=coordall(in1,:);coordref(:,3)=coordref(:,3)+k*p.dz;
                    coordtar=coordall(in2,:);coordtar(:,3)=coordtar(:,3)+(k+1)*p.dz;
    %                 coordref=[locs.xnm(in1),locs.ynm(in1),locs.znm(in1)+k*p.dz];
    %                 coordtar=[locs.xnm(in2),locs.ynm(in2),locs.znm(in2)+(k+1)*p.dz];
                    shiftsl(k,:)=findshift(coordref,coordtar,p); %shift between two consecutive frames
                end  
                shiftsl
                shiftslc=cumsum(shiftsl,1);
            end
            for ks=1:length(slices)
                sh=slices(ks);
                indslice=find(S==sh);
                for kl=1:length(indslice) 
                    if p.onlynominal
                        shifth=[0 0 0];
                    else
                        shifth=squeeze(shift(kl,:,ks))+squeeze(shiftslc(ks));
                    end
                    in2f=obj.locData.loc.filenumber==filenumbers(indslice(kl));
                    obj.locData.loc.xnm(in2f)=obj.locData.loc.xnm(in2f)+shifth(1); %instead we need to use cumulative shifts here!
                    obj.locData.loc.ynm(in2f)=obj.locData.loc.ynm(in2f)+shifth(2);
                    obj.locData.loc.znm(in2f)=obj.locData.loc.znm(in2f)+shifth(3)+slices(ks)*p.dz;
                end
            end
            
            
%             [zpossort,filesortind]=sort(zpos);
%             for k=1:numfiles-1
%                 in1f=(locs.filenumber==filesortind(k));
%                 in2f=(locs.filenumber==filesortind(k+1));
%                 inabove=(locs.znm>p.dz/2-p.overlap/2)& locs.znm<p.dz/2+p.overlap/2;
%                 inbelow=(locs.znm<-p.dz/2+p.overlap/2) & (locs.znm>-p.dz/2-p.overlap/2);
%                 in1=in1f&inabove;
%                 in2=in2f&inbelow;
%                 
%                 coordref=[locs.xnm(in1),locs.ynm(in1),locs.znm(in1)+zpossort(k)*p.dz];
%                 coordtar=[locs.xnm(in2),locs.ynm(in2),locs.znm(in2)+zpossort(k+1)*p.dz];
%                 shift(k,:)=findshift(coordref,coordtar,p); %shift between two consecutive frames
%             end
%             shiftc=cumsum(shift,1);
%             in2f=obj.locData.loc.filenumber==filesortind(1);
%             obj.locData.loc.znm(in2f)=obj.locData.loc.znm(in2f)+zpossort(1)*p.dz;
%             for k=1:numfiles-1
%                 in2f=obj.locData.loc.filenumber==filesortind(k+1);
%                 obj.locData.loc.xnm(in2f)=obj.locData.loc.xnm(in2f)+shiftc(k,1); %instead we need to use cumulative shifts here!
%                 obj.locData.loc.ynm(in2f)=obj.locData.loc.ynm(in2f)+shiftc(k,2);
%                 obj.locData.loc.znm(in2f)=obj.locData.loc.znm(in2f)+zpossort(k+1)*p.dz+shiftc(k,3);
%             end
            obj.locData.regroup;
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end

function shift=findshift(coordref,coordtar,p)
pixrec=10;
wind=6;
ploton=true;
maxshift=10*wind;
slicewidth=150 ; %nm
pixrecz=3;
cmin=min(min(coordref,[],1,'omitnan'),min(coordtar,[],1,'omitnan'));
cmax=max(max(coordref,[],1,'omitnan'),max(coordtar,[],1,'omitnan'));
rangex=cmin(1):pixrec:cmax(1);
rangey=cmin(2):pixrec:cmax(2);
slicex=cmin(1):slicewidth:cmax(1);
slicey=cmin(2):slicewidth:cmax(2);

zbin=cmin(3):pixrecz:cmax(3);

hr=histcounts2(coordref(:,1),coordref(:,2),rangex,rangey);
ht=histcounts2(coordtar(:,1),coordtar(:,2),rangex,rangey);

if isfield(p,'axxy')
    drawnow
    ploton=p.axxy;
end
% f=figure(99);
[dx,dy,abg]=getShiftCorr(hr,ht,ploton,maxshift,true,wind);
% figure(98);plotaxis=gca;
coordtarc=coordtar+[dx dy 0]*pixrec;
zpos=finddisplacementZ2(coordref,coordtarc,slicex,slicey,zbin,5,p.axz);

shift=[dx.*pixrec,dy.*pixrec,zpos];
end

function pard=guidef(obj)


% p(1).value=0; p(1).on={}; p(1).off={'texta','drift_timepoints','text1','drift_pixrec','text2','drift_window','text3','drift_maxdrift','drift_maxpixelst','drift_maxpixels'};
% p(2).value=1; p(2).on=p(1).off; p(2).off={};
%             
% pard.correctxy.object=struct('String','Correct xy-drift','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
% pard.correctxy.position=[1,1];
% 
pard.dzt.object=struct('String','z-distance (nm)','Style','text');
pard.dzt.position=[1,1];
pard.dzt.Width=1;
% 
pard.dz.object=struct('String','300','Style','edit');
% pard.drift_timepoints.object.TooltipString=sprintf('whole data is divided into timepoints individual \n blocks. Range: 7-20');
pard.dz.position=[1,2];
pard.dz.Width=0.5;
% 
pard.overlapt.object=struct('String','overlap (nm)','Style','text');
pard.overlapt.position=[2,1];
pard.overlapt.Width=1;
pard.overlap.object=struct('String','150','Style','edit');
pard.overlap.position=[2,2];
pard.overlap.Width=0.5;


pard.onlynominal.object=struct('String','Only shift by nominal distance','Style','checkbox');
pard.onlynominal.position=[1,3];
pard.onlynominal.Optional=true;
pard.onlynominal.Width=2;
% 
% pard.drift_pixrec.object=struct('String','10','Style','edit');
% pard.drift_pixrec.position=[2,2.65];
% % pard.drift_pixrec.isnumeric=1;
% pard.drift_pixrec.object.TooltipString=sprintf('pixel size (nm) for reconstruction. \n Smaller for well defined peak. But slower \n Range: 10-25');
% pard.drift_pixrec.Optional=true;
% pard.drift_pixrec.Width=0.25;
% 
% pard.text2.object=struct('String','window pix','Style','text');
% pard.text2.position=[3,1];
% pard.text2.Optional=true;
% pard.text2.Width=0.75;
% 
% pard.drift_window.object=struct('String','7','Style','edit');
% pard.drift_window.position=[3,1.65];
% % pard.drift_window.isnumeric=1;
% pard.drift_window.object.TooltipString=sprintf('size of region for peakfinding (ellipt. Gaussian). \n should be small, but cover clear maximum. \n Range: 7-15');
% pard.drift_window.Optional=true;
% pard.drift_window.Width=0.25;
% 
% pard.text3.object=struct('String','maxdrift nm','Style','text');
% pard.text3.position=[4,1];
% pard.text3.Optional=true;
% 
% pard.drift_maxdrift.object=struct('String','1000','Style','edit');
% pard.drift_maxdrift.position=[4,2];
% pard.drift_maxdrift.Width=0.9;
% 
% pard.drift_maxdrift.object.TooltipString=sprintf('Maximum drift expected. \n Smaller if data is sparse and wrong peak found. \n larger if no clear peak found. \n Range 250-2000');
% pard.drift_maxdrift.Optional=true;
% 
% pard.drift_maxpixelst.object=struct('String','max size (pix)','Style','text');
% pard.drift_maxpixelst.position=[5,1];
% pard.drift_maxpixelst.Optional=true;
% 
% pard.drift_maxpixels.object=struct('String','4096','Style','edit');
% pard.drift_maxpixels.position=[5,2];
% pard.drift_maxpixels.object.TooltipString=sprintf('Maximum size of the reconstructed images. Smaller for speed and lower memory consumption, larger for noisy signal. 128-4096');
% pard.drift_maxpixels.Optional=true;
% pard.drift_maxpixels.Width=0.9;
% 
% p(1).value=0; p(1).on={}; p(1).off={'textaz','drift_timepointsz','drift_pixreczt','drift_pixrecz','drift_windowzt','drift_windowz','zranget','zrange','slicewidtht','slicewidth'};
% p(2).value=1; p(2).on=p(1).off; p(2).off={};
% pard.correctz.object=struct('String','Correct z-drift','Style','checkbox','Value',0,'Callback',{{@obj.switchvisible,p}});
% pard.correctz.position=[1,3];
% 
% pard.textaz.object=struct('String','timepoints z','Style','text','Visible','off');
% pard.textaz.position=[2,3];
% pard.textaz.Optional=true;
% pard.textaz.Width=.75;
% 
% pard.drift_timepointsz.object=struct('String','10','Style','edit','Visible','off');
% pard.drift_timepointsz.object.TooltipString=sprintf('whole data is divided into timepoints individual \n blocks. Range: 10-40');
% pard.drift_timepointsz.position=[2,3.65];
% pard.drift_timepointsz.Optional=true;
% pard.drift_timepointsz.Width=.25;
% 
% pard.drift_pixreczt.object=struct('String','z binwidth nm','Style','text','Visible','off');
% pard.drift_pixreczt.position=[2,4];
% pard.drift_pixreczt.Optional=true;
% pard.drift_pixreczt.Width=.75;
% 
% pard.drift_pixrecz.object=struct('String','5','Style','edit','Visible','off');
% pard.drift_pixrecz.position=[2,4.65];
% pard.drift_pixrecz.isnumeric=1;
% pard.drift_pixrecz.object.TooltipString=sprintf('pixel size (nm) for reconstruction. \n Smaller for well defined peak. But slower \n Range: 10-25');
% pard.drift_pixrecz.Optional=true;
% pard.drift_pixrecz.Width=.25;
% 
% pard.drift_windowzt.object=struct('String','z fit window pix','Style','text','Visible','off');
% pard.drift_windowzt.position=[3,3];
% pard.drift_windowzt.Optional=true;
% pard.drift_windowzt.Width=.75;
% 
% pard.drift_windowz.object=struct('String','9','Style','edit','Visible','off');
% pard.drift_windowz.position=[3,3.65];
% pard.drift_windowz.isnumeric=1;
% pard.drift_windowz.object.TooltipString=sprintf('size of region for peakfinding (ellipt. Gaussian). \n should be small, but cover clear maximum. \n Range: 7-15');
% pard.drift_windowz.Optional=true;
% pard.drift_windowz.Width=.25;
% 
% 
% pard.zranget.object=struct('String','zrange nm','Style','text','Visible','off');
% pard.zranget.position=[4,3];
% pard.zranget.Optional=true;
% 
% pard.zrange.object=struct('String','-400 400','Style','edit','Visible','off');
% pard.zrange.position=[4,4];
% pard.zrange.Optional=true;
% pard.zrange.Width=0.9;
% 
% pard.slicewidtht.object=struct('String','slice width nm','Style','text','Visible','off');
% pard.slicewidtht.position=[5,3];
% pard.slicewidtht.Optional=true;
% 
% pard.slicewidth.object=struct('String','200','Style','edit','Visible','off');
% pard.slicewidth.position=[5,4];
% pard.slicewidth.Optional=true;
% pard.slicewidth.Width=0.90;
% 
% pard.smoothmode.object=struct('String',{{'smoothing cubic spline','linear'}},'Style','popupmenu');
% pard.smoothmode.position=[6,2];
% pard.smoothmode.Optional=true;
% pard.smoothmode.Width=2.;
% pard.smoothpar.object=struct('String','','Style','edit');
% pard.smoothpar.position=[6,4];
% pard.smoothpar.Optional=true;
% pard.smoothpar.Width=.5;
% pard.smoothpar.object.TooltipString=sprintf('Parameter for cubic splien interpolation. \n leave empty for automatic determination. \n 0.01 for little smoothing, 10 for strong smoothing.');
% 
% 
% pard.drift_reference.object=struct('String','reference is last frame','Style','checkbox');
% pard.drift_reference.position=[7,3];
% pard.drift_reference.Optional=true;
% pard.drift_reference.object.TooltipString=sprintf('If checked, drift at end of data set is set to zero. \n Useful for sequential acquisition, use this for first data set.');
% pard.drift_reference.Width=2;
% pard.drift_reference.Optional=true;
% 
% pard.drift_mirror2c.object=struct('String',{{'no mirror', '2 Channels, mirrored, vertical split', '2 Channels, mirrored, horizontal split'}},'Style','popupmenu');
% pard.drift_mirror2c.position=[7,1];
% pard.drift_mirror2c.Optional=true;
% pard.drift_mirror2c.object.TooltipString=sprintf('Vertical split: next to each other, horizontal split: below each other.');
% pard.drift_mirror2c.Width=2;
% pard.drift_mirror2c.Optional=true;
% 
% 
% % pard.drift_individual.object=struct('String','correct every file individually','Style','checkbox','Value',1);
% % pard.drift_individual.position=[8,1];
% % pard.drift_individual.Width=2;
% % pard.drift_individual.Optional=true;
% pard.drift_whatfiles.object=struct('String',{{'visible','all files'}},'Style','popupmenu','Value',1);
% pard.drift_whatfiles.position=[8,1];
% pard.drift_whatfiles.Width=1.5;
% pard.drift_whatfiles.Optional=true;
% 
% pard.drift_ask.object=struct('String','?','Style','checkbox','Value',0);
% pard.drift_ask.position=[8,2.6];
% pard.drift_ask.Width=.4;
% pard.drift_ask.Optional=true;
% 
% pard.save_dc.object=struct('String','Save driftcorrected SML','Style','checkbox','Value',1);
% pard.save_dc.position=[8,3];
% pard.save_dc.Width=2;
% pard.save_dc.Optional=true;

pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description={''};
end