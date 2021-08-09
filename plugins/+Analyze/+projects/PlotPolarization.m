classdef PlotPolarization<interfaces.DialogProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a
    % model of choice. Flat: step function convolved with Gaussian
    % (=Erf). Disk: Projection of a homogeneously filled disk, convolved
    % with Gaussian. Ring: Projection of a ring, convolved with
    % Gaussian. Distance: Two Gaussians in a distance d.
    methods
        function obj=PlotPolarization(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
%             obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)           
            locs=obj.locData.getloc({'xnm','ynm','P','inMask'},'layer',1,'Position','roi');
            img=obj.getPar('sr_image');
            ax=obj.initaxis('out');
            hold(ax,'off')
             imagesc(ax,img.rangex,img.rangey,img.image)
            hold(ax,'on')
            
            Ph=min(max(locs.P,-1),1);
            switch p.orient.selection
                case '\'
                    orient=acos(Ph)/2+0*pi/2;
                case '/'
                    orient=acos(-Ph)/2+pi/2;
            end
            d=p.lengthline;
            u=d*cos(orient);
            v=d*sin(orient);   
            
            if p.numberlines<length(orient)
                indg=ceil(rand(p.numberlines,1)*length(orient));
                valid=false(length(orient),1);
                valid(indg)=true;
            else
                valid=true(length(orient),1);
            end
            
            if p.inmask
                valid=valid & locs.inMask;
            end
            colors=hsv(256);
            cind=floor((Ph+1)/2*255)+1;
            for k=1:length(locs.xnm)
                if ~valid(k)
                    continue
                end
                x=locs.xnm(k)+u(k)*[-1 1];
                y=locs.ynm(k)+v(k)*[-1 1];
                line(ax,x/1000,y/1000,'Color',colors(cind(k),:))
                hold(ax,'on')
            end
            axis(ax,'equal')
                out=[];
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end




function pard=guidef(obj)
% 
%  p(1).value=0; p(1).on={}; p(1).off={'binwidth'};
% p(2).value=1; p(2).on={'binwidth'}; p(2).off={};
% pard.setbinwidth.object=struct('String','set binwidth (nm) (otherwise: pixelsize):','Style','checkbox','Value',1,'Callback',{{@obj.switchvisible,p}});
% pard.setbinwidth.position=[1,1];
% pard.setbinwidth.Width=3;

pard.lengthlinet.object=struct('String','length line (nm)','Style','text');
pard.lengthlinet.position=[1,1];
pard.lengthlinet.Width=1;
pard.lengthline.object=struct('String','25','Style','edit');
pard.lengthline.position=[1,2];
pard.lengthline.Width=.5;


pard.orientt.object=struct('String','orientation','Style','text');
pard.orientt.position=[2,1];
pard.orientt.Width=1;
pard.orient.object=struct('String',{{'/','\'}},'Style','popupmenu');
pard.orient.position=[2,2];
pard.orient.Width=.5;

pard.numberlinest.object=struct('String','number of lines','Style','text');
pard.numberlinest.position=[1,3];
pard.numberlinest.Width=1;
pard.numberlines.object=struct('String','inf','Style','edit');
pard.numberlines.position=[1,4];
pard.numberlines.Width=.5;

pard.inmask.object=struct('String','use only localizations in mask','Style','checkbox');
pard.inmask.position=[3,1];
pard.inmask.Width=2;
% pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
% pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;
% 
% pard.text2.object=struct('String','fitmodel:','Style','text');
% pard.text2.position=[2,1];
% 
% pard.fitmodel.object=struct('String','Gauss|Flat|Disk|Ring|Distance','Style','popupmenu');
% pard.fitmodel.position=[2,2];
% 
% pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
% pard.restrictsigma.position=[2,3];
% 
% 
% p(1).value=0; p(1).on={}; p(1).off={'linelength'};
% p(2).value=1; p(2).on={'linelength'}; p(2).off={};
% pard.linelengthcheck.object=struct('String','set length (nm)','Style','checkbox','Callback',{{@obj.switchvisible,p}});
% pard.linelengthcheck.position=[3,1];
% 
% pard.linelength.object=struct('String','250','Style','edit');
% pard.linelength.position=[3,2];
% pard.linelength.TooltipString='This overrides the length of the ROI and uses a well-defined ROI. Useful for direct comparison.';
% pard.linelengthcheck.TooltipString=pard.linelength.TooltipString;
% pard.plugininfo.name='Line profiles';
pard.plugininfo.description='';
pard.plugininfo.type='ProcessorPlugin';
end