classdef lineprofile<interfaces.DialogProcessor
    % LINEPROFILE Calculates profiles along a linear ROI and fits it with a model of choice
    methods
        function obj=lineprofile(varargin)        
            obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_layerson','linewidth_roi','znm_min','znm_max','sr_pixrec','layernames'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)           
            results=make_lineprofiles(obj.locData,p);
            out.clipboard=results;
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef
pard.setbinwidth.object=struct('String','set binwidth (nm) (otherwise: pixelsize):','Style','checkbox','Value',1);
pard.setbinwidth.position=[1,1];
pard.setbinwidth.Width=3;

pard.binwidth.object=struct('String','2','Style','edit');
pard.binwidth.position=[1,4];
pard.binwidth.TooltipString='Binwidth for profiles. If not checked, use pixel size of reconstruction';
pard.setbinwidth.TooltipString=pard.binwidth.TooltipString;

pard.text2.object=struct('String','fitmodel:','Style','text');
pard.text2.position=[2,1];

pard.fitmodel.object=struct('String','Gauss|Flat|Disk|Ring|Distance','Style','popupmenu');
pard.fitmodel.position=[2,2];

pard.restrictsigma.object=struct('String','sigma=<locp>','Style','checkbox');
pard.restrictsigma.position=[2,3];

pard.linelengthcheck.object=struct('String','set length (nm)','Style','checkbox');
pard.linelengthcheck.position=[3,1];

pard.linelength.object=struct('String','250','Style','edit');
pard.linelength.position=[3,2];
pard.linelength.TooltipString='This overrides the length of the ROI and uses a well-defined ROI. Useful for direct comparison.';
pard.linelengthcheck.TooltipString=pard.linelength.TooltipString;
pard.plugininfo.name='Line profiles';
pard.plugininfo.description=sprintf('Calculates profiles along a linear ROI and fits it with a model of choice. \n Flat: step function convolved with Gaussian (=Erf). \n Disk: Projection of a homogeneously filled disk, convolved with Gaussian. \n Ring: Projection of a ring, convolved with Gaussian. \n Distance: Two Gaussians in a distance d.');
pard.plugininfo.type='ProcessorPlugin';
end