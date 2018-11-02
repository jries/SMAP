classdef RegisterLocs<recgui.DialogProcessor
    properties
        isz=0;
        transformation=[];
    end
    methods
        function obj=RegisterLocs(varargin)        
                obj@recgui.DialogProcessor;
            if nargin>0
                obj.handle=varargin{1};
            end
        end
        
        function out=run(obj,p)
            p.isz=obj.isz;
            obj.transformation=   transform_locs(obj.locData,p);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function attachLocData(obj,locData)
            attachLocData@recgui.GuiProcessor(obj,locData);
%             addlistener(obj.locData,'synchronizeGui',@obj.synchronizeGui);
            addlistener(obj.locData,'loaded',@obj.updateGui);
        end
        function updateGui(obj,event,data)
            ff=obj.locData.files.file;
            str={};
            for k=1:length(ff)
                if 0%ff(k).istiff
                    str{end+1}=' ';
                else
                    str{end+1}=ff(k).name;
                end
            end
            obj.guihandles.targetselect.String=str;
            obj.guihandles.refselect.String=str;
            %PSF range, z rnage change
            zt={'zmin','tzmin','zmax'};
            pt={'PSFmax','tPSFmax'};
            hgui=obj.guihandles;
            
            if isfield(obj.locData.loc,'znm')
                reccontrolVisibility(hgui,zt,'on');
                reccontrolVisibility(hgui,pt,'off');
                obj.isz=1;
            else
                reccontrolVisibility(hgui,zt,'off');
                reccontrolVisibility(hgui,pt,'on');       
                obj.isz=0;
            end
            
        end
        function makeGui(obj)
            makeGui@recgui.DialogProcessor(obj);
            obj.guihandles.savebutton.Callback=@obj.savebutton;
        end
        function savebutton(obj,a,b)
            if isempty(obj.transformation)
                errordlg('first calculate a transformation')
            else
                fn=obj.guihandles.Tfile.String;
                [f,path]=uiputfile(fn,'Save last transformation as transformation file _T.mat');
                if f
                    obj.guihandles.Tfile.String=[path f];
                    transformation=obj.transformation; 
                    save([path f],'transformation');
                end 
            end
        end 
    end
    methods(Static)
        function info=info(obj)
            info.name='RegisterLocs';
            info.class=@RegisterLocs;
            info.tag='RegisterLocs';
%             obj.info=info;
        end

    end
end




function pard=guidef
pard.texta.object=struct('String','target','Style','text');
pard.texta.position=[1,2];
pard.textb.object=struct('String','reference','Style','text');
pard.textb.position=[1,1];
pard.targetselect.object=struct('Style','popupmenu','String','File');
pard.targetselect.position=[2,2];
pard.targetselect.object.TooltipString='file for target localizations';

pard.refselect.object=struct('Style','popupmenu','String','File');
pard.refselect.position=[2,1];
pard.refselect.object.TooltipString='file for reference localizations';

% pard.targetpart.object=struct('Style','popupmenu','String','all|left|right|top|bottom');
% pard.targetpart.position=[3,2];
% pard.targetpart.object.TooltipString='part of image to use as target';
pard.refpart.object=struct('Style','popupmenu','String','all|left-right|top-bottom');
pard.refpart.position=[3,1.5];
pard.refpart.object.TooltipString='part of image to use as reference';

pard.targetmirror.object=struct('Style','popupmenu','String','no mirror|left-right|up-down|both');
pard.targetmirror.position=[4,2];
pard.targetmirror.object.TooltipString='mirror target part';
% 
% pard.refmirror.object=struct('Style','popupmenu','String','no mirror|left-right|top-bottom');
% pard.refmirror.position=[4,1];
% pard.refmirror.object.TooltipString='mirror reference part';

pard.Tfile.object=struct('Style','edit','String','settings/temp/temp_T.mat');
pard.Tfile.position=[8,1];
pard.Tfile.Width=3;
pard.Tfile.object.TooltipString='default file for transformation matrix. You can select new file after transformation has been calculated.';

pard.savebutton.object=struct('Style','pushbutton','String','save');
pard.savebutton.position=[8,4];
pard.savebutton.object.TooltipString='Save the newly calculated transformation matrix.';


pard.transform.object=struct('Style','popupmenu','String','projective|affine|similarity|polynomial|lwm|pwl');
pard.transform.position=[6,1];
pard.transform.object.TooltipString='select one of Matlabs transformations. Not all might work.';

pard.transformparam.object=struct('Style','edit','String','3');
pard.transformparam.position=[6,2];
pard.transformparam.object.TooltipString='Parameter for lwm and polynomial';

pard.tlp.object=struct('Style','text','String','maxlocp');
pard.tlp.position=[1,3];

pard.maxlocprec.object=struct('Style','edit','String','15');
pard.maxlocprec.position=[1,4];
pard.maxlocprec.object.TooltipString='maximum localization precision (nm) of locs used to calculate transformation';
pard.tlp.object.TooltipString=pard.maxlocprec.object.TooltipString;

pard.tzmin.object=struct('Style','text','String','zmin, zmax');
pard.tzmin.position=[2,3];

pard.zmin.object=struct('Style','edit','String','-200');
pard.zmin.position=[2,4];
pard.zmin.Width=0.5;

pard.zmin.object.TooltipString='zrange (nm) used to calculate transformation. ';
pard.tzmin.object.TooltipString=pard.zmin.object.TooltipString;
% pard.tzmax.object=struct('Style','text','String','zmax');
% pard.tzmax.position=[3,4.5];

pard.zmax.object=struct('Style','edit','String','200');
pard.zmax.position=[2,4.5];
pard.zmax.Width=0.5;
pard.zmax.object.TooltipString=pard.zmin.object.TooltipString;

pard.tPSFmax.object=struct('Style','text','String','PSF max');
pard.tPSFmax.position=[2,3];

pard.PSFmax.object=struct('Style','edit','String','160');
pard.PSFmax.position=[2,4];
pard.PSFmax.object.TooltipString='maximum size of PSF (nm) of locs used to calculate transformation';
pard.tPSFmax.object.TooltipString=pard.PSFmax.object.TooltipString;

pard.pixst.object=struct('Style','text','String','pixelsize nm');
pard.pixst.position=[3,3];


pard.pixrec.object=struct('Style','edit','String','10');
pard.pixrec.position=[3,4];
pard.pixrec.object.TooltipString='image correlation: size of pixel for reconstruction (nm). Typical: 10-50';
pard.pixst.object.TooltipString=pard.pixrec.object.TooltipString;


pard.maxshiftt.object=struct('Style','text','String','maxshift nm');
pard.maxshiftt.position=[4,3];
pard.maxshift.object=struct('Style','edit','String','5000');
pard.maxshift.position=[4,4];

pard.maxshiftt.object.TooltipString=sprintf('correlation: maximum shift between halves (nm). Typical 1000-10000. \n Increase if no good maximum is found in shiftcorr. \n Decrease if wrong maximum is found');
pard.maxshift.object.TooltipString=pard.maxshiftt.object.TooltipString;

pard.maxlocsusedt.object=struct('Style','text','String','maxlocs used');
pard.maxlocsusedt.position=[5,3];
pard.maxlocsused.object=struct('Style','edit','String','10000');
pard.maxlocsused.position=[5,4];
pard.maxlocsused.object.TooltipString=sprintf('number of localizations used to calculate correlation. Typical >10000. \n reduce if too slow, increase if too imprecise.');
pard.maxlocsusedt.object.TooltipString=pard.maxlocsused.object.TooltipString;

pard.maxshiftregistert.object=struct('Style','text','String','max shift register (nm)');
pard.maxshiftregistert.position=[6,3];
pard.maxshiftregister.object=struct('Style','edit','String','250');
pard.maxshiftregister.position=[6,4];
pard.maxshiftregister.object.TooltipString=sprintf('Max shift (nm) between corresponding loclaizations in shifted images. \n Could be due to aberrations or limited loclaization precision \n Typical 200-500');
pard.maxshiftregistert.object.TooltipString=pard.maxshiftregister.object.TooltipString;


pard.midpcheck.object=struct('Style','checkbox','String','set midpoint (pix) to');
pard.midpcheck.position=[7,3];
pard.midp.object=struct('Style','edit','String','256');
pard.midp.position=[7,4];
pard.midp.Width=0.5;

pard.midpcheck.object.TooltipString=sprintf('midpoint in camera pixels. Typically 256. If not checked, the middle of the camera frame will be used. \n Use this to bring the correct correlation maximum into the middle in order to reject wrong maxima.');
pard.midp.object.TooltipString=pard.midpcheck.object.TooltipString;


pard.useroi.object=struct('Style','checkbox','String','use ROI');
pard.useroi.position=[7,2];
pard.useroi.object.TooltipString='if checked, use only localizations in ROI and layer 1, make sure to select both halves';
end