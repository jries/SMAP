classdef BlinkingMoviePresentation<interfaces.DialogProcessor
%     BlinkingMoviePresentation renders a movie with the left side showing
%     the camera images and the right side showing sngle molecule
%     localizations slowly building up a superresolution image
    methods
        function obj=BlinkingMoviePresentation(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={'sr_size','sr_pos','sr_pixrec'};
            obj.showresults=true;
        end
        
        function out=run(obj,p)         
            out=[];
            file=obj.locData.files.file(p.dataselect.Value);
            locs=obj.locData.getloc({'xnm','ynm','frame','locprecnm'},'layer',1,'position','roi');
            makeBlinkMovie(locs,file,p);

        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.textx.object=struct('Style','text','String','filter settings from layer 1');
pard.textx.position=[1,3];
pard.textx.Width=2;

pard.dataselect.object=struct('Style','popupmenu','String','File');
pard.dataselect.position=[1,1];

pard.texta.object=struct('Style','text','String','number of frames');
pard.texta.position=[2,1];


pard.numberOfFrames.object=struct('Style','edit','String','1000');
pard.numberOfFrames.position=[2,2];

pard.textb.object=struct('Style','text','String','first frame');
pard.textb.position=[3,1];

pard.frame_min.object=struct('Style','edit','String','1');
pard.frame_min.position=[3,2];

pard.removedark.object=struct('Style','checkbox','String','remove dark frames','Value',1);
pard.removedark.position=[4,1];
pard.removedark.Width=2;
pard.removedark.TooltipString='removes frames without localiaztions';
%sum
pard.textc.object=struct('Style','text','String','frame rate fps');
pard.textc.position=[5,1];

pard.framerate.object=struct('Style','edit','String','20');
pard.framerate.position=[5,2];


pard.outputFormat.object=struct('Style','popupmenu','String',{{'MPEG-4','Uncompressed AVI','Motion JPEG 2000'}});
pard.outputFormat.position=[6,1];
pard.outputFormat.Width=2;
pard.syncParameters={{'filelist_short','dataselect',{'String'}}};
pard.plugininfo.name='BlinkingMoviePresentation';
pard.plugininfo.type='ProcessorPlugin';
end